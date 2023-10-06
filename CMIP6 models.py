"""Validate historical CMIP6 models

1. Collect CMIP6 model files (.nc)
2. Collect GSOD data

Install necessary packages in terminal:

    conda install -f requirements.txt

"""
import datetime  # used to convert time from gregorian to POSIXct format
import os
import concurrent.futures
from functools import partial

import netCDF4  # to read .nc files
import numpy as np
import pandas as pd
from scipy.stats import pearsonr


def main():
    cmip6_data_folder = "/ouce-home/data/model/cmip6/CMIP6/"

    # Upload the historical CMIP6 model netCDF files
    # Read netCDF precipitation files (r1i1p1f1 Variant level used for 4 models
    # (CanESM5, IPSL-CM 6A-LR, MIROC6 and MRI-ESM2-0)
    historical_filepaths = {
        "miroc6": (
            "CMIP/MIROC/MIROC6/historical/r1i1p1f1/Amon/pr/gn/latest/"
            "pr_Amon_MIROC6_historical_r1i1p1f1_gn_185001-201412.nc"
        ),
        "awi": (
            "CMIP/AWI/AWI-CM-1-1-MR/historical/r1i1p1f1/Amon/pr/gn/latest/"
            "pr_Amon_AWI-CM-1-1-MR_historical_r1i1p1f1_gn_185001-201412.nc"
        ),
        "cesm2": (
            "CMIP/NCAR/CESM2-WACCM/historical/r1i1p1f1/Amon/pr/gn/latest/"
            "pr_Amon_CESM2-WACCM_historical_r1i1p1f1_gn_185001-201412.nc"
        ),
        "fgoals": (
            "CMIP/CAS/FGOALS-g3/historical/r1i1p1f1/Amon/pr/gn/latest/"
            "pr_Amon_FGOALS-g3_historical_r1i1p1f1_gn_185001-201612.nc"
        ),
    }

    # Load future CMIP6 model data
    future_filepaths = {}
    for ssp in ("ssp126", "ssp245", "ssp585"):
        future_filepaths[ssp, "miroc6"] = (
            f"ScenarioMIP/MIROC/MIROC6/{ssp}/r1i1p1f1/Amon/pr/gn/latest/"
            f"pr_Amon_MIROC6_{ssp}_r1i1p1f1_gn_201501-210012.nc"
        )
        future_filepaths[ssp, "awi"] = (
            f"ScenarioMIP/AWI/AWI-CM-1-1-MR/{ssp}/r1i1p1f1/Amon/pr/gn/latest/"
            f"pr_Amon_AWI-CM-1-1-MR_{ssp}_r1i1p1f1_gn_201501-210012.nc"
        )
        future_filepaths[ssp, "cesm2"] = (
            f"ScenarioMIP/NCAR/CESM2-WACCM/{ssp}/r1i1p1f1/Amon/pr/gn/latest/"
            f"pr_Amon_CESM2-WACCM_{ssp}_r1i1p1f1_gn_201501-229912.nc"
        )
        future_filepaths[ssp, "fgoals"] = (
            f"ScenarioMIP/CAS/FGOALS-g3/{ssp}/r1i1p1f1/Amon/pr/gn/latest/"
            f"pr_Amon_FGOALS-g3_{ssp}_r1i1p1f1_gn_201501-210012.nc"
        )

    # For testing purposes use GSOD data for Uganda already collected
    # consider using CRU data
    observed_data = pd.read_csv(
        "/Users/pamelaacheng/Library/CloudStorage/OneDrive-Nexus365/DPhil/2022/Landslides/Data/Precipitation/gsod_uganda.csv"
    )

    process = partial(
        process_model_data,
        cmip6_data_folder=cmip6_data_folder,
        observed_data=observed_data,
    )

    # Bias correction for each CMIP6 model
    with concurrent.futures.ProcessPoolExecutor() as executor:
        model, correlation_coefficient = executor.map(
            process, historical_filepaths.items()
        )
        print(f"Processed {model=}")
        print(f"Correlation coefficient = {correlation_coefficient[0]}")


def process_model_data(model, filepath, cmip6_data_folder=None, observed_data=None):
    model_data = read_dataset(os.path.join(cmip6_data_folder, filepath))
    corrected_model_data = bias_correct(observed_data, model_data)
    # Calculate the correlation between observed and model data
    correlation_coefficient = pearsonr(observed_data, corrected_model_data)
    return model, correlation_coefficient


def bias_correct(observed_data, model_data):
    # Calculate the bias correction factor
    bias_correction_factor = np.nanmean(observed_data) / np.nanmean(model_data)
    # Apply bias correction to the model data
    return model_data * bias_correction_factor


def read_dataset(filename):
    # Open the dataset
    nc_p = netCDF4.Dataset(filename, "r")

    # Get the list of variable names
    variable_names = nc_p.variables.keys()

    # Print the variable names
    print(variable_names)

    # CHECK THE FORMAT OF THE TIME VARIABLE (nc files usually have time in UTC,
    # Gregorian, etc format which needs to be changed to POSIXct object)

    # Access the time variable
    time_variable = nc_p.variables["time"]

    print("Time variable attributes:")
    print("Units:", time_variable.units)
    print("Calendar:", time_variable.calendar)
    # Time variable attributes:
    # Units: days since 0001-01-01 00:00:00
    # Calendar: noleap

    # CONVERT TIME FROM PROLEPTIC GREGORIAN TO POSIXct FORMAT
    # Access the time variable
    time_variable = nc_p.variables["time"]
    # time_datetime = nc_p.variables['time']
    # Get the time values
    time_values = time_variable[:]
    # Convert time values to datetime objects

    base_date = datetime.datetime(1, 1, 1)
    time_datetime = [base_date + datetime.timedelta(days=int(t)) for t in time_values]

    # Get latitude and longitude values
    lat_values = nc_p.variables["lat"][:]
    lon_values = nc_p.variables["lon"][:]

    # Get precipitation variable
    precipitation_variable = nc_p.variables["pr"]
    precipitation_data = []

    # Write data rows
    for i, t in enumerate(time_datetime):
        year, month, day = t.year, t.month, t.day

        for lat in lat_values:
            for lon in lon_values:
                lat_index = (lat_values == lat).nonzero()[0][0]
                lon_index = (lon_values == lon).nonzero()[0][0]

                precipitation_value = precipitation_variable[i, lat_index, lon_index]
                Pr_mmd = (
                    precipitation_value * 86400
                )  # convert precipitation from kg m-2 s-1 to mm/day

                row_data = {
                    "lat": lat,
                    "lon": lon,
                    "t": t,
                    "year": year,
                    "month": month,
                    "day": day,
                    "precipitation_value": precipitation_value,
                    "Pr_mmd": Pr_mmd,
                }
                precipitation_data.append(row_data)

    # Close the NetCDF file
    nc_p.close()

    # Prcp is in kg m-2 s-1

    # Print the precipitation dataframe
    awi1 = pd.DataFrame(precipitation_data)
    print(awi1.head(-6))

    # Check the years of the future precipitation data (data validation check)
    prec_years = awi1["year"].unique()
    print(prec_years)
    return precipitation_data

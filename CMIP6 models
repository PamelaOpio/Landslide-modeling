# Validate historical CMIP6 models
#1. Collect CMIP6 model files (.nc)
#2. Collect GSOD data


#Install necessary packages in terminal 
#conda install requests
import pandas as pd
#conda install -c conda-forge sklearn-pandas
#conda install -c conda-forge matplotlib
#conda install sklearn.neighbors
#conda install scikit-learn
#conda install -c conda-forge scipy #for nearest neighbor algorithms    
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt #to plot the maps
import cartopy.crs as crs
import netCDF4 #to convert .nc files to .csv files
import csv
import ee
import requests
import io

#Upload the historical CMIP6 model netCDF files
#Read netCDF precipitation files (r1i1p1f1 Variant level used for 4 models (CanESM5, IPSL-CM 6A-LR, MIROC6 and MRI-ESM2-0)
#miroc6 = "/Users/pamelaacheng/Library/CloudStorage/OneDrive-Nexus365/DPhil/2022/HDM4/CMIP6 data/Historical /pr_Amon_MIROC6_historical_r10i1p1f1_gn_195001-201412.nc"
awi = "/Users/pamelaacheng/Library/CloudStorage/OneDrive-Nexus365/DPhil/2022/HDM4/CMIP6 data/Historical /pr_Amon_AWI-CM-1-1-MR_historical_r1i1p1f1_gn_200801-200812.nc"
#cesm2 = "/Users/pamelaacheng/Library/CloudStorage/OneDrive-Nexus365/DPhil/2022/HDM4/CMIP6 data/Historical /pr_day_CESM2-WACCM_historical_r1i1p1f1_gn_20100101-20150101.nc"
#fgoals = ""

#Read the variable names
#Open the dataset
nc_p = netCDF4.Dataset(awi, "r")

# Get the list of variable names
variable_names = nc_p.variables.keys()

# Print the variable names
print(variable_names)
#dict_keys(['time', 'time_bnds', 'lat', 'lat_bnds', 'lon', 'lon_bnds', 'pr'])

#CHECK THE FORMAT OF THE TIME VARIABLE (nc files usually have time in UTC, Gregorian, etc format which needs to be changed to POSIXct object)
# Access the time variable
time_variable = nc_p.variables['time']

print("Time variable attributes:")
print("Units:", time_variable.units)
print("Calendar:", time_variable.calendar)
#Time variable attributes:
#Units: days since 0001-01-01 00:00:00
#Calendar: noleap

#Convert the PRCP nc file to csv
# Open the NetCDF file
csv_p = "netcdf_p.csv" #define the output file

#CONVERT TIME FROM PROLEPTIC GREGORIAN TO POSIXct FORMAT
# Access the time variable
time_variable = nc_p.variables['time']
#time_datetime = nc_p.variables['time']
# Get the time values
time_values = time_variable[:]
# Convert time values to datetime objects
import datetime #used to convert time from gregorian to POSIXct format
base_date = datetime.datetime(1, 1, 1)
time_datetime = [base_date + datetime.timedelta(days=int(t)) for t in time_values]
# Extract year, month, and day from the converted time
years = [t.year for t in time_datetime]
months = [t.month for t in time_datetime]
days = [t.day for t in time_datetime]

# Get latitude and longitude values
lat_values = nc_p.variables['lat'][:]
lon_values = nc_p.variables['lon'][:]

# Get precipitation variable
precipitation_variable = nc_p.variables['pr']

# Prepare CSV file
with open(csv_p, "w", newline="") as csvfile:
    csv_writer = csv.writer(csvfile)

    # Write header
    header = ["Lat", "Lon", "Time", "Year", "Month", "Day", "Precipitation", "Pr_mmd"]
    csv_writer.writerow(header)

    # Write data rows
    for i, t in enumerate(time_datetime):
        year, month, day = t.year, t.month, t.day
        
        for lat in lat_values:
            for lon in lon_values:
                lat_index = (lat_values == lat).nonzero()[0][0]
                lon_index = (lon_values == lon).nonzero()[0][0]
                
                precipitation_value = precipitation_variable[i, lat_index, lon_index]
                Pr_mmd = precipitation_value * 86400 #convert precipitation from kg m-2 s-1 to mm/day
                
                row_data = [lat, lon, t, year, month, day, precipitation_value, Pr_mmd]
                csv_writer.writerow(row_data)

# Close the NetCDF file
nc_p.close()

#Prcp is in kg m-2 s-1

#Print the precipitation dataframe
awi1 = pd.read_csv('netcdf_p.csv')
print(awi1.head(-6))

#Check the years of the future precipitation data (data validation check)
prec_years = awi['Year'].unique()
print(prec_years)

#For testing purposes use GSOD data for Uganda already collected
gsod_Ug = pd.read_csv("/Users/pamelaacheng/Library/CloudStorage/OneDrive-Nexus365/DPhil/2022/Landslides/Data/Precipitation/gsod_uganda.csv")

#Bias correction for each CMIP6 model
import numpy as np

# Load observed and model data (repeat this process for each of the models)
observed_data = gsod_Ug #consider using CRU data
model_data1 = awi1
#model_data2 = miroc6
#model_data3 = cesm2
#model_data4 = fgoals

# Calculate the bias correction factor
bias_correction_factor1 = np.nanmean(observed_data) / np.nanmean(model_data1)
#bias_correction_factor2 = np.nanmean(observed_data) / np.nanmean(model_data2)
#bias_correction_factor3 = np.nanmean(observed_data) / np.nanmean(model_data3)
#bias_correction_factor4 = np.nanmean(observed_data) / np.nanmean(model_data4)

# Apply bias correction to the model data
corrected_model_data1 = model_data1 * bias_correction_factor1
#corrected_model_data2 = model_data1 * bias_correction_factor2
#corrected_model_data3 = model_data1 * bias_correction_factor3
#corrected_model_data4 = model_data1 * bias_correction_factor4

#Validation
from scipy.stats import pearsonr


# Calculate the correlation between observed and model data (repeat process for each of the models)
correlation_coefficient1 = pearsonr(observed_data, corrected_model_data1)
#correlation_coefficient2 = pearsonr(observed_data, corrected_model_data2)
#correlation_coefficient3 = pearsonr(observed_data, corrected_model_data3)
#correlation_coefficient4 = pearsonr(observed_data, corrected_model_data4)

# Print the correlation coefficient
print(f"Correlation Coefficient: {correlation_coefficient1[0]}")


#Modelling future rainfall
# Load future CMIP6 model data 
future_model_data_miroc6 = "/Users/pamelaacheng/Library/CloudStorage/OneDrive-Nexus365/DPhil/2022/HDM4/CMIP6 data/MIROC6/SSP126/pr_day_MIROC6_ssp126_r1i1p1f1_gn_20250101-20341231.nc" 
#future_model_data_fgoals = # Load future data for FGOALS model
#future_model_data_cesm2 = # Load future data for CESM2 model
#future_model_data_awi = # Load future data for AWI model

# Create an ensemble by averaging the projections
ensemble_future_data = (future_model_data_miroc6 + future_model_data_fgoals + future_model_data_cesm2 + future_model_data_awi) / 4.0

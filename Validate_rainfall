#Install necessary libraries
import pandas as pd
#conda install -c conda-forge sklearn-pandas
#conda install -c conda-forge matplotlib
#conda install sklearn.neighbors
#conda install scikit-learn
#conda install -c conda-forge scipy #for nearest neighbor algorithms    
import numpy as np
import matplotlib.pyplot as plt 
import ee

#Access the rainfall data from google earth engine
#Authenticate and initialize GEE
ee.Authenticate()
ee.Initialize()

# Load Reanalysis precipitation data
#IMERG rain in mm/hr
imerg = ee.ImageCollection('NASA/GPM_L3/IMERG_V06') \
    #.filterBounds(uganda) \
    #.filterDate(start_date, end_date) \
    #.select('precipitationCal')
#CHIRPS rain in mm/day
chirps = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY')\
  #.filterBounds(uganda)\
  #.filterDate(start_date, end_date) \
  #.select('precipitation');
#ERA5 #rain in m/day
era5 = ee.ImageCollection("ECMWF/ERA5/DAILY")\
  #.filterBounds(uganda)\
  #.filterDate(start_date, end_date) \
  #.select('total_precipitation');
#GSMaP #rain in mm/hr
#gsmap = ee.ImageCollection('JAXA/GPM_L3/GSMaP/v6/reanalysis')\
  #.filterBounds(uganda)\
  #.filterDate(start_date, end_date) \
  #.select('hourlyPrecipRateGC');
# Create the image collections for the two date ranges
gsmap_reanalysis = ee.ImageCollection("JAXA/GPM_L3/GSMaP/v6/reanalysis") \
    .filterDate('2000-01-01', '2014-12-31') \
    .select('hourlyPrecipRateGC')

gsmap_operational = ee.ImageCollection("JAXA/GPM_L3/GSMaP/v6/operational") \
    .filterDate('2014-01-01', '2023-12-31') \
    .select('hourlyPrecipRateGC')
# Merge the two image collections
gsmap = gsmap_reanalysis.merge(gsmap_operational)

##Collect GSOD data for Uganda
import requests
import io
# Define the URL to retrieve GSOD data for multiple stations
station_codes = ['63602099999', '63702099999', '63705099999', '63630099999', '63726099999', '63630499999', '63654099999', '63658099999', '63674099999', '63680099999', '63682099999', '63684099999']  # station codes for Uganda
station_codes_str = ','.join(station_codes)  # Join the station codes with a comma

url = f"https://www.ncei.noaa.gov/access/services/data/v1?dataset=global-summary-of-the-day&stations={station_codes_str}&startDate=2010-01-01&endDate=2022-12-31&dataTypes=PRCP,LATITUDE,LONGITUDE,NAME&format=csv"

# Send a GET request to retrieve the data
response = requests.get(url)

# Load the response content into a DataFrame
gsod_Ug = pd.read_csv(io.StringIO(response.text))

# Save the filtered data to a CSV file
gsod_Ug.to_csv("/Users/pamelaacheng/Library/CloudStorage/OneDrive-Nexus365/DPhil/2022/Landslides/Data/Precipitation/gsod_uganda.csv", index=False)

#Run-time for entire data set taking too long so split data by year to iterate over each gsod location and retrieve reanalysis data 
#Split data into different years to reduce processing time
# Preprocess the GSOD data 
#gsod_Ug = gsod_Ug.dropna(subset=['PRCP'])  # Remove rows with missing precipitation values

# Convert 'DATE' column to datetime type
gsod_Ug['DATE'] = pd.to_datetime(gsod_Ug['DATE'])

# Filter gsod_Ug to include only dates range of interest
gsod_Ug = gsod_Ug.loc[(gsod_Ug['DATE'].dt.year >= 2018) & (gsod_Ug['DATE'].dt.year <= 2018)]

# Reset the index of the DataFrame
gsod_Ug = gsod_Ug.reset_index(drop=True)

gsod_data = gsod_Ug

# Define your CHIRPS image collection 1981-2023
def collect_chirps_data(date, latitude, longitude):
    # Convert date to Earth Engine format
    ee_date = ee.Date(date)

    # Filter Chirps collection based on date and location
    chirps_collection = chirps.filterDate(ee_date, ee_date.advance(1, 'day'))

    # Check if Chirps collection has any images
    if chirps_collection.size().getInfo() > 0:
    # Get the first image from the collection
        chirps_image = ee.Image(chirps_collection.first())

    # Create a point geometry for the location
        point = ee.Geometry.Point(longitude, latitude)

    # Reduce the image to the point location
        value = chirps_image.reduceRegion(ee.Reducer.first(), point, 1000)

    # Get the precipitation value
        precipitation_value = value.get('precipitation')

    # Check if the precipitation value is available
        if precipitation_value.getInfo() is not None:
        # Multiply the value by 24 to convert from m/day to mm/day
            value_mmd = ee.Number(precipitation_value)

        # Print the chirps rainfall value
            return value_mmd.getInfo()
        else:
        # Print 0 when no data is available
            return 0
    else:
    # Print 0 when no data is available
        return 0

# Define your ERA5 image collection from 1979-2020
def collect_era5_data(date, latitude, longitude):
    # Convert date to Earth Engine format
    ee_date = ee.Date(date)
    
    # Filter ERA5 collection based on date and location
    era5_collection = era5.filterDate(ee_date, ee_date.advance(1, 'day'))

    # Check if era5 collection has any images
    if era5_collection.size().getInfo() > 0:
    # Get the first image from the collection
        era5_image = ee.Image(era5_collection.first())

    # Create a point geometry for the location
        point = ee.Geometry.Point(longitude, latitude)

    # Reduce the image to the point location
        value = era5_image.reduceRegion(ee.Reducer.first(), point, 1000)

    # Get the precipitation value
        precipitation_value = value.get('total_precipitation')

    # Check if the precipitation value is available
        if precipitation_value.getInfo() is not None:
        # Multiply the value by 24 to convert from m/day to mm/day
            value_mmd = ee.Number(precipitation_value).multiply(1000)

        # Print the ERA5 rainfall value
            #return value_mmd.getInfo()
            return precipitation_value.getInfo()
        else:
        # Print 0 when no data is available
            return 0
    else:
    # Print 0 when no data is available
        return 0

# Define your GSMaP image collection #2000-2014 and 2014-2023
def collect_GSMaP_data(date, latitude, longitude):
    # Convert date to Earth Engine format
    ee_date = ee.Date(date)

    # Filter GSMaP collection based on date and location
    gsmap_collection = gsmap.filterDate(ee_date, ee_date.advance(1, 'day'))

    # Check if GSMaP collection has any images
    if gsmap_collection.size().getInfo() > 0:
    # Get the first image from the collection
        gsmap_image = ee.Image(gsmap_collection.first())

    # Create a point geometry for the location
        point = ee.Geometry.Point(longitude, latitude)

    # Reduce the image to the point location
        value = gsmap_image.reduceRegion(ee.Reducer.first(), point, 1000)

    # Get the precipitation value
        precipitation_value = value.get('hourlyPrecipRateGC')

    # Check if the precipitation value is available
        if precipitation_value.getInfo() is not None:
        # Multiply the value by 24 to convert from mm/hr to mm/day
            #value_mmd = ee.Number(precipitation_value).multiply(24)

        # Print the GSMaP rainfall value
            #return value_mmd.getInfo()
            return precipitation_value.getInfo()
        else:
        # Print 0 when no data is available
            return 0
    else:
    # Print 0 when no data is available
        return 0

#Function to collect imerg images
def collect_imerg_data(date, latitude, longitude):
    # Convert date to Earth Engine format
    ee_date = ee.Date(date)
    
        # Filter imerg collection based on date and location
    imerg_collection = imerg.filterDate(ee_date, ee_date.advance(1, 'day'))

    # Check if imerg collection has any images
    if imerg_collection.size().getInfo() > 0:
    # Get the first image from the collection
        imerg_image = ee.Image(imerg_collection.first())

    # Create a point geometry for the location
        point = ee.Geometry.Point(longitude, latitude)

    # Reduce the image to the point location
        value = imerg_image.reduceRegion(ee.Reducer.first(), point, 1000)

    # Get the precipitation value
        precipitation_value = value.get('precipitationCal')

    # Check if the precipitation value is available
        if precipitation_value.getInfo() is not None:
        # Multiply the value by 24 to convert from mm/hr to mm/day
            #value_mmd = ee.Number(precipitation_value).multiply(24)

        # Print the imerg rainfall value
            return precipitation_value.getInfo()
            #return value_mmd.getInfo()
        else:
        # Print 0 when no data is available
            return 0
    else:
    # Print 0 when no data is available
        return 0

# Create an empty DataFrame to store the validation data
validation_data = pd.DataFrame(columns=['latitude', 'longitude', 'name', 'date', 'gsod_rain', 'chirps_rain', 'era5_rain', 'imerg_rain', 'GSMaP_rain'])

print(validation_data)

# Iterate GSOD over each row in the data
for index, row in gsod_data.iterrows():
    latitude = row['LATITUDE']
    longitude = row['LONGITUDE']
    date = row['DATE']
    name = row['NAME']
    gsod_rainfall = row['PRCP']

    # Collect CHIRPS data
    chirps_rainfall = collect_chirps_data(date, latitude, longitude)

    # Collect ERA5 data
    era5_rainfall = collect_era5_data(date, latitude, longitude)

    # Collect IMERG data
    imerg_rainfall = collect_imerg_data(date, latitude, longitude)

    # Collect GSMaP data
    GSMaP_rainfall = collect_GSMaP_data(date, latitude, longitude)

    # Add the validation data to the DataFrame
    row_data = {
        'latitude': latitude,
        'longitude': longitude,
        'date': date,
        'name': name,
        'gsod_rain': gsod_rainfall,
        'chirps_rain': chirps_rainfall,
        'era5_rain': era5_rainfall,
        'imerg_rain': imerg_rainfall,
        'GSMaP_rain': GSMaP_rainfall
    }
    validation_data = pd.concat([validation_data, pd.DataFrame([row_data])])

# Convert rainfall data to mm/day for all collections
# Multiply era5  by 1000 to convert m/day to mm/day
validation_data['era5_rain'] = validation_data['era5_rain'] * 1000

# Multiply imerg and GSMaP by 24 to convert from mm/hr to mm/day
validation_data['imerg_rain'] = validation_data['imerg_rain'] * 24
validation_data['GSMaP_rain'] = validation_data['GSMaP_rain'] * 24

# Convert PRCP data FROM INCHES/DAY TO MM/DAY
validation_data['gsod_rain'] = validation_data['gsod_rain'] * 25.4

print(validation_data)

#Combine teh different csv files into one for validation
#Collect and combine the validation csv files into one
import glob

# Get a list of all CSV files in the directory
val_files = glob.glob('/Users/pamelaacheng/Library/CloudStorage/OneDrive-Nexus365/DPhil/2022/Landslides/Data/Precipitation/Validation files/*.csv')

# Initialize an empty list to store individual DataFrames
dataframes = []

# Read each CSV file into a DataFrame and append it to the list
for csv_file in val_files:
    df = pd.read_csv(csv_file)
    dataframes.append(df)

# Concatenate all DataFrames into a single DataFrame
combined_val = pd.concat(dataframes, ignore_index=True)

print(combined_val.head(2460))

#Correlation analysis
import numpy as np
# Exclude rows with missing or non-numeric values
#valid_data = validation_data.dropna(subset=['chirps_rain', 'gsod_rain', 'era5_rain', 'imerg_rain', 'GSMaP_rain'], how='any')
gsod_rain = combined_val['gsod_rain'].astype(float)
chirps_rain = combined_val['chirps_rain'].astype(float)
era5_rain = combined_val['era5_rain'].astype(float)
imerg_rain = combined_val['imerg_rain'].astype(float)
GSMaP_rain = combined_val['GSMaP_rain'].astype(float)

# Calculate correlation coefficient
corr_chirps = np.corrcoef(chirps_rain, gsod_rain)[0, 1]
corr_era5 = np.corrcoef(era5_rain, gsod_rain)[0, 1]
corr_imerg = np.corrcoef(imerg_rain, gsod_rain)[0, 1]
corr_GSMaP = np.corrcoef(GSMaP_rain, gsod_rain)[0, 1]
print('Correlation CHIRPS:', corr_chirps)
print('Correlation ERA5:', corr_era5)
print('Correlation IMERG:', corr_imerg)
print('Correlation GSMaP:', corr_GSMaP)

#Root Mean Square Error analysis
from sklearn.metrics import mean_squared_error
import numpy as np
import math

#Comparative contour map plot of CHIRPS and GSOD
#Necessary libraries
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

#Time series 

# Convert the 'date' column to datetime format
combined_val['date'] = pd.to_datetime(combined_val['date'])

# Group the validation data by date and calculate the total/mean precipitation for CHIRPS and GSOD

gsod_tt= combined_val.groupby('date')['gsod_rain'].mean()
chirps_tt = combined_val.groupby('date')['chirps_rain'].mean()
era5_tt = combined_val.groupby('date')['era5_rain'].mean()
imerg_tt = combined_val.groupby('date')['imerg_rain'].mean()
GSMaP_tt = combined_val.groupby('date')['GSMaP_rain'].mean()

# Calculate the coefficient of determination (R-squared)
rsqu_chirps = corr_chirps ** 2
rsqu_era5 = corr_era5 ** 2
rsqu_imerg = corr_imerg ** 2
rsqu_GSMaP = corr_GSMaP ** 2

# Convert R-squared to a percentage with one decimal place
rsq_chirps_per = round(rsqu_chirps * 100, 1)
rsq_era5_per = round(rsqu_era5 * 100, 1)
rsq_imerg_per = round(rsqu_imerg * 100, 1)
rsq_GSMaP_per = round(rsqu_GSMaP * 100, 1)

# Create a time series plot
figure, axes = plt.subplots(nrows=5, sharex=True, figsize=(10, 6))
gsod_ax, chirps_ax, era5_ax, imerg_ax, gsmap_ax = axes

# Plot GSOD precipitation
gsod_ax.plot(gsod_tt.index, gsod_tt, label='GSOD', color='red', linewidth=0.8)
gsod_ax.text(0.4, 0.7, 'GSOD', transform=gsod_ax.transAxes, fontsize=12, fontweight='bold')
#gsod_ax.set_title('GSOD Precipitation (mm)', loc='left')  # Add y-axis label as a heading for GSOD

# Plot CHIRPS precipitation
chirps_ax.plot(chirps_tt.index, chirps_tt, label='CHIRPS', color='black', linewidth=0.8)
chirps_ax.text(0.4, 0.7, 'CHIRPS', transform=chirps_ax.transAxes, fontsize=12, fontweight='bold')
chirps_ax.text(0.8, 0.7, f'Confidence: {rsq_chirps_per}%', transform=chirps_ax.transAxes, fontsize=8)

# Plot ERA5 precipitation
era5_ax.plot(era5_tt.index, era5_tt, label='ERA5', color='green', linewidth=0.8)
era5_ax.text(0.4, 0.7, 'ERA5', transform=era5_ax.transAxes, fontsize=12, fontweight='bold')
era5_ax.text(0.8, 0.7, f'Confidence: {rsq_era5_per}%', transform=era5_ax.transAxes, fontsize=8)

# Plot IMERG precipitation
imerg_ax.plot(imerg_tt.index, imerg_tt, label='IMERG', color='grey', linewidth=0.8)
imerg_ax.text(0.4, 0.7, 'IMERG', transform=imerg_ax.transAxes, fontsize=12, fontweight='bold')
imerg_ax.text(0.8, 0.7, f'Confidence: {rsq_imerg_per}%', transform=imerg_ax.transAxes, fontsize=8)

# Plot GSMaP precipitation
gsmap_ax.plot(GSMaP_tt.index, GSMaP_tt, label='GSMaP', color='blue', linewidth=0.8)
gsmap_ax.text(0.4, 0.7, 'GSMaP', transform=gsmap_ax.transAxes, fontsize=12, fontweight='bold')
gsmap_ax.text(0.8, 0.7, f'Confidence: {rsq_GSMaP_per}%', transform=gsmap_ax.transAxes, fontsize=8)

# Set the y-axis limits to 0 and 200
for ax in axes:
    ax.set_ylim(0, 200)

# Add title and labels
figure.suptitle('Temporal Variations of Precipitation')
gsmap_ax.set_xlabel('Date')
era5_ax.set_ylabel('Total Precipitation (mm)')

# Save the plot as an image file
figure.savefig('time_series_precipitation_plot.png', dpi=300)  # Saved to working directory /Users/pamelaacheng/Landslides

# Show the plot
plt.show()

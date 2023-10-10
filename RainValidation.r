# Load necessary packages
library(rnoaa)
library(ncdf4)
library(tidyverse)

# STEP 01: EXTRACT GROUND OBSERVATIONS FROM GSOD
# Set years of interest
years <- 2014:2020
options(noaakey = 'yEUSoHmWxHlaFqJSDpSKtltllrJKjMLC') # API key to access NOAA data

# Define function to extract GSOD rainfall data for Uganda
get_gsod_data <- function(years, country) {
  
  # Get GSOD data for Uganda and years of interest
  gsod_data <- ncdc(datasetid = "GSOD", 
                    stationid = "", 
                    startdate = paste(years[1], "0101", sep = ""), 
                    enddate = paste(years[length(years)], "1231", sep = ""),
                    limit = 1000,
                    datatypeid = "PRCP",
                    locationid = paste("COUNTRY:UG", sep = ""))
  
  # Filter for rainfall data and convert to data frame
  rain_data <- gsod_data$data %>% 
    select(station, date, value) %>% 
    as.data.frame()
  
  # Rename columns
  colnames(rain_data) <- c("station_id", "date", "gsod_rainfall_in")
  
  # Convert date column to Date format
  rain_data$date <- as.Date(rain_data$date, format = "%Y-%m-%d")
  
  # Convert GSOD rainfall from inches to meters
  rain_data$gsod_rainfall_m <- rain_data$gsod_rainfall_in * 0.0254
  
  # Convert rain_data to a data frame
  rain_data <- as.data.frame(rain_data)
  
  # Return rainfall data
  return(rain_data)
}


#STEP 02: EXTRACT PRECIPITATION FROM SATELLITE-BASED PACKAGES (ERA5, CHIRPS, IMERG)
# Define function to extract ERA5 rainfall data
get_era5_data <- function(years, country) {
  
  # Set base URL for ERA5 data
  base_url <- "https://cds.climate.copernicus.eu/api/v2/resources/era5-land-monthly"
  
  # Set up CDS API request
  era5_request <- list(
    "variable" = "total_precipitation",
    "year" = years,
    "month" = "1/12",
    "format" = "netcdf",
    "area" = paste("box(", 
                   as.numeric(country$lon) - 2, ",",
                   as.numeric(country$lat) - 2, ",",
                   as.numeric(country$lon) + 2, ",",
                   as.numeric(country$lat) + 2, ")"),
    "grid" = "0.25/0.25"
  )
  
  # Send CDS API request and download data
  era5_data <- tempfile()
  curl_download(url = paste(base_url, "?", names(era5_request), "=", era5_request, collapse = "&", sep = ""), 
                destfile = era5_data)
  
  # Open the NetCDF file
  nc <- nc_open(era5_data)
  
  # Extract rainfall data from NetCDF file
  era5_rainfall <- ncvar_get(nc, "tp")
  
  # Extract time data from NetCDF file
  era5_time <- ncvar_get(nc, "time")
  
  # Convert time data to Date format
  era5_dates <- as.Date("1900-01-01") + era5_time/24/60/60
  
  # Close the NetCDF file
  nc_close(nc)
  
  # Convert rainfall data to data frame
  era5_data <- data.frame(date = era5_dates, era5_rainfall) %>% 
    mutate(era5_rainfall = era5_rainfall/1000) %>% 
    filter(year(date) %in% years) %>% 
    group_by(date) %>% 
    summarize(era5_rainfall = sum(era5_rainfall))
  
  # Rename columns
  colnames(era5_data) <- c("date", "era5_rainfall")
  
  # Convert date column to Date format
  era5_data$date <- as.Date(era5_data$date)
  
  # Display first 10 rows of data
  head(era5_data, 10)
}

# STEP 03: COMPARE THE GROUND AND SATELLITE-BASED RAINFALL DATA
# Get rainfall data from GSOD and ERA5 for Uganda
uganda <- list(lon = 32.5, lat = 1)

gsod_data <- get_gsod_data(years, uganda)

era5_data <- get_era5_data(years, uganda)
  
  # Join ERA5 and GSOD data on date
  prcp_data <- inner_join(rain_data, era5_data, by = "date")
  
  # Calculate correlation coefficient
  correlation <- cor(prcp_data$gsod_rainfall, prcp_data$era5_rainfall)
  
  # Calculate RMSE
  rmse <- sqrt(mean((gsod_data$gsod_rainfall - prcp_data$era5_rainfall)^2))
  
  # Display results
  cat("Correlation coefficient:", correlation, "\n")
  cat("RMSE:", rmse, "\n")
  
  

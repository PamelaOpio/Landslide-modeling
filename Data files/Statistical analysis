# 1. Load the raster library
library(raster)
install.packages("rasterVis")
library(rasterVis)
library(ggplot2)
library(stars)

# 2. DESCRIPTIVE STATISTICS

# 2.1 Minimum, mean and maximum values
# Read in the raster file
raster_file <- raster("path/to/raster_file.tif")
# Get the minimum value
setMinMax(raster_DEM_Ug) #obtains the min and max for the elevations
setMinMax(raster_slope) #obtains the min and max of slope
setMinMax(raster_aspect) #obtains the min and max of slope
# Get the mean value
mean_rainfall <- cellStats(raster_rainfall, stat = "mean")
print(mean_rainfall)

mean_elevation <- cellStats(raster_DEM_Ug, stat = "mean")
print(mean_elevation)

# 2.2 Get the standard deviation
std_dev_rain <- cellStats(raster_rainfall, stat = "sd")
print(std_dev_rain)


# 2.3 Create histogram
# Load the raster file into R #raster_DEM_Ug <- raster("/Users/pamelaacheng/Library/CloudStorage/OneDrive-Nexus365/DPhil/2022/Landslides/Data/DEM/DEM_Uganda.tif")
# Extract the values from the raster #elevation_values <- values(raster_DEM_Ug)

# a. Remove NA values from elevation
elevation_values <- elevation_values[!is.na(elevation_values)]
# Plot the histogram
ggplot(data.frame(elevation_values), aes(x=elevation_values)) + 
  geom_histogram(binwidth=100, color="black", fill="white") + 
  ggtitle("Histogram of Elevation Values") + 
  xlab("Elevation (m)") + 
  ylab("Frequency")

# b. Remove NA values from slope
slope_values <- slope_values[!is.na(slope_values)]
# Plot the histogram
ggplot(data.frame(slope_values), aes(x=slope_values)) + 
  geom_histogram(binwidth=3000, color="black", fill="white") + 
  ggtitle("Histogram of Slope Values") + 
  xlab("Slope (m)") + 
  ylab("Frequency")


# 2.4 Visualize data
#a. Level plot of DEM data
plot(raster_slope)
levelplot(raster_slope)
levelplot(raster_aspect)
levelplot(raster_elevation)

# b. Soil data and lithology
levelplot(raster_soil)

# c. LULC data
levelplot(LULC_1992)

#d. Rainfall data
levelplot(raster_rain)

# 3. Landslide data
landslides = gpd.read_file('/Users/pamelaacheng/Library/CloudStorage/OneDrive-Nexus365/DPhil/2022/Landslides/Data/Inventory/Ug_Landslides_updated.csv', 
                           layer='KGS_landslide_inventory_data')
landslides = landslides[~landslides['Failure_Type'].isin(['rockfall', 'rockslide'])]
landslides = landslides.dropna(subset=['FailureDate'])
landslides['FailureDate'] = pd.to_datetime(landslides['FailureDate'])
landslides.sort_values('FailureDate')
landslides = landslides[landslides['FailureDate'].between('2000-01-01', '2021-01-01')]
landslides = landslides.rename(columns={'FailureDate': 'time', 'Longitude83': 'lon', 'Latitude83': 'lat'})
landslides

# 2. INFERENTIAL STATISTICS
#The rasterVis package in R provides several functions for plotting inferential statistics of raster data, 
#such as levelplot.raster and facetRaster.
#For example, if you have a raster dataset named raster_data and a corresponding set of p-values named 
# p_values, you can use the levelplot.raster function to create a color-coded plot of the p-values with 
# a significance threshold, such as a p-value of 0.05:

library(rasterVis)
levelplot(p_values, cuts=c(0,0.05,1),col.regions=c("green","red"),main="p-values")

#n this example, values of p_values less than or equal to 0.05 are colored green and greater than 0.05 
#are colored red.
#You can also use the facetRaster function to create a faceted plot of the raster data and the 
# corresponding p-values:

facetRaster(raster_data, p_values, main="Raster Data and P-Values")

#This will create a plot with two panels, one showing the raster data and the other showing the p-values.
# In both examples, you can customize the plot using various arguments, such as col.regions for the 
# color scheme, cuts for the number of levels, and xlab and ylab for the axis labels.
# In addition, you can also use the package latticeExtra to overlay the inferential statistics on the 
# raster plot using the layer function.
library(latticeExtra)
levelplot(raster_data) + layer(sp.points(...,pch=16,col="red"))



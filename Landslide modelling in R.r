# 1. Load required libraries
install.packages("randomForest")
library(randomForest)
library(raster)
install.packages("rgdal")
library(rgdal)
library(sp)
install.packages("gridExtra")
library(gridExtra)
install.packages("readr")
library(readr)
library(stars)
install.packages("sf") #to read landslide shape file
library(sf)
library(ggplot2)
install.packages("dplyr")
library(dplyr)
install.packages("caret") #to check model performance
library(caret)
#library(geosphere)


# 2. Load the static (conditioning factors) raster data into R
raster_DEM_Ug <- raster("/Users/pamelaacheng/Library/CloudStorage/OneDrive-Nexus365/DPhil/2022/Landslides/Data/DEM/DEM_Uganda.tif")
raster_aspect <- raster("/Users/pamelaacheng/Library/CloudStorage/OneDrive-Nexus365/DPhil/2022/Landslides/Data/Slope&Aspect/Ug_Aspect.tif")
raster_plan <- raster("")
raster_profile <- raster("")
raster_slope <- raster("/Users/pamelaacheng/Library/CloudStorage/OneDrive-Nexus365/DPhil/2022/Landslides/Data/Slope&Aspect/Ug_Slope.tif")
raster_lithology <- raster("")
raster_soil <- raster("/Users/pamelaacheng/Library/CloudStorage/OneDrive-Nexus365/DPhil/2022/Landslides/Data/Soil&Lithology/Ug_Soil.tif")
#fault
#rivers
#roads

# 3. EXTRACTING DATA FROM TEH RASTER FILES
# 3.1. Extracting data from DEM file (extract the elevation values, slope values, aspect values, x and y coordinates, cell size, and projection)
elevation_values <- values(raster_DEM_Ug)
slope_values <- values(raster_slope) #slope_values <- terrain(raster_DEM_Ug, opt = "slope")
aspect_values <- values(raster_aspect) #aspect_values <- terrain(raster_DEM_Ug, opt = "aspect")
x_coords <- xFromCol(raster_DEM_Ug)
y_coords <- yFromRow(raster_DEM_Ug)
projection <- projection(raster_DEM_Ug)
LULC_2015_values <- values(LULC_2015) 

# 3.2. Extracting data from the LULC raster file with several years covered in different bands
# Extract the data from all bands and create a raster stack
bands_LULC_Ug <- stack("/Users/pamelaacheng/Library/CloudStorage/OneDrive-Nexus365/DPhil/2022/Landslides/Data/Land Use & cover/Ug_LULC_ESACCI_1992_2015.tif")
# Get the number of bands in the raster
nbands_LULC_Ug <- nlayers(bands_LULC_Ug)
# Extract each band
LULC_1992 <- bands_LULC_Ug[[1]]
LULC_1993 <- bands_LULC_Ug[[2]]
LULC_2015 <- bands_LULC_Ug[[24]]
# Plot the data
plot(LULC_1992)
plot(LULC_1993)
plot(LULC_2015)

# Load the CSV file containing the legend information
csv_LULC <- read.csv("/Users/pamelaacheng/Library/CloudStorage/OneDrive-Nexus365/DPhil/2022/Landslides/Data/Land Use & cover/ESACCI-LC-Legend.csv")

# Create legend for raster file
#legend_LULC <- levelplot(LULC_1992, 
#main = "Raster File Legend", 
#colorkey = TRUE, 
#scales = list(draw = TRUE))

# Add csv legend to raster legend
#legend_grid <- rbind(legend_LULC, csv_LULC, size = "last")
# Plot raster file with legend
#grid.arrange(legend_grid, LULC_1992, ncol = 2)


# 4. LOAD THE DYNAMIC (TRIGGER) RASTER DATA
raster_LULC_Ug <- raster("/Users/pamelaacheng/Library/CloudStorage/OneDrive-Nexus365/DPhil/2022/Landslides/Data/Land Use & cover/Ug_LULC_ESACCI_1992_2015.tif")
raster_rainfall <-raster("")

# 5. SET DATA WITH LANDSLIDES AND NO LANDSLIDES
# 5.1 Load the landslide inventory data (csv file)
#landslide_data <-read.csv("/Users/pamelaacheng/Library/CloudStorage/OneDrive-Nexus365/DPhil/2022/Landslides/Data/Inventory/Ug_Landslides_updated.csv")
landslides <- st_read("/Users/pamelaacheng/Library/CloudStorage/OneDrive-Nexus365/DPhil/2022/Landslides/Data/Inventory/UgandaLandslide.shp")
raster_landslide <-raster("/Users/pamelaacheng/Library/CloudStorage/OneDrive-Nexus365/DPhil/2022/Landslides/Data/Inventory/Historical landslides/Ug_Landslide_COG-1980_2018.tif")
landslide_df <-read.csv("/Users/pamelaacheng/Library/CloudStorage/OneDrive-Nexus365/DPhil/2022/Landslides/Data/Inventory/Ug_Landslides_13Feb23.csv")
#landslides <- readOGR("/Users/pamelaacheng/Library/CloudStorage/OneDrive-Personal/QGIS/Ug_Landslides.shp")

class(landslide_df) #This checks the data class
str(landslide_df) #Obtains names of columns
str(landslide_df$landslide_location)# Check the type of data in the column
head(landslide_df$landslide_location) #Check the type of data in the column

# 5.2 Create a SpatialPointsDataFrame from the landslides inventory data
landslides_spdf <- SpatialPointsDataFrame(coords = cbind(landslide_df$long, landslide_df$lat),
                                          data = landslide_df,
                                          proj4string = CRS("+proj=longlat +datum=WGS84"))
# 5.3 Define the spatial range for the random points using the minimum and maximum values for the x and y coordinates of Uganda
longmin <- 29.57
longmax <- 35.03
latmin <- -1.48
latmax <- 4.24
# Create a bounding box object representing Uganda
bbox <- matrix(c(longmin, latmin, longmax, latmax), ncol=2, byrow=TRUE)
#colnames(bbox) <- c("longmin", "longmax")
#rownames(bbox) <- c("latmin", "latmax")
# Generate random points within the extent of Uganda
random_points <- spsample(bbox, n=1000, type="random")
# Convert the random points to a spatial points data frame
#random_points_spdf <- SpatialPointsDataFrame(random_points, data=data.frame(id=1:1000))
# Convert the random points to a data frame
random_points_df <- data.frame(long = coordinates(random_points)[,1],
                               lat = coordinates(random_points)[,2])
# Add a column to the random points data frame indicating that no landslide occurred at these locations
random_points_df$landslide <- 0
# Merge the random points data frame with the landslides inventory data frame
merged_df <- rbind(landslide_df, random_points_df)

# 5.4 Extract values from raster for the landslide and no landlside locations
xy_coords <- merged_df[, c("long", "lat")]

model_slopes_values <- extract(raster_slope, xy_coords)
merged_df$slope <- model_slopes_values

model_elevation_values <- extract(raster_DEM_Ug, xy_coords)
merged_df$elevation <- model_elevation_values

model_soil_values <- extract(raster_soil, xy_coords)
merged_df$soil <- model_soil_values

#do the same for lithology, aspect, profile and plan
#need to be creative with rainfall data
#what to do with distance to faults
#what to do with distance to roads
#what to do with distance to rivers


# 6. MODEL LANDSLIDES

# 6.1 First check that the class and CRS of the landslide data and other raster files match
class(landslides)#Check class of landlside data
class(raster_slope)#Check class of landlside data
crs(raster_slope)#Check CRS of the slope raster
crs(landslides)#Check CRS of the landslides data
landslides_reproj <- spTransform(landslides, crs(raster_slope))# Reproject landslides inventory CRS to match with slope if need be

# 6.2 Split the data into training and testing sets
set.seed(123)
ind <- sample(2, nrow(merged_df), replace=TRUE, prob=c(0.7, 0.3))
train <- merged_df[ind==1,]
test <- merged_df[ind==2,]

# 6.3 Build the random forest model
train <- na.omit(train)
model <- randomForest(landslide ~ ., data=train)

# 6.4 Make predictions using the model
predictions <- predict(model, newdata=test)
#predictions <- predict(model, test)
confusionMatrix(predictions, test$triggered) #Evaluate the model performance


# 6.5 Calculate the model's accuracy
accuracy <- mean(predictions == test$landslide)
print(accuracy)

# 6.6 Save the model to a file
saveRDS(model, file = "landslide_model.rds")



# 7. USE MODEL TO PREDICT FUTURE EVENTS
# 7.1 Load future climate and LULC data

# 7.2 Use future climate and LULC data to predict Landslide events






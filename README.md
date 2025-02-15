# Introduction

As climate change increases in severity, so does the occurrence of natural disasters (Vernick, 2025). In the Canadian context, due to high amounts of fuel from forests and landscapes, we are at increasing risk of wildfire. As recently as 2023, British Columbia experienced its worst fire year on record (“2023 Is Now Officially the Most Expensive, Most Destructive Wildfire Season on Record in B.C.,” 2023). 

To improve wildfire predictions, we can use climate indicators to assess upcoming wildfire seasons. Our study focuses on British Columbia, a region with diverse climate zones. In this project, we analyze snow depth patterns and their correlation with wildfire density. Snow accumulation, measured in centimeters, can impact soil moisture, vegetation growth, and fuel availability, all of which influence wildfire potential (Westerling et al., 2006). 

This tutorial explores the hypothesis that higher snow depth in winter leads to lower wildfire density in the following summer. Research suggests that increased snow accumulation delays the onset of spring, shortening the fire season (Holden et al., 2018). Conversely, in years with less snow accumulation, an early spring can lead to drier vegetation and increased wildfire risk (Holden et al., 2018).

Using R, we will conduct various statistical and spatial analyses to assess this relationship. By the end of this tutorial, users will be able to perform:
- Data cleaning and aggregation
- Descriptive and inferential statistics
- Point pattern analysis
- Spatial autocorrelation & interpolation
- Geographically Weighted Regression (GWR)

Throughout this process, users will generate maps, tables, and charts to visualize findings and contribute to wildfire risk assessment efforts.

# Setting Up Your Code
Firstly we need to prepare our R project workspace, we do this by installing packages we will need for the project so the program has the tools and functions they provide. We will also load in a BC shapefile to properly visualize our maps. Then we will call our packages into the project with the library function. Next we will set up the directory so the project knows the folder we are working within and defining file paths to any necessary data. Lastly we will call in and set our boundary shapefile and set a projection. In this case since our project is centered around British Columbia and we would like to preserve area we will be using the Albers projection.

#Getting Your Data
Our data for this example came from multiple sources. The wildfire points were downloaded from provided by Professor Bone and I believe were originally sourced from the BC Wildfire Service. The climate data was downloaded from the PCIC (Pacific Climate Impacts Consortium) and the filtered to use the BC Hydro weather stations within a time frame.
```r
# Install necessary packages
install.packages("gt")
install.packages("moments")
install.packages("dplyr")
install.packages("lubridate")
install.packages("sf")
install.packages("ggplot2")
install.packages("gstat")
install.packages("sp")
install.packages("raster")
install.packages("tmap")
install.packages("viridis")
install.packages("bcmaps")
install.packages("nortest")
install.packages("spatstat")
install.packages("spatstat.geom")
install.packages("nngeo")
install.packages("RANN")
install.packages("spatstat.explore")
install.packages("spdep")
install.packages("MASS")

# Install BC Map (from GitHub)
remotes::install_github("bcgov/bcmaps")  

# Load in libraries for project
library(gt)
library(moments)
library(dplyr)
library(lubridate)
library(sf)
library(ggplot2)
library(gstat)
library(sp)
library(raster)
library(tmap)
library(viridis)
library(bcmaps)  
library(nortest)
library(spatstat)
library(spatstat.geom)
library(nngeo)  
library(RANN)    
library(spatstat.explore)
library(spdep)
library(MASS) 

# Set Project Directory
setwd("C:/Users/admin/Desktop/Project/FinalProject/")

# Define file paths
csv_folder <- "C:/Users/admin/Desktop/Project/FinalProject/BCH"
metadata_file <- "C:/Users/admin/Desktop/Project/FinalProject/station-metadata-by-history.csv"
output_csv <- "Processed_Snow_Depth.csv"
final_shapefile <- "Cleaned_SnowDepth_Stations.shp"
final_csv <- "Cleaned_SnowDepth_With_Locations.csv"

#Load boundary and create shapefile
bc_boundary <- bc_bound() 
# Put boundary in albers Projection & Check it worked
st_crs(bc_boundary) <- 3005 
print(st_crs(bc_boundary)) 

# Save the BC Boundary Shapefile for reference
st_write(bc_boundary, "C:/Users/admin/Desktop/Project/FinalProject/BC_Boundary.shp", delete_dsn = TRUE)
```

# Cleaning Your Data
The next step in the project is 'cleaning' our data to make it usable for us. This step makes sure that the snow depth data from our weather stations is formatted correctly. We then extract the data from each individual station CSV. The data is then filtered to remove irrelevant stations. Then we aggregate it together into a single CSV. We also set our study period by choosing only winter months of September - March when snow may be present to increase efficiency and average it for the entire season. 

```
#list files within the folder
csv_files <- list.files(path = csv_folder, pattern = "\\.csv$", full.names = TRUE)

#Defines a function and prints the name of the file input
process_snow_depth <- function(file) {
  cat("Processing file:", file, "\n")
  
  # Read CSV and checks for errors
  station_data <- tryCatch({
    read.csv(file, skip = 1, header = TRUE, na.strings = c("None", "NA", ""))
  }, error = function(e) {
    cat("Error reading file:", file, "\n")
    return(NULL)
  })
  
  # Checks if snow depth column exists
  if (!"SNOW_ON_THE_GROUND" %in% colnames(station_data)) {
    cat("Skipping file (no snow depth data):", file, "\n")
    return(NULL)
  }
  
  # Convert times and filters data
  station_data$time <- ymd_hms(station_data$time)
  station_data$SNOW_ON_THE_GROUND <- as.numeric(station_data$SNOW_ON_THE_GROUND)
  station_data <- station_data %>% filter(SNOW_ON_THE_GROUND >= 0 & SNOW_ON_THE_GROUND <= 200)
  
  # Select winter months (Sept–March) for seasonal analysis
  winter_data <- station_data %>%
    filter(month(time) %in% c(9, 10, 11, 12, 1, 2, 3)) %>%
    filter(!is.na(SNOW_ON_THE_GROUND))  
 
   #if there is the value of snow on the ground is 0 it skips that station
  if (nrow(winter_data) == 0) {
    cat("Skipping file (no valid winter data):", file, "\n")
    return(NULL)
  }
  
  # calculates seasonal snow depth
  seasonal_snow <- mean(winter_data$SNOW_ON_THE_GROUND, na.rm = TRUE)
  station_id <- tools::file_path_sans_ext(basename(file))
  
  return(data.frame(
    Native.ID = station_id,
    SEASONAL_SNOW_DEPTH = seasonal_snow
  ))
}
#combines snow depth from all csvs into single data frame
snow_depth_data <- do.call(rbind, lapply(csv_files, process_snow_depth))

#Checks if output is successful
if (!is.null(snow_depth_data)) {
  write.csv(snow_depth_data, file = output_csv, row.names = FALSE)
  cat("\nSaved snow depth data.\n")
} else {
  cat("\nNo snow depth data found.\n")
}
```
# Mapping Climate Data
Now we will map our newly aggregated snow depth data so we can see the stations average seasonal snow depth. This This code will set a boundary, load the data, clip it, check for errors that may cause issue and then generate a map with customizations for legibility. Lastly it will save the map as an output. 

![Snow Depth Map](https://github.com/JColalillo/Spatial-Analysis-Final/blob/main/Snow_Depth_Points_Map_Fixed.png?raw=true)

Now that we have those mapped out we can begin our analysis through Inverse Distance Weighted Interpolation. This is a spatial interpolation technique which can estimate the values across the entire province by using the data at the points of each station. It assumes closer points will have a greater effect on the prediction than those further away. The weighting of each point will decrease as the distance to an unknown area grows. IDW is calculated by using a weighted average of close points using this equation

![IDW Equation](![image](https://github.com/user-attachments/assets/2a05ac2b-fb49-444d-b08a-d67511228506)
```
# Load Data
snow_points <- st_read("Cleaned_SnowDepth_Stations.shp")
bc_boundary <- st_read("BC_Boundary.shp")  #Load BC boundary shapefile

# Double check CRS again
target_crs <- 3005
if (st_crs(snow_points) != target_crs) {
  snow_points <- st_transform(snow_points, crs = target_crs)
}
if (st_crs(bc_boundary) != target_crs) {
  bc_boundary <- st_transform(bc_boundary, crs = target_crs)
}

# Rename Column
colnames(snow_points)[colnames(snow_points) == "SEASONA"] <- "SEASONAL_SNOW_DEPTH"

# Checks Numeric Snow Depth Values & Remove NAs
snow_points <- snow_points %>% filter(!is.na(SEASONAL_SNOW_DEPTH))
snow_points$SEASONAL_SNOW_DEPTH <- as.numeric(snow_points$SEASONAL_SNOW_DEPTH)

# Define Grid Resolution for Interpolation (10 km)
grid_res <- 10000

# Create Grid
bbox <- st_bbox(bc_boundary)  #Use BC boundary extent, NOT snow points
grid <- expand.grid(
  x = seq(bbox["xmin"], bbox["xmax"], by = grid_res),
  y = seq(bbox["ymin"], bbox["ymax"], by = grid_res)
)
coordinates(grid) <- ~x + y
proj4string(grid) <- CRS(st_crs(target_crs)$proj4string)

# Convert Grid to SF Format
grid_sf <- st_as_sf(grid, coords = c("x", "y"), crs = target_crs)

# Perform IDW Interpolation
idw_model <- gstat(
  formula = SEASONAL_SNOW_DEPTH ~ 1, 
  locations = as(snow_points, "Spatial"), 
  nmax = 15
)
idw_result <- predict(idw_model, grid)

# Convert IDW Result to Raster
idw_raster <- rasterFromXYZ(as.data.frame(idw_result)[, c("x", "y", "var1.pred")])
crs(idw_raster) <- target_crs

# Convert BC boundary to Spatial
bc_boundary_sp <- as(bc_boundary, "Spatial")

# Mask & Crop
idw_raster_clipped <- mask(idw_raster, bc_boundary_sp)  #Clip to BC boundary
idw_raster_clipped <- crop(idw_raster_clipped, extent(bc_boundary_sp))  # final cropping

# Convert NA values to transparent
idw_raster_clipped[is.na(idw_raster_clipped)] <- NA

# Save the Clipped Raster 
writeRaster(idw_raster_clipped, "IDW_SnowDepth_Clipped.tif", format = "GTiff", overwrite = TRUE)

# Convert Raster to DataFrame for Mapping
idw_df <- as.data.frame(idw_raster_clipped, xy = TRUE) %>% filter(!is.na(var1.pred))
colnames(idw_df) <- c("x", "y", "snow_depth")
idw_sf <- st_as_sf(idw_df, coords = c("x", "y"), crs = target_crs)

# Save as Shapefile
st_write(idw_sf, "IDW_SnowDepth_Clipped.shp", delete_dsn = TRUE)

# Generate Final Snow Depth Map 
snow_depth_map <- ggplot() +
  geom_raster(data = idw_df, aes(x = x, y = y, fill = snow_depth)) +  # IDW Raster
  geom_sf(data = bc_boundary, fill = NA, color = "white", linewidth = 0.8) +  # BC Boundary
  scale_fill_viridis_c(
    option = "plasma",
    name = "Snow Depth (cm)",
    breaks = scales::pretty_breaks(n = 5)
  ) +
  labs(
    title = "Interpolated Snow Depth (IDW)",
    subtitle = "Inverse Distance Weighted (IDW) Interpolation",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    plot.background = element_rect(fill = "black"),  # Set background
    panel.background = element_rect(fill = "black"),  # Darken panel
    plot.title = element_text(face = "bold", size = 16, color = "white"),  #title
    plot.subtitle = element_text(size = 12, color = "white"),  # White subtitle
    axis.title = element_text(size = 12, color = "white"),  #axis labels
    axis.text = element_text(size = 10, color = "white"),  #axis text
    legend.text = element_text(size = 10, color = "white"),  #legend text
    legend.title = element_text(size = 12, face = "bold", color = "white")  #legend title
  )

# Save Updated Map
ggsave("IDW_SnowDepth_Final.png", plot = snow_depth_map, width = 10, height = 8, dpi = 300)

cat("\nIDW Snow Depth Map Created\n")
```
#IDW Output Map
The map below provides the estimated snow depth across the entire province through a continuous representation. The areas which are estimated to receieve more snow are in yellow & orange and the lesser areas in purple.

![IDW Map of Estimated Snow Accumulations in British Columbia 2023](https://github.com/JColalillo/Spatial-Analysis-Final/blob/main/IDW_SnowDepth_Final.png?raw=true)

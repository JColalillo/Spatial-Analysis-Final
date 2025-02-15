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

![IDW Equation](https://github.com/user-attachments/assets/2a05ac2b-fb49-444d-b08a-d67511228506)

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
# IDW Output Map

The map below provides the estimated snow depth across the entire province through a continuous representation. The areas which are estimated to receieve more snow are in yellow & orange and the lesser areas in purple.

![IDW Map of Estimated Snow Accumulations in British Columbia 2023](https://github.com/JColalillo/Spatial-Analysis-Final/blob/main/IDW_SnowDepth_Final.png?raw=true)

# Descriptive Statistics
Now we will look at our event data which is the 2024 British Columbia wildfire data. It is stored in a shapefile with multiple fields that will allow us to calculate descriptive statistics and interpret the data. In this case will examine the size of the fires which is measured by the hectare.

```
# File Path
wildfire_shapefile <- "C_FIRE_PNT_point.shp"

# Load Shapefile
wildfire_data <- st_read(wildfire_shapefile)

# Displays all column names
print(colnames(wildfire_data))  # Displays all column names

# Convert SIZE_HA to Numeric
wildfire_data$SIZE_HA <- as.numeric(wildfire_data$SIZE_HA)

# Remove NAs
wildfire_data <- wildfire_data[!is.na(wildfire_data$SIZE_HA), ]

# Caculate table 1 stats
basic_stats <- data.frame(
  mean = mean(wildfire_data$SIZE_HA, na.rm = TRUE),
  sd = sd(wildfire_data$SIZE_HA, na.rm = TRUE),
  median = median(wildfire_data$SIZE_HA, na.rm = TRUE),
  mode = as.numeric(names(sort(table(wildfire_data$SIZE_HA), decreasing = TRUE)[1]))  # Most Frequent Value
)

# Calculate table 2 stats
advanced_stats <- data.frame(
  skewness = skewness(wildfire_data$SIZE_HA, na.rm = TRUE),
  kurtosis = kurtosis(wildfire_data$SIZE_HA, na.rm = TRUE),
  CoV = sd(wildfire_data$SIZE_HA, na.rm = TRUE) / mean(wildfire_data$SIZE_HA, na.rm = TRUE),
  normality = shapiro.test(wildfire_data$SIZE_HA[1:5000])$p.value  # Shapiro-Wilk Normality Test
)

# Print stats
print(basic_stats)
print(advanced_stats)

# Table 1 image
basic_table <- basic_stats %>%
  gt() %>%
  tab_header(
    title = "Table 1: Descriptive Statistics for Wildfire Size (Ha)",
    subtitle = "Summary of wildfire sizes in BC"
  ) %>%
  cols_label(
    mean = "Mean",
    sd = "SD",
    median = "Median",
    mode = "Mode"
  ) %>%
  fmt_number(columns = everything(), decimals = 2) %>%
  tab_options(
    table.width = px(600)
  )

# Table 2 image
advanced_table <- advanced_stats %>%
  gt() %>%
  tab_header(
    title = "Table 2: Additional Statistics for Wildfire Size (Ha)",
    subtitle = "Distribution properties of wildfire sizes in BC"
  ) %>%
  cols_label(
    skewness = "Skewness",
    kurtosis = "Kurtosis",
    CoV = "CoV",
    normality = "Normality"
  ) %>%
  fmt_number(columns = everything(), decimals = 2) %>%
  tab_options(
    table.width = px(600)
  )

# Print Tables
print(basic_table)
print(advanced_table)
```
We can now interpret our output tables to look for trends within the size of the fires. The mean size of the fires is 640 Ha but we can tell this heavily skewed by larger fires. We know this because the standard deviation of 13 231 Ha indicates very high variability in the size of the fires, this means there are lots of fires either much smaller or larger than the mean. The median being only .1 Ha tells us the majority of wildfires are quite small and are contained or die before they spread. Lastly for our basic statistics the mode is only .01 Ha which also indicates the majority of fires are very small. When considering the skew the value of 34 is highly positive and a small amount of large fires a right skew. The kurtosis value of 1305 is exceptionally high, this tells us most values are near 0 but the outliers are huge. The coefficient of variation  of 20 tells us that the size varies greatly fire to fire. Lastly a normality gives us our p value of 0 which means our size distribution is not normal and heavily skewed.

![Descriptive Statistics on Wildfire Size](https://github.com/JColalillo/Spatial-Analysis-Final/blob/main/WildfireTable1.png?raw=true)
![Descriptive Statistics on Wildfire Size](https://github.com/JColalillo/Spatial-Analysis-Final/blob/main/WildfireTable2.png?raw=true)

#Mapping Our Fire Events
In this section we will map our wildfire events so we can see their distribution across the province. As you can see from the map below the entire province is subject to risk of wildfire but there is a higher frequency as you move east in the province.
```
# Load Wildfire Points Data
C_FIRE_PNT_point <- st_read("C_FIRE_PNT_point.shp")  

# Load BC boundary using bcmaps
bc_boundary <- bc_bound()  

# checks for CRS Consistency (Albers Projection)
target_crs <- 3005  

# Changes CRS as necessary
if (st_crs(C_FIRE_PNT_point) != target_crs) {
  C_FIRE_PNT_point <- st_transform(C_FIRE_PNT_point, crs = target_crs)
}
if (st_crs(bc_boundary) != target_crs) {
  bc_boundary <- st_transform(bc_boundary, crs = target_crs)
}

# Creates Map
fire_points_map <- ggplot() +
  geom_sf(data = bc_boundary, fill = NA, color = "white", linewidth = 1) +  # BC boundary in white for contrast
  geom_sf(data = C_FIRE_PNT_point, color = "red", size = 0.7, alpha = 0.8) +  # Fire points with improved visibility
  theme_minimal(base_size = 14) +  # Increase base font size
  labs(
    title = "Wildfire Locations in British Columbia During 2024",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme(
    plot.background = element_rect(fill = "black"),  # Set background to black
    panel.background = element_rect(fill = "black"),  # Darken panel
    plot.title = element_text(face = "bold", size = 16, color = "white"),  # Title in white
    axis.title = element_text(size = 12, color = "white"),  # Axis labels in white
    axis.text = element_text(size = 10, color = "white"),  # Axis text in white
    panel.grid.major = element_line(color = "gray40", size = 0.3)  # Adjust gridlines for readability
  )

# Save Fire Points Map
ggsave("Fire_Points_Map_Improved.png", plot = fire_points_map, width = 10, height = 8, dpi = 300)
```
![Wildfire Locations 2024](https://github.com/JColalillo/Spatial-Analysis-Final/blob/main/Fire_Points_Map_Improved.png?raw=true)

# Point Pattern Analysis Using Nearest Neighbor
In this section of code we will conduct a point pattern analysis on our wildfire locations using nearest neighbour analysis. This analysis will tell us about the type of distribution of the wildfire events, they can be random, clustered or dispersed in distribution.
```
# Load Wildfire Data
wildfire_shapefile <- "C_FIRE_PNT_point.shp"
wildfire_data <- st_read(wildfire_shapefile)

# Load BC Boundary
bc_boundary_shapefile <- "C:/Users/admin/Desktop/Project/FinalProject/BCShape.shp"
bc_boundary <- st_read(bc_boundary_shapefile)

# Check CRS 
target_crs <- 3005
if (st_crs(wildfire_data) != st_crs(target_crs)) {
  wildfire_data <- st_transform(wildfire_data, crs = target_crs)
}
if (st_crs(bc_boundary) != st_crs(target_crs)) {
  bc_boundary <- st_transform(bc_boundary, crs = target_crs)
}

# Extract wildfire coordinates
wildfire_coords <- st_coordinates(wildfire_data)

# calculate NND
nearest_neighbors <- nngeo::st_nn(wildfire_data, wildfire_data, k = 2, progress = FALSE)

# Get 1st NND
nearest_distances <- sapply(nearest_neighbors, function(x) {
  if (length(x) > 1) {
    st_distance(wildfire_data[x[1], ], wildfire_data[x[2], ])
  } else {
    NA
  }
})

# Convert to DataFrame
wildfire_nn_stats <- data.frame(
  Fire_ID = wildfire_data$FIRE_ID,
  Nearest_Neighbor_Distance = as.numeric(nearest_distances)
)

# Calculate Summ stats
nn_summary <- data.frame(
  Mean_NN = mean(wildfire_nn_stats$Nearest_Neighbor_Distance, na.rm = TRUE),
  Median_NN = median(wildfire_nn_stats$Nearest_Neighbor_Distance, na.rm = TRUE),
  Min_NN = min(wildfire_nn_stats$Nearest_Neighbor_Distance, na.rm = TRUE),
  Max_NN = max(wildfire_nn_stats$Nearest_Neighbor_Distance, na.rm = TRUE),
  SD_NN = sd(wildfire_nn_stats$Nearest_Neighbor_Distance, na.rm = TRUE)
)

print(nn_summary)

# Makes a histogram of NND
nn_plot <- ggplot(wildfire_nn_stats, aes(x = Nearest_Neighbor_Distance)) +
  geom_histogram(bins = 30, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Nearest Neighbor Distances of Wildfires",
       x = "Distance (meters)",
       y = "Frequency") +
  theme_minimal()

# Save the histogram
ggsave("Wildfire_Nearest_Neighbor_Histogram.png", plot = nn_plot, width = 10, height = 8, dpi = 300)

# Calculate NNI
observed_nn <- mean(wildfire_nn_stats$Nearest_Neighbor_Distance, na.rm = TRUE)
area_bc <- as.numeric(st_area(bc_boundary))
num_fires <- nrow(wildfire_data)
expected_nn <- 1 / (2 * sqrt(num_fires / area_bc))
nn_index <- observed_nn / expected_nn  

# Calculate Z-Score
z_score <- (observed_nn - expected_nn) / (0.26136 / sqrt(num_fires))

# Determine Spatial Pattern
spatial_pattern <- ifelse(nn_index < 1, "Clustered",
                          ifelse(nn_index > 1, "Dispersed", "Random"))

# Create NN summary Table
nn_table <- data.frame(
  Observed_Distance = observed_nn,
  Expected_Distance = expected_nn,
  NN_Index = nn_index,
  Z_Score = z_score,
  Pattern = spatial_pattern
)

# Display Table
nn_table %>%
  gt() %>%
  tab_header(
    title = "Nearest Neighbor Summary for Wildfires",
    subtitle = "Analysis of Wildfire Spatial Distribution"
  ) %>%
  fmt_number(columns = c(Observed_Distance, Expected_Distance, NN_Index, Z_Score), decimals = 2) %>%
  cols_label(
    Observed_Distance = "Observed NN Distance (m)",
    Expected_Distance = "Expected NN Distance (m)",
    NN_Index = "Clark & Evans Index",
    Z_Score = "Z-Score",
    Pattern = "Spatial Pattern"
  ) %>%
  tab_options(
    table.width = px(600)
  )
```
The code will print our results in a table. As you can see our Observed NN distance is 7503 which tells us the average distance between each fire is 7503 meters. The value of for 11 666 for expected NN distance tells us if the fires were random  then the expected average distance between them would be 11 666m. lastly our Nearest Neighbour Index value (C&E Index) is .64 we interpret this as clustered. a value = 1 would be random, less than 1 clustered and above 1 dispersed. A value of .64 tells us the fires are very clustered by nature. Lastly our very negative z score of -653 756 tells us that this did not occur by accident.

![Nearest Neighbor Data](https://github.com/JColalillo/Spatial-Analysis-Final/blob/main/NearestNeighborWildfire.png?raw=true)

# Descriptive Statistics Using Kernel Density Estimates (KDE)
In this section we will complete descriptive statistics using KDE. KDE creates a surface that highlights areas of higher or lower concentration of event occurences, in this case wildfires. It estimates the probability of density of a dataset and visualizes it on a continous surface instead of points. It gives higher values to area with more events
and uses a bandwidth to determine a radius of how far to smooth the values. The scale on this map shows relative density. As you can see the South East of the province has the highest density of fires but there are other notable clusters around the province. We can also note that the northern areas of the province have significantly lower fire density. 

```
# Load Wildfire Data
wildfire_shapefile <- "C_FIRE_PNT_point.shp"
wildfire_data <- st_read(wildfire_shapefile)

# Load BC Boundary
bc_boundary_shapefile <- "BC_Boundary.shp"
bc_boundary <- st_read(bc_boundary_shapefile)

# Ensure Correct CRS (BC Albers Projection)
target_crs <- 3005

if (st_crs(wildfire_data) != target_crs) {
  wildfire_data <- st_transform(wildfire_data, crs = target_crs)
}

if (st_crs(bc_boundary) != target_crs) {
  bc_boundary <- st_transform(bc_boundary, crs = target_crs)
}

# Extract wildfire coordinates
wildfire_coords <- st_coordinates(wildfire_data)

# Calculate Kernel Density Estimate (KDE)
bbox <- st_bbox(bc_boundary)

kde_result <- kde2d(
  x = wildfire_coords[, 1], 
  y = wildfire_coords[, 2], 
  n = 200,  # Higher resolution for better visualization
  lims = c(bbox["xmin"], bbox["xmax"], bbox["ymin"], bbox["ymax"])
)

# Convert KDE result to raster
kde_raster <- raster(kde_result)
crs(kde_raster) <- st_crs(target_crs)$proj4string

# Clip KDE raster to BC boundary
kde_raster_clipped <- mask(crop(kde_raster, bc_boundary), bc_boundary)

# Normalize KDE for better visualization
kde_raster_clipped <- calc(kde_raster_clipped, function(x) x / max(x, na.rm = TRUE))

# Convert raster to dataframe for plotting
kde_df <- as.data.frame(rasterToPoints(kde_raster_clipped))
colnames(kde_df) <- c("x", "y", "density")

# Create KDE Map with Fixed Scale Bar
kde_plot <- ggplot() +
  geom_raster(data = kde_df, aes(x = x, y = y, fill = density), alpha = 0.85) +  # Adjust alpha for better contrast
  geom_sf(data = bc_boundary, fill = NA, color = "white", linewidth = 1.5) +  # Thicker boundary for visibility
  scale_fill_viridis_c(
    name = "Wildfire Density",
    option = "magma",  # Higher contrast color scale
    breaks = scales::pretty_breaks(n = 5),  # Ensures evenly spaced values
    guide = guide_colorbar(barwidth = 12, barheight = 0.8)  # Adjusts scale bar size
  ) +
  labs(
    title = "Wildfire Density in British Columbia",
    subtitle = "Kernel Density Estimation (KDE)",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal(base_size = 16) +  # Increase text size for readability
  theme(
    legend.position = "bottom",
    legend.key.width = unit(2, "cm"),  # Makes legend bar wider
    legend.key.height = unit(0.5, "cm"),  # Reduces height
    legend.text = element_text(size = 12, color = "white"),  # Larger legend labels
    legend.title = element_text(size = 14, face = "bold", color = "white"),  
    plot.title = element_text(face = "bold", size = 18, color = "white"),
    plot.subtitle = element_text(size = 14, color = "white"),
    axis.title = element_text(size = 14, color = "white"),
    axis.text = element_text(size = 12, color = "white"),
    plot.background = element_rect(fill = "black"),  
    panel.background = element_rect(fill = "black"),
    panel.grid.major = element_line(color = "gray40", size = 0.3)
  )

# Save the improved KDE map
ggsave("Wildfire_KDE_Map_Fixed_Scale.png", plot = kde_plot, width = 12, height = 8, dpi = 300)

# Display KDE plot
print(kde_plot)
```
![Wildfire KDE of British Columbia 2024](https://github.com/JColalillo/Spatial-Analysis-Final/blob/main/Wildfire_KDE_Map_Fixed_Scale.png?raw=true)

# Ordinary Least Squares Regression (OLS)
OLS allows us to to see any intiial indication of a relationship between snow depth and fire density. It will also let us know if further modelling such as geographical weighted regression (GWR) is necessary. OLS allows to measure how changes in snow depth (our independent variable) affect our wildfire density (dependent variable). It uses the following equation 
![OLS Equation](https://github.com/JColalillo/Spatial-Analysis-Final/blob/main/OLSEQ.png?raw=true)
To interpret our results a negative coefficient indicates higher snow depth results in lower wildfire density and a positive coefficient meaning more snow would lead to more fires. The R value if closer to 1 indicates a stronger relationship and 0 indicates little to no correlation. lastly the P-value tests for significance and a value under 0.05 would indicate this is the case.
We can complete our OLS with the following code
```
#Load Data
final_data_sf <- st_read("final_data_with_residuals.shp")
bc_boundary <- bc_bound()  # Load BC boundary from bcmaps

#Check CRS
target_crs <- 3005  # BC Albers Projection
if (st_crs(final_data_sf) != target_crs) {
  final_data_sf <- st_transform(final_data_sf, crs = target_crs)
}
if (st_crs(bc_boundary) != target_crs) {
  bc_boundary <- st_transform(bc_boundary, crs = target_crs)
}

# Define Breaks for Legend
residual_range <- range(final_data_sf$residuals, na.rm = TRUE)
residual_breaks <- seq(residual_range[1], residual_range[2], length.out = 5)  # Generate 5 breaks
residual_labels <- c("Strong Underprediction", "Slight Underprediction", "Neutral", "Slight Overprediction", "Strong Overprediction")

# Create Residuals Map
residuals_map <- ggplot() +
  geom_sf(data = final_data_sf, aes(fill = residuals), color = NA) +  # Rasterized residuals
  geom_sf(data = bc_boundary, fill = NA, color = "white", linewidth = 1.2) +  # Overlay BC boundary
  scale_fill_viridis_c(
    option = "plasma",
    name = "Residuals (Fire Density - Predicted)",
    breaks = residual_breaks, 
    labels = residual_labels
  ) +
  theme_minimal(base_family = "Arial") +  #font
  labs(
    title = "Residuals from OLS Regression: Snow Depth vs. Fire Density",
    subtitle = "Negative = Model Underpredicted, Positive = Model Overpredicted",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 16, color = "white"),  # Title
    plot.subtitle = element_text(size = 12, color = "white"),  # Subtitle
    axis.title = element_text(size = 12, color = "white"),  # Axis labels
    axis.text = element_text(size = 10, color = "white"),  # Axis text
    legend.text = element_text(size = 10, color = "white"),  # Legend text
    legend.title = element_text(size = 12, face = "bold", color = "white"),  # Legend title
    panel.background = element_rect(fill = "black"),  # Background
    plot.background = element_rect(fill = "black")  # Set full plot background
  )

# Save the Residuals Map
ggsave("Residuals_Map.png", plot = residuals_map, width = 10, height = 8, dpi = 300, bg = "black")
```
Now lets interpret the results of our output. 


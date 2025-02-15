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
```

```![Snow_Depth_Points_Map_Fixed](https://github.com/user-attachments/assets/7ee31eef-351d-4e93-9393-59d14f6e1c19)

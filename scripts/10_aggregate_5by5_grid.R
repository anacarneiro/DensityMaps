## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Mapping the global distribution of seabird populations
## R script to aggregate results into a 5x5 degree grid
## Ana Carneiro and Anne-Sophie Bonnet-Lebrun
## July 2018
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# This code uses as input a 10 km x 10 km raster containing the predicted number of birds in each cell. This raster is 
# first projected into WGS84 (longitude and latitude), with a resolution roughly equivalent to the initial 10 km2 resolution, 
# i.e. 1.25x1.25 degrees. Since 1.25x1.25 degree cells at high latitudes do not have the same area as cells at low latitudes, 
# a spatially heterogeneous rescaling is necessary to avoid distorting the values (number of predicted birds per cell). 
# That rescaling is done by dividing the value of the cell by the area of the original cell (10 km2) and multiplying it by the 
# area of the new cell (which depends on the latitude of the cell). The last step is then to aggregate results into 5x5 degree cells, 
# by summing the values (number of birds) in blocks of 16 1.25x1.25 degree cells.
# 
# To calculate the area of each grid cell, we consider a spherical model of the Earth. The radius of this sphere is the Earth's authalic radius (R = 6371007.2m).
# The curved surface area of a grid cell spanning longitudes lon0 to lon1 (lon1 > lon0) and latitudes lat0 to lat1 (lat1 > lat0) is:
# A = (lon1 - lon0)*(sin(lat1) - sin(lat0)*R^2
# (see answer 2 - https://gis.stackexchange.com/questions/29734/how-to-calculate-area-of-1-x-1-degree-cells-in-a-raster) 
# So for a grid of resolution res_lon and res_lat, a cell centred on lon and lat (all latitudes and resolutions in radians) has an area of:
# A = (lon + res_lon/2 - (lon - res_lon/2))*(sin(lat + res_lat/2) - sin(lat - res_lat/2))*R^2 = res_lon*(sin(lat + res_lat/2) - sin(lat - res_lat/2))


#####################################################
################ LOADING PACKAGES ###################
#####################################################

library(raster)
library(rgdal)



#####################################################
######### GENERAL DIRECTIONS AND FILES ##############
#####################################################

## GENERAL DIR
dir <- "C:/Users/ana.carneiro/Documents/GEF analysis/supplements"  ## copy and paste here your working directory

## DIRECTION TO YOUR RASTERS (ALL DEM CLASSES COMBINED AND BY YEAR QUARTER)
dir_demClasses <- "C:/Users/ana.carneiro/Documents/GEF analysis/supplements/scripts_results/09_sum_demClasses"

## DIRECTION TO YOUR RESULTS
dir_5by5 <- "C:/Users/ana.carneiro/Documents/GEF analysis/supplements/scripts_results/10_aggregate_5by5_grid"


#####################################################
####### CONVERT INTO A 5X5 DEGREE RESOLUTION ########
#####################################################

## resolution and projection
r <- raster(xmn=-180, xmx=180, ymn=-90, ymx=0, ncols=288, nrows=72) 
projection(r)
res(r) # 1.25 x 1.25 degree raster

RES = res(r) # the resolution of the raster (in degrees)
# res_lon = RES[1]*pi/180 (in radians) and res_lat = RES[2]*pi/180 (in radians)
R = 6371007.2 # the Earth's authalic radius (in meters)
lat = yFromRow(r, 1:nrow(r)) # latitude of the centroid of each cell (in degrees, need to be converted in radians)
area = (sin(pi/180*(lat + RES[2]/2)) - sin(pi/180*(lat - RES[2]/2))) * (RES[1] * pi/180) * R^2
r_area <- setValues(r, rep(area, each=ncol(r))) # gives the area of each grid cell in meters 


files <- list.files(dir_demClasses, full.names = TRUE)

df <- data.frame()
result <- c()

for (i in 1:length(files)){

  a <- raster(files[i])
  a[is.na(a)] <- 0 
  
  b <- sum(getValues(a)) # total number of birds in the 10km x 10km raster
  
  
  ## reprojecting and rescaling
  a_proj <- projectRaster(a, r, method = "bilinear")
  a_proj2 <- a_proj * r_area / 100000000 # rescaling the values in each cell
  a_proj2[is.na(a_proj2)] <- 0 
  c <- sum(getValues(a_proj2)) # total number of birds in the 1.25 x 1.25 degree raster
  
  ## aggregating results into 5 x 5 degree grids (by summing values in blocks of 16 cells)
  a_proj3 <- aggregate(a_proj2, fact=4, fun = sum) 
  d <- sum(getValues(a_proj3)) # total number of birds in the 5 x 5 degree raster
  
  ## exporting results
  raster_name_1 <- gsub(dir_demClasses, "", files[i])
  raster_name_2 <- paste0(dir_5by5, "/", raster_name_1)
  writeRaster(a_proj3, filename=raster_name_2, format="GTiff", overwrite=TRUE)
  
  ## check results (checking that total number of predicted birds is unaffected by the projection and aggregation)
  result <- cbind(b,c,d)        
  df <- rbind(df, as.data.frame(result))
}
  
  
write.csv(df, paste0(dir_5by5, "/check-results_rasters.csv"))  

  
  

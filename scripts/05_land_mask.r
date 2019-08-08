## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Mapping the global distribution of seabird populations
## R script to clip the kernels to the land
## Ana Carneiro & Lizzie Pearmain
## Aug 2018
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



#####################################################
################ LOADING PACKAGES ###################
#####################################################

library(rgdal)
library(raster)



#####################################################
######### GENERAL DIRECTIONS AND FILES ##############
#####################################################

## GENERAL DIR
dir <- "C:/Users/ana.carneiro/Documents/GEF analysis/supplements"  ## copy and paste here your working directory

## PROJECTIONS
DgProj <- CRS("+proj=laea +lon_0=111.686552278329 +lat_0=-88.6447576802893 +ellps=WGS84") ## based on the null grid coords

## YOUR LANDMAsK
mask <- readOGR(dsn= paste0(dir, "/baselayer"), layer = "world-dissolved")   ## copy and paste here the direction of your worldmap (shapefile, WGS84)
mask_proj <- spTransform(mask, DgProj)   ## changing projection
mask_proj_pol <- as(mask_proj, "SpatialPolygons")   ## converting SpatialPolygonsDataFrame to SpatialPolygons

## DIRECTION TO YOUR UDs
dir_UDs <- "C:/Users/ana.carneiro/Documents/GEF analysis/supplements/scripts_results/04_combining_selecting_renaming/renaming"

## SAVING RESULTS
## Create a folder in your computer to save results  
UDs_clip_land <- "C:/Users/ana.carneiro/Documents/GEF analysis/supplements/scripts_results/05_land_mask/"



#####################################################
#################### LOADING FILES ##################
#####################################################

all_files <- list.files(dir_UDs, full.names = TRUE)


for (i in 1:length(all_files)){
  print(i)
  
  ## read in files
  rast <- raster(all_files[i])
  
  ## set to NA cells that overlap mask (land)
  rast_mask <- mask(rast, mask_proj_pol, inverse = TRUE)
  rast_mask[is.na(rast_mask)] <- 0 

  ## rescale to 1
  rast_mask_sum1 <- rast_mask/sum(getValues(rast_mask))

  ## saving results
  rast_name_final <- sub("C:/Users/ana.carneiro/Documents/GEF analysis/supplements/scripts_results/04_combining_selecting_renaming/renaming/", "", all_files[i])
  rast_name_final <- paste0(UDs_clip_land, rast_name_final)
  writeRaster(rast_mask_sum1, filename=rast_name_final, format="GTiff", overwrite = TRUE)
}



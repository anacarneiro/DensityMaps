## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Mapping the global distribution of seabird populations
## R script to run kernel analysis (per individual and then merged) 
## Ana Carneiro
## May 2018
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



#####################################################
################ LOADING PACKAGES ###################
#####################################################

lu=function (x=x) length(unique(x))
library(rgeos)
library(rgdal)
library(sp)
library(geosphere)
library(adehabitatHR)
library(raster)


#####################################################
######### GENERAL DIRECTIONS AND FILES ##############
#####################################################

## GENERAL DIR
dir <- "C:/Users/ana.carneiro/Documents/GEF analysis/supplements"  ## copy and paste here your working directory

## PROJECTIONS
land <- readOGR(dsn="C:/Users/ana.carneiro/Documents/GEF analysis/supplements/baselayer", layer = "world-dissolved")  ## copy and paste here the direction of your worldmap (shapefile, WGS84)
DgProj <- CRS("+proj=laea +lon_0=111.686552278329 +lat_0=-88.6447576802893 +ellps=WGS84") ## based on the null grid coords

## TO SAVE KERNEL RESULTS
## Create a folder in your computer to save kernel results (for breeding_stage_v1 and breeding_stage_v2)  
dir_kernels_v1 <- "C:/Users/ana.carneiro/Documents/GEF analysis/supplements/scripts_results/03_kernels/breed_stage_v1"  
dir_kernels_v2 <- "C:/Users/ana.carneiro/Documents/GEF analysis/supplements/scripts_results/03_kernels/breed_stage_v2"


#####################################################
################# LOADING SPP DATA ##################
#####################################################

species <- "Tristan Albatross"

## read species file after cleaning (results from the script 02_cleaning_data)
df <- read.csv(paste0(dir, "/scripts_results/02_cleaning_data/", species, ".csv"))  
df$dtime <- as.POSIXlt(df$dtime, format = "%Y-%m-%d %H:%M:%S ")
str(df)


#######################################################
##### CREATING A NULL GRID (SOUTHERN HEMISPHERE) ######
#######################################################

so.grid <- expand.grid(LON = seq(-177.5, 177.5, 5), LAT = seq(-87.5, -2.5, 5))
coordinates(so.grid) <- ~LON+LAT
proj4string(so.grid) <- proj4string(land)
so.grid.proj <- spTransform(so.grid, CRS=DgProj)
coords <- so.grid.proj@coords

c <- min(coords[,1])   ## to check my min lon
d <- max(coords[,1])   ## to check my max lon

e <- min(coords[,2])   ## to check my min lat
f <- max(coords[,2])   ## to check my max lat

## by= 10000 means my cell size
a= seq(c, d, by=10000)
b= seq(e, f, by=10000)
null.grid <- expand.grid(x=a,y=b)
coordinates(null.grid) <- ~x+y
gridded(null.grid) <- TRUE
class(null.grid)
plot(null.grid)


#######################################################
################## KERNEL ANALYSIS ####################
#######################################################

all_data <- df
head(all_data)
str(all_data)



## BASED ON BREEDING STAGE V1
## per individual


## this is based on metadata created on 2_cleaning_data script
meta <- read.csv(paste0(dir, "/metadata_files/02_species_to_run_kernels/metadata_", species, "_breed_stage_v1.csv")) 

all_data$breed_stage_v1 <- factor(all_data$breed_stage_v1)
all_data$bird_id <- factor(all_data$bird_id)
all_data$colony_name <- factor(all_data$colony_name)
all_data$device_comb <- factor(all_data$device_comb)
all_data$age <- factor(all_data$age)


for (i in 1:nrow(meta)){
  print(i)
  tracks <- all_data[(all_data$common_name==meta$common_name[i]) & (all_data$colony_name==meta$colony_name[i]) & (all_data$device_comb==meta$device_comb[i]) &
                       (all_data$age==meta$age[i]) & (all_data$breed_stage_v1==meta$breed_stage_v1[i]),]
  coordinates(tracks) <- ~coords_lon+coords_lat
  proj4string(tracks) <- DgProj
  tracks$bird_id <- factor(tracks@data$bird_id)
  tracks$device_comb <- factor(tracks$device_comb)
  
  if(tracks$device[1]=="GLS"){
    kerns_gls <- kernelUD(tracks[,11], grid = null.grid, h = 200000)  ## smoothing factor equals 200 km for GLS data
    stk_gls <- stack(estUDm2spixdf(kerns_gls))
    sum_all_gls_raw <- overlay(stk_gls, fun = mean)
    sum_all_gls_sum <- sum_all_gls_raw/sum(getValues(sum_all_gls_raw))
    
    KDERasName_gls_sum = paste(dir_kernels_v1, "/", meta$common_name[i], "_", meta$colony_name[i], "_", meta$device_comb[i], "_", 
                               meta$age[i], "_", meta$breed_stage_v1[i], "_sum", ".tif", sep="")
    writeRaster(sum_all_gls_sum, filename=KDERasName_gls_sum, format="GTiff", overwrite = TRUE)} else {
      
      kerns_gps_ptt <- kernelUD(tracks[,11], grid = null.grid, h = 50000)   ## smoothing factor equals 50 km for PTT and GPS data
      stk_gps_ptt <- stack(estUDm2spixdf(kerns_gps_ptt))
      sum_all_gps_ptt_raw <- overlay(stk_gps_ptt, fun = mean)
      sum_all_gps_ptt_sum <- sum_all_gps_ptt_raw/sum(getValues(sum_all_gps_ptt_raw))
      
      KDERasName_gps_ptt_sum = paste(dir_kernels_v1, "/", meta$common_name[i], "_", meta$colony_name[i], "_", meta$device_comb[i], "_", 
                                     meta$age[i], "_", meta$breed_stage_v1[i], "_sum", ".tif", sep="")
      writeRaster(sum_all_gps_ptt_sum, filename=KDERasName_gps_ptt_sum, format="GTiff", overwrite = TRUE)}
}



## BASED ON BREEDING STAGE V2

## this is based on metadata created on 2_cleaning_data script
meta <- read.csv(paste0(dir, "/metadata_files/02_species_to_run_kernels/metadata_", species, "_breed_stage_v2.csv"))   

meta$colony_name <- factor(meta$colony_name)
meta$device_comb <- factor(meta$device_comb)
meta$age <- factor(meta$age)
meta$breed_stage_v2 <- factor(meta$breed_stage_v2)

## subsetting tracking data to  include only data which breed_stage_v2 equals "non-breeding_Q2_Q3" and "non-breeding_Q4_Q1"
all_data_v2 <- subset(all_data, breed_stage_v2 %in% c("non-breeding_Q2_Q3", "non-breeding_Q4_Q1"))   
all_data <- all_data_v2
all_data$breed_stage_v2 <- factor(all_data$breed_stage_v2)
all_data$bird_id <- factor(all_data$bird_id)
all_data$colony_name <- factor(all_data$colony_name)
all_data$device_comb <- factor(all_data$device_comb)
all_data$age <- factor(all_data$age)


for (i in 1:nrow(meta)){
  print(i)
  tracks <- all_data[(all_data$common_name==meta$common_name[i]) & (all_data$colony_name==meta$colony_name[i]) & (all_data$device_comb==meta$device_comb[i]) &
                       (all_data$age==meta$age[i]) & (all_data$breed_stage_v2==meta$breed_stage_v2[i]),]
  coordinates(tracks) <- ~coords_lon+coords_lat
  proj4string(tracks) <- DgProj
  tracks$bird_id <- factor(tracks@data$bird_id)
  tracks$device_comb <- factor(tracks$device_comb)
  
  if(tracks$device[1]=="GLS"){
    kerns_gls <- kernelUD(tracks[,11], grid = null.grid, h = 200000)
    stk_gls <- stack(estUDm2spixdf(kerns_gls))
    sum_all_gls_raw <- overlay(stk_gls, fun = mean)
    sum_all_gls_sum <- sum_all_gls_raw/sum(getValues(sum_all_gls_raw))

    KDERasName_gls_sum = paste(dir_kernels_v2, "/", meta$common_name[i], "_", meta$colony_name[i], "_", meta$device_comb[i], "_", 
                               meta$age[i], "_", meta$breed_stage_v2[i], "_sum", ".tif", sep="")
    writeRaster(sum_all_gls_sum, filename=KDERasName_gls_sum, format="GTiff", overwrite = TRUE)} else {
      
      kerns_gps_ptt <- kernelUD(tracks[,11], grid = null.grid, h = 50000)
      stk_gps_ptt <- stack(estUDm2spixdf(kerns_gps_ptt))
      sum_all_gps_ptt_raw <- overlay(stk_gps_ptt, fun = mean)
      sum_all_gps_ptt_sum <- sum_all_gps_ptt_raw/sum(getValues(sum_all_gps_ptt_raw))

      KDERasName_gps_ptt_sum = paste(dir_kernels_v2, "/", meta$common_name[i], "_", meta$colony_name[i], "_", meta$device_comb[i], "_", 
                                     meta$age[i], "_", meta$breed_stage_v2[i], "_sum", ".tif", sep="")
      writeRaster(sum_all_gps_ptt_sum, filename=KDERasName_gps_ptt_sum, format="GTiff", overwrite = TRUE)}
}



## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Mapping the global distribution of seabird populations
## R script to clean and prepare files (STD) prior to kernel analysis 
## Ana Carneiro
## May 2018
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



#####################################################
################ LOADING PACKAGES ###################
#####################################################

library(maps)
library(geosphere)
library(trip)
library(adehabitatHR)
library(maptools)
library(rgdal)
library(sp)
library(rgeos)
library(geosphere)
library(mapdata)
library(tidyverse)
library(dplyr)
library(rworldmap)
data(countriesLow)
lu=function (x=x) length(unique(x))



#####################################################
######### GENERAL DIRECTIONS AND FILES ##############
#####################################################

## GENERAL DIR
dir <- "C:/Users/ana.carneiro/Documents/GEF analysis/supplements"  ## copy and paste here your working directory

## PROJECTIONS
land <- readOGR(dsn="C:/Users/ana.carneiro/Documents/GEF analysis/supplements/baselayer", layer = "world-dissolved")  ## copy and paste here the direction of your worldmap (shapefile, WGS84)
DgProj <- CRS("+proj=laea +lon_0=111.686552278329 +lat_0=-88.6447576802893 +ellps=WGS84") ## based on the null grid coords

## DIRECTION OF THE ORIGINAL SPECIES FILES (AS DOWNLOADED FROM THE SEABIRD TRACKING DATABASE)
data_std <- "C:/Users/ana.carneiro/Documents/GEF analysis/supplements/original_datasets"



#####################################################
################# LOADING SPP DATA ##################
#####################################################

## Analysis for this manuscript included data for all the species below
## For illustrative pusposes we are running the scripts with Tristan Albatross data only

#species <- "Amsterdam Albatross"
#species <- "Antipodean Albatross"
#species <- "Atlantic Yellow-nosed Albatross"
#species <- "Black Petrel"
#species <- "Black-browed Albatross"
#species <- "Buller's Albatross"
#species <- "Chatham Albatross"
#species <- "Grey Petrel"
#species <- "Grey-headed Albatross"
#species <- "Light-mantled Albatross"
#species <- "Northern Giant Petrel" 
#species <- "Northern Royal Albatross"
#species <- "Salvin's Albatross"
#species <- "Sooty Albatross"
#species <- "Southern Giant Petrel"
#species <- "Spectacled Petrel"
#species <- "Wandering Albatross"  
#species <- "Westland Petrel"
#species <- "White-capped Albatross"
#species <- "White-chinned Petrel"
#species <- "Indian Yellow-nosed Albatross"


species <- "Tristan Albatross"

df <- read.csv(paste0(data_std, "/", species, ".csv"))  



#####################################################
###### CHECKING DATA AND CREATING DTIME COLUMN ######
#####################################################

df$bird_id <- factor(df$bird_id)
df$track_id <- factor(df$track_id)
df$dtime <- as.POSIXct(strptime(paste(df$date_gmt, df$time_gmt, sep=" "), "%Y-%m-%d %H:%M:%S ")) 
head(df)
str(df)

## PLOTTING TO CHECK RESULTS 
xlow<-min(df$longitude)+0.2
xup<-max(df$longitude)-0.2
yup<-max(df$latitude)-0.5
ylow<-min(df$latitude)+0.5
windows(600,400)
plot(latitude~longitude, data=df, pch=16, cex=0.3, col=factor(dataset_id), asp=1, xlim=c(xlow,xup), ylim=c(ylow,yup), main="", frame=F, axes=F, xlab="", ylab="")
plot(countriesLow, col='darkgrey', add=T)
legend("topright", legend=levels(factor(df$dataset_id)), pch=16, col=unique(factor(df$dataset_id)))



#####################################################
################## DATA CLEANING ####################
############# SPEED FILTER FOR PTT DATA #############
#####################################################

## REMOVE LOCATIONS WITH > 90 KM/H

x2 <- data.frame()
x3 <- data.frame()
id <- unique(as.factor(df$track_id))

for(i in 1:nlevels(id)){
  a <- subset(df, track_id == id[i])
  if(a$device[1]=="PTT" & nrow(a) > 5){
    # Order the rows by ID, then by time
    b <- a[order(a$dtime), ]
    # Remove completely-duplicated rows
    b <- b[!duplicated(b),]
    b$dtime <- adjust.duplicateTimes(b$dtime, b$track_id)
    # Change times to hours since first fix
    b$hours <- as.numeric(difftime(b$dtime, min(b$dtime),units = "hours"))
    x2 <- rbind(x2, as.data.frame(b))} else {
      x3 <- rbind(x3, as.data.frame(a))
    }
}
 
x1 <- x2

## APPLY MCCONNELL SPEED FILTER IN TRIP PACKAGE TO REMOVE ERRONEOUS FIXES
x2 <- data.frame(lat = x1$latitude, lon = x1$longitude, DateTime = x1$dtime, id = x1$track_id)

## CREATE COORDINATE VARIABLE
coordinates(x2) <- c("lon","lat")

## CREATE TRIP OBJECT
tr <- trip(x2, c("DateTime","id"))

## MCCONNELL SPEED FILTER; ignore coordinates warning as data are lonlat
x1$Filter <- speedfilter(tr, max.speed = 90)

## REMOVE FILTERED COORDINATES
x1 <- subset(x1, x1$Filter==TRUE)
x1$hours <- NULL
x1$Filter <- NULL

## COMBINE PTT FILTERED FILE BACK INTO DF FILE (with GLS and GPS data)
df <- rbind(x1,x3)

## PLOT DATA TO CHECK IT LOOK OK
xlow<-min(df$longitude)+0.2
xup<-max(df$longitude)-0.2
yup<-max(df$latitude)-0.5
ylow<-min(df$latitude)+0.5
windows(600,400)
plot(latitude~longitude, data=df, pch=16, cex=0.3, col=factor(dataset_id), asp=1, xlim=c(xlow,xup), ylim=c(ylow,yup), main="", frame=F, axes=F, xlab="", ylab="")
plot(countriesLow, col='darkgrey', add=T)
legend("topright", legend=levels(factor(df$dataset_id)), pch=16, col=unique(factor(df$dataset_id)))



#####################################################
################ LINEAR INTERPOLATION ###############
################## PTT AND GPS DATA #################
#####################################################

x4 <- df

## BECAUSE THERE ARE DATASETS FROM NEW ZEALAND - BETTER TO PROJECT THE DATA BEFORE INTERPOLATION
coordinates(x4) <- ~longitude+latitude
proj4string(x4) <- CRS(proj4string(land))
## changing projection
x4_laea <- spTransform(x4, DgProj)
x4 <- as.data.frame(x4_laea)


## INTERPOLATION
x5 <- data.frame()
x6 <- data.frame()
dset_id <- unique(as.factor(x4$dataset_id))

x4$id_stage <- paste0(x4$track_id, "_", x4$breed_stage)   ## had to create that to run the analysis because some tracks start as breeding and end up as non-breeding

for(i in 1:nlevels(dset_id)){
  print(dset_id[i])
  tracks <- x4[x4$dataset_id==dset_id[i],]
  tracks$device <- factor(tracks$device)

    if(tracks$device[1]!="GLS"){
    tracks$id_stage <- factor(tracks$id_stage)
    tracks$track_time <- as.double(tracks$dtime)
    tracks$track_time <- adjust.duplicateTimes(tracks$track_time, tracks$id_stage)
    traj <- as.ltraj(xy=data.frame(tracks$longitude, tracks$latitude), date=as.POSIXct(tracks$track_time, origin="1970/01/01", tz="GMT"), id=tracks$id_stage, 
                     typeII = TRUE)
    

    ## Rediscretization every 3600 seconds
    tr <- adehabitatLT::redisltraj(traj, 3600, type="time")
    
    ## Convert output into a data frame
    tracks.intpol <- data.frame()
    for (l in 1:length(unique(tracks$id_stage))){
      # print(tracks$id_stage[l])
      out <- tr[[l]]
      out$ID <- as.character(attributes(tr[[l]])[4])				
      tracks.intpol <- rbind(tracks.intpol, out)} 
    
    ### re-insert
    
    tracks.intpol$dataset_id <- tracks$dataset_id[match(tracks.intpol$ID,tracks$id_stage)]
    tracks.intpol$scientific_name <- tracks$scientific_name[match(tracks.intpol$ID,tracks$id_stage)]
    tracks.intpol$common_name <- tracks$common_name[match(tracks.intpol$ID,tracks$id_stage)]
    tracks.intpol$site_name <- tracks$site_name[match(tracks.intpol$ID,tracks$id_stage)]
    tracks.intpol$colony_name <- tracks$colony_name[match(tracks.intpol$ID,tracks$id_stage)]
    tracks.intpol$lat_colony <- tracks$lat_colony[match(tracks.intpol$ID,tracks$id_stage)]
    tracks.intpol$lon_colony <- tracks$lon_colony[match(tracks.intpol$ID,tracks$id_stage)]
    tracks.intpol$device <- tracks$device[match(tracks.intpol$ID,tracks$id_stage)]
    tracks.intpol$bird_id <- tracks$bird_id[match(tracks.intpol$ID,tracks$id_stage)]
    tracks.intpol$track_id <- tracks$track_id[match(tracks.intpol$ID,tracks$id_stage)]
    tracks.intpol$original_track_id <- tracks$original_track_id[match(tracks.intpol$ID,tracks$id_stage)]
    tracks.intpol$age<-tracks$age[match(tracks.intpol$ID,tracks$id_stage)]
    tracks.intpol$sex<-tracks$sex[match(tracks.intpol$ID,tracks$id_stage)]
    tracks.intpol$breed_stage <- tracks$breed_stage[match(tracks.intpol$ID,tracks$id_stage)]
    tracks.intpol$breed_status <- tracks$breed_status[match(tracks.intpol$ID,tracks$id_stage)]
    tracks.intpol$track_time <- as.double(tracks.intpol$date)
   
    tracks.intpol <- tracks.intpol[order(tracks.intpol$ID, tracks.intpol$track_time),]
    tracks.intpol$track_time <- adjust.duplicateTimes(tracks.intpol$track_time, tracks.intpol$ID)
    
    #### combines all data
    x5 <- rbind(x5, as.data.frame(tracks.intpol))} else {
      x6 <- rbind(x6, as.data.frame(tracks))}
}
 
## renaming x5 and combining x5 and x6
head(x5)
head(x6) 
 
names(x5)[names(x5) == "date"] <- "dtime"
names(x5)[names(x5) == "x"] <- "longitude"
names(x5)[names(x5) == "y"] <- "latitude"
x5$ID <- NULL
x5$dx <- NULL
x5$dy <- NULL
x5$dist <- NULL
x5$dt <- NULL
x5$R2n <- NULL
x5$abs.angle <- NULL
x5$rel.angle <- NULL

x6$argos_quality <- NULL
x6$date_gmt <- NULL
x6$time_gmt <- NULL
x6$track_time <- as.double(x6$dtime)
x6$id_stage <- NULL

df <- rbind(x5,x6)


## PLOTTING TO CHECK RESULTS
xlow<-min(df$longitude)+0.2
xup<-max(df$longitude)-0.2
yup<-max(df$latitude)-0.5
ylow<-min(df$latitude)+0.5
windows(600,400)
plot(latitude~longitude, data=df, pch=16, cex=0.3, col=factor(dataset_id), asp=1, xlim=c(xlow,xup), ylim=c(ylow,yup), main="", frame=F, axes=F, xlab="", ylab="")
legend("topright", legend=levels(factor(df$dataset_id)), pch=16, col=unique(factor(df$dataset_id)))



#####################################################
############### REMOVE INLAND POSITION ##############
#####################################################

## Exclude all locations within 5 km (GPS) or 15 km (PTT) of the colony, and retain all the rest. Colony here is a single lat/lon position.

x7 <- df

x7$colony_name <- factor(x7$colony_name)
coordinates(x7) <- ~longitude+latitude
proj4string(x7) <- DgProj
land_p <-  spTransform(land, DgProj)

## Creating a buffer around colonies
gps_df <- data.frame()
ptt_df <- data.frame()
gls_df <- data.frame()

col <- unique(x7$colony_name)

for(i in 1:nlevels(col)){
  print(col[i])
  sub_col <- x7[x7$colony_name==col[i],]
  df_col <- data.frame(cbind(lon=sub_col$lon_colony[1], lat=sub_col$lat_colony[1]))
  coordinates(df_col) <- ~lon+lat
  ## assigning projection
  proj4string(df_col) <- CRS(proj4string(land))
  ## changing projection
  col_laea <- spTransform(df_col, DgProj)
  ## buffers
  buf_small <- gBuffer(col_laea, width = 5000)
  buf_large <- gBuffer(col_laea, width = 15000)
  ## removing land overlap
  sub_col@data$device <- factor(sub_col@data$device)
  dev <- unique(sub_col@data$device)
  for (j in 1:nlevels(dev)){
    print(dev[j])
    a <- sub_col[sub_col$device==dev[j],]
    if (a$device[1]=="GPS"){
      gps <- a[is.na(over(a, buf_small)),]
      gps_df <- rbind(gps_df, as.data.frame(gps))
      } else if (a$device[1]=="PTT"){
        ptt <- a[is.na(over(a, buf_large)),]
        ptt_df <- rbind(ptt_df, as.data.frame(ptt))
      } else {
          gls_df <- rbind(gls_df, as.data.frame(a))
      }
  }
}
 
df_all <- rbind(gps_df, ptt_df, gls_df)
df_all$coords_lon <- df_all$longitude
df_all$coords_lat <- df_all$latitude
coordinates(df_all) <- ~longitude+latitude
proj4string(df_all) <- CRS(proj4string(col_laea))
## changing projection
df_all_wgs <- spTransform(df_all, CRS(proj4string(land)))
df_all_wgs <- as.data.frame(df_all_wgs)
df <- df_all_wgs


## PLOTTING TO CHECK RESULTS
xlow<-min(df$longitude)+0.2
xup<-max(df$longitude)-0.2
yup<-max(df$latitude)-0.5
ylow<-min(df$latitude)+0.5
windows(600,400)
plot(latitude~longitude, data=df, pch=16, cex=0.3, col=factor(dataset_id), asp=1, xlim=c(xlow,xup), ylim=c(ylow,yup), main="", frame=F, axes=F, xlab="", ylab="")
legend("topright", legend=levels(factor(df$dataset_id)), pch=16, col=unique(factor(df$dataset_id)))
plot(countriesLow, col='darkgrey', add=T)



#######################################################
################### YEAR QUARTERS #####################
# Q1(Jan-Mar), Q2(Apr-Jun), Q3(July-Sep), Q4(Oct-Dec) #
#######################################################
#######################################################

## THIS IS TO CREATE NON-BREEDING (NB) DISTRIBUTIONS BASED ON YEAR QUARTER
## AND SEASONAL (SUMMER: Q4_Q1, WINTER: Q2_Q3)

df$quarter <- quarters(df$dtime)
str(df)
df$quarter <- factor(df$quarter)
levels(df$quarter)


## PREPARING DISTRIBUTIONS FOR KERNEL ANALYSIS (PREPARING A COLUMN TO GUIDE DATA GROUP ANALYSIS)
## CREATING NB DISTRIBUTIONS PER QUARTER (v1) AND SEASON (v2)
df$breed_stage_v1 <-  df$breed_stage
df$breed_stage_v2 <-  df$breed_stage

## DATASETS WITH NO NB DATA WILL PRODUCE AN ERROR - IGNORE
x1 <- subset(df, breed_stage == "non-breeding")
x1$quarter <- factor(x1$quarter)
x2 <- subset(df, !breed_stage == "non-breeding")
q <- unique(factor(x1$quarter))
A_summer <- data.frame()
A_winter <- data.frame()

for(i in 1:nlevels(q)){
  print(q)
  a <- subset(x1, quarter == q[i])
  a$quarter <- factor(a$quarter)
  if(a$quarter[1] %in% c("Q1", "Q4")){
    a$breed_stage_v1 <- paste0(a$breed_stage[1], "_", a$quarter[1])
    a$breed_stage_v2 <- paste0(a$breed_stage[1], "_", "Q4_Q1")
    A_summer <- rbind(A_summer, as.data.frame(a))} else {
      a$breed_stage_v1 <- paste0(a$breed_stage[1], "_", a$quarter[1])
      a$breed_stage_v2 <- paste0(a$breed_stage[1], "_", "Q2_Q3")
      A_winter <- rbind(A_winter, as.data.frame(a))}
}

NB <- rbind(A_summer, A_winter)        
df <- rbind(NB, x2)



#######################################################
############## SPECIES METADATA FOR KERNELS ###########
#######################################################

## Create a folder in your computer to export the metadata to run kernels
df$device_comb <- df$device
levels(df$device_comb)[levels(df$device_comb)=="GPS"] <- "GPSorPTT"
levels(df$device_comb)[levels(df$device_comb)=="PTT"] <- "GPSorPTT"


## metadata based on non-breeding stages split per year quarter
num_birds_v1 <- aggregate(df$bird_id, list(df$common_name, df$site_name, df$colony_name, df$device_comb, df$age, df$breed_stage_v1), lu)
colnames(num_birds_v1) <- c("common_name", "site_name", "colony_name", "device_comb", "age", "breed_stage_v1", "n_birds")
num_birds_v1

## metadata based on non-breeding stages split per seasonal combinations
num_birds_v2 <- aggregate(df$bird_id, list(df$common_name, df$site_name, df$colony_name, df$device_comb, df$age, df$breed_stage_v2), lu)
colnames(num_birds_v2) <- c("common_name", "site_name", "colony_name", "device_comb", "age", "breed_stage_v2", "n_birds")
num_birds_v2_sub <- subset(num_birds_v2, grepl("*non-breeding", breed_stage_v2))
num_birds_v2_sub

write.csv(num_birds_v1, paste0(dir, "/metadata_files/2_species_to_run_kernels/", "metadata_", num_birds_v1$common_name[1], "_breed_stage_v1.csv"), row.names = FALSE)
write.csv(num_birds_v2_sub, paste0(dir, "/metadata_files/2_species_to_run_kernels/", "metadata_", num_birds_v2_sub$common_name[1], "_breed_stage_v2.csv"), row.names = FALSE)



#######################################################
######## REMOVING GROUPS WITH LESS THAN 5 ROWS ########
#######################################################

## THIS IS A REQUIREMENT TO RUN KERNELS IN THE ADEHABITATHR PACKAGE
## BASED ON BREEDING STAGE V1 (AS IS THE MORE LIMITING. V2 COMBINES TWO NB QUARTERS)

df$breed_stage_v1 <- factor(df$breed_stage_v1)
df$bird_id <- factor(df$bird_id)
df$age <- factor(df$age)
num_birds_v1$age <- factor(num_birds_v1$age)


## REMOVING GROUPS WITH LESS THAN 5 ROWS
file_to_run <- data.frame()

for (i in 1:nrow(num_birds_v1)){
  print(i)
  birds_row <- df[(df$common_name==num_birds_v1$common_name[i]) & (df$colony_name==num_birds_v1$colony_name[i]) & 
                          (df$device_comb==num_birds_v1$device_comb[i]) & (df$age==num_birds_v1$age[i]) & 
                          (df$breed_stage_v1==num_birds_v1$breed_stage_v1[i]),]
  birds_row$bird_id <- factor(birds_row$bird_id)
  id <- unique(birds_row$bird_id)
  for (b in 1:nlevels(id)){
    sub_bird <- birds_row[(birds_row$bird_id==id[b]),]
    if (nrow(sub_bird)>5){
      file_to_run <- rbind(file_to_run, as.data.frame(sub_bird))}
  }
}

df <- file_to_run
head(df)



#######################################################
########### REORDER COLUMNS (ALL SAME FORMAT) #########
#######################################################

x <- df
x1 <- x[c("dtime", "track_id", "dataset_id",	"scientific_name",	"common_name",	"site_name",
          "colony_name",	"lat_colony", "lon_colony",	"device",	"bird_id",	"original_track_id",
          "age",	"sex",	"breed_stage",	"breed_status",	"track_time",	"coords_lon",
          "coords_lat",	"longitude",	"latitude",	"quarter",	"breed_stage_v1",	"breed_stage_v2",
          "device_comb")]
df <- x1
str(df)



#######################################################
############ EXPORT RESULTS READY TO KERNELS #########
#######################################################

## Create a folder in your computer to export tracking data already cleaned and formatted
write.csv(df, paste0(dir, "/scripts_results/02_cleaning_data/", species, ".csv"), row.names = FALSE)



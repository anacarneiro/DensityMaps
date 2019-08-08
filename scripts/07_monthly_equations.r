## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Mapping the global distribution of seabird populations
## R script to combine UDs according to equations
## Ana Carneiro & Lizzie Pearmain
## June 2018
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#####################################################
################ LOADING PACKAGES ###################
#####################################################

kpacks <- c('plyr','dplyr', 'raster')
new.packs <- kpacks[!(kpacks %in% installed.packages()[,"Package"])]
if(length(new.packs)) install.packages(new.packs)
lapply(kpacks, require, character.only=T)
remove(kpacks, new.packs)


#####################################################
######### GENERAL DIRECTIONS AND FILES ##############
#####################################################

## GENERAL DIR
dir <- "C:/Users/ana.carneiro/Documents/GEF analysis/supplements"  ## copy and paste here your working directory

## DIRECTION TO YOUR METADATA FILES (metadata spreadsheet for raster weightings)
dir_metadata <- "C:/Users/ana.carneiro/Documents/GEF analysis/supplements/scripts_results/06_combining_metadata_files"

## DIRECTION TO YOUR UDs 
dir_UDs <- "C:/Users/ana.carneiro/Documents/GEF analysis/supplements/scripts_results/05_land_mask"

## DIRECTION TO YOUR RESULTS
dir_monthly <- "C:/Users/ana.carneiro/Documents/GEF analysis/supplements/scripts_results/07_monthly_equations"


  
#####################################################
############### MONTHLY EQUATIONS ###################
#####################################################

# Read in metadata for the raster weightings 
# metadata is available for all species
meta_all <- read.csv(paste0(dir_metadata, "/metadata_monthly_equations.csv"))
head(meta_all)

# UDs available for a test case is for the Tristan Albatross
# Subsetting the metadata to run scripts
# BUT loop already created for when having multiple populations
meta_sub <- subset(meta_all, Species %in% "Tristan Albatross")
rownames(meta_sub) <- seq(length=nrow(meta_sub))
meta_all <- meta_sub


WEIGHTS <- data.frame()

meta_all$Species <- factor(meta_all$Species)  ## [change Species to factor so we can loop over it]
list_1_species <- levels(meta_all$Species)
for (i in 1:length(levels(meta_all$Species))){
  
  ## loop over each species -> meta_s1
  tmp_1_species <- levels(meta_all$Species)[i] ## keep track of which species we're on
  meta_s1 <- subset(meta_all, Species == tmp_1_species) ## subset meta to this species
  meta_s1$IslandGroup <- factor(meta_s1$IslandGroup)  ## [change IslandGroup to factor so we can loop over it]
  list_2_islandgroups <- levels(meta_s1$IslandGroup) ## what colonies are there for this species?
  
  ## loop over each IslandGroup -> meta_s2
  for(j in 1:length(levels(meta_s1$IslandGroup))){
    tmp_2_islandgroups <- levels(meta_s1$IslandGroup)[j] ## which IslandGroup we're on
    meta_s2 <- subset(meta_s1, IslandGroup == tmp_2_islandgroups) ## subset to this IslandGroup
    meta_s2$Month <- factor(meta_s2$Month) ## [change month to factor so we can loop over it]
    list_3_months <- levels(meta_s2$Month) ## what months are there for this IslandGroup?
    
    ## loop over each month -> meta_s3
    for(k in 1:length(levels(meta_s2$Month))){
      tmp_3_month <- levels(meta_s2$Month)[k] ## which month we're on
      meta_s3 <- subset(meta_s2, Month == tmp_3_month) ## subset to this month
      meta_s3$DemClass <- factor(meta_s3$DemClass)  ## [change DemClass to factor so we can loop over it]
      list_4_demclasses <- levels(meta_s3$DemClass) ## what demclasses are there in this month?
      
      ## loop over each DemClass -> meta_s4
      for(l in 1:length(levels(meta_s3$DemClass))){
        tmp_4_demclass <- levels(meta_s3$DemClass)[l]
        meta_s4 <- subset(meta_s3, DemClass == tmp_4_demclass)
        
        ## create weights -> wgt
        wgt <- data.frame()
        for(m in 1:nrow(meta_s4)){
          sub_row <-  meta_s4[m,]
          
          if(sub_row$HalfOnNest=="TRUE"){
            sub_row$Weight <- (sub_row$NumDays/sub_row$MonthLength) * 0.5
            sub_row$Check <- (sub_row$Weight * 2)
          } else if(sub_row$HalfOnNest=="FALSE"){
            sub_row$Weight <- (sub_row$NumDays/sub_row$MonthLength)
            sub_row$Check <- (sub_row$Weight)
          }
          
          ## Add in raster name to the wgt dataframe
          sub_row$Raster <- paste0(tmp_1_species, "_", tmp_2_islandgroups, "_", 
                                   sub_row$Device, "_", sub_row$RepAge, "_", sub_row$RepBreedStage, "_sum.tif")
          
          wgt <- rbind(wgt, as.data.frame(sub_row))
        }
        
        WEIGHTS <- rbind(WEIGHTS, as.data.frame(wgt))  ## to save results
        
          ## CHECKS ----
        ## Check that sum(NumDays) = MonthLength. - All days of the month are accounted for for that DemClass.
        stopifnot(isTRUE(all.equal(sum(wgt$NumDays), wgt$MonthLength[1])))
        ## Check that sum(Weights) = 1. Remember to multiply Weight by 2 IF HalfOnNest = TRUE, so that the weights will sum to 1. 
        stopifnot(isTRUE(all.equal(sum(wgt$Check), 1)))
        
        ## Read rasters from the wgt dataframe
        raster_stack <- stack()   
        for(n in 1:nrow(wgt)){
          sub_row_wgt <-  wgt[n,]
          
          read_rasters <- raster(paste0(dir_UDs, "/", sub_row_wgt$Raster))
          tmp_raster_weight <- read_rasters * sub_row_wgt$Weight     ## mutiply by the weight it represents
          raster_stack <- stack(raster_stack, tmp_raster_weight)}    ## combine results of each month for each demographic class
        sum_raster_stack <- sum(raster_stack)                        ## add rasters 
        times_demClass <- sum_raster_stack * wgt$value[1]            ## combine with the number of individuals
        
        ## Export results
        name_raster_1 <- paste0(dir_monthly, "/", 
                                tmp_1_species, "_", tmp_2_islandgroups, "_", tmp_3_month,  "_", tmp_4_demclass)
        
        writeRaster(times_demClass, filename=name_raster_1, format="GTiff", overwrite = TRUE)
      }
    }
  }
}
  
write.csv(WEIGHTS, paste0(dir_monthly, "/07_monthly_weights_to_check.csv"), row.names=FALSE)  


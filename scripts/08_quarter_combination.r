## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Mapping the global distribution of seabird populations
## R script to combine (average) rasters according to year quarters
## Ana Carneiro & Lizzie Pearmain
## June 2018
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



#####################################################
################ LOADING PACKAGES ###################
#####################################################

library(raster)



#####################################################
######### GENERAL DIRECTIONS AND FILES ##############
#####################################################

## GENERAL DIR
dir <- "C:/Users/ana.carneiro/Documents/GEF analysis/supplements"  ## copy and paste here your working directory

## DIRECTION TO YOUR METADATA FILE (metadata spreadsheet for combining months into year quarters)
dir_metadata <- "C:/Users/ana.carneiro/Documents/GEF analysis/supplements/scripts_results/07_monthly_equations"

## DIRECTION TO YOUR RASTERS (PER MONTH AND DEMCLASS)
dir_RASTER <- "C:/Users/ana.carneiro/Documents/GEF analysis/supplements/scripts_results/07_monthly_equations"

## DIRECTION TO YOUR RESULTS
dir_QUARTER <- "C:/Users/ana.carneiro/Documents/GEF analysis/supplements/scripts_results/08_quarter_combination"



#####################################################
########### AVERAGE RASTER PER YEAR QUARTER #########
#####################################################

# Read in metadata for the raster combination
# Test case is from Tristan Albatross
meta_all <- read.csv(paste0(dir_metadata, "/07_monthly_weights_to_check.csv"))
head(meta_all)


meta_all$Species <- factor(meta_all$Species) ## [change Species to factor so we can loop over it]
list_1_species <- levels(meta_all$Species)

for (i in 1:length(levels(meta_all$Species))){
  
  ## loop over each species -> meta_s1
  tmp_1_species <- levels(meta_all$Species)[i] ## keep track of which species we're on
  print(tmp_1_species)
  meta_s1 <- subset(meta_all, Species == tmp_1_species) ## subset meta to this species
  meta_s1$IslandGroup <- factor(meta_s1$IslandGroup)  ## [change IslandGroup to factor so we can loop over it]
  list_2_islandgroups <- levels(meta_s1$IslandGroup) ## what IslandGroup are there for this species?
  
  ## loop over each IslandGroup -> meta_s2
  for(j in 1:length(levels(meta_s1$IslandGroup))){
    tmp_2_islandgroups <- levels(meta_s1$IslandGroup)[j] ## which IslandGroup we're on
    print(tmp_2_islandgroups)
    meta_s2 <- subset(meta_s1, IslandGroup == tmp_2_islandgroups) ## subset to this IslandGroup
    meta_s2$DemClass <- factor(meta_s2$DemClass)  ## [change DemClass to factor so we can loop over it]
    list_3_demclasses <- levels(meta_s2$DemClass) ## what demclasses are there in this IslandGroup?
   
     ## loop over each DemClass -> meta_s3
    for(k in 1:length(levels(meta_s2$DemClass))){
      tmp_3_demclass <- levels(meta_s2$DemClass)[k]
      print(tmp_3_demclass)
      meta_s3 <- subset(meta_s2, DemClass == tmp_3_demclass)
      
      ## create patterns to match raster files (per year quarter)
      pattern_files_Q1 <- paste0(tmp_1_species, "_", tmp_2_islandgroups, "_", c("(1|2|3)"), "_", 
                                 tmp_3_demclass, ".tif")
      
      pattern_files_Q2 <- paste0(tmp_1_species, "_", tmp_2_islandgroups, "_", c("(4|5|6)"), "_", 
                                 tmp_3_demclass, ".tif")
      
      pattern_files_Q3 <- paste0(tmp_1_species, "_", tmp_2_islandgroups, "_", c("(7|8|9)"), "_", 
                                 tmp_3_demclass, ".tif")
      
      pattern_files_Q4 <- paste0(tmp_1_species, "_", tmp_2_islandgroups, "_", c("(10|11|12)"), "_", 
                                 tmp_3_demclass, ".tif")

      ## list files based on the previous patterns
      files_Q1 <- list.files(path = dir_RASTER, pattern = pattern_files_Q1, full.names = TRUE)
      files_Q2 <- list.files(path = dir_RASTER, pattern = pattern_files_Q2, full.names = TRUE)
      files_Q3 <- list.files(path = dir_RASTER, pattern = pattern_files_Q3, full.names = TRUE)
      files_Q4 <- list.files(path = dir_RASTER, pattern = pattern_files_Q4, full.names = TRUE)
      
      ## creating year quarter stacks
      stack_Q1 <- stack(files_Q1)
      stack_Q2 <- stack(files_Q2)
      stack_Q3 <- stack(files_Q3)
      stack_Q4 <- stack(files_Q4)
      
      ## averaging year quarters
      mean_Q1 <- overlay(stack_Q1, fun = mean)
      mean_Q2 <- overlay(stack_Q2, fun = mean)
      mean_Q3 <- overlay(stack_Q3, fun = mean)
      mean_Q4 <- overlay(stack_Q4, fun = mean)
      
      ## rasters names
      name_Q1 <- paste0(dir_QUARTER, "/", tmp_1_species, "_", 
                        tmp_2_islandgroups, "_", tmp_3_demclass, "_Q1")
      
      name_Q2 <- paste0(dir_QUARTER, "/", tmp_1_species, "_", 
                        tmp_2_islandgroups, "_", tmp_3_demclass, "_Q2")
      
      name_Q3 <- paste0(dir_QUARTER, "/", tmp_1_species, "_", 
                        tmp_2_islandgroups, "_", tmp_3_demclass, "_Q3")
      
      name_Q4 <- paste0(dir_QUARTER, "/", tmp_1_species, "_", 
                        tmp_2_islandgroups, "_", tmp_3_demclass, "_Q4")
      
      ## export results
      
      writeRaster(mean_Q1, filename=name_Q1, format="GTiff", overwrite = TRUE)
      writeRaster(mean_Q2, filename=name_Q2, format="GTiff", overwrite = TRUE)
      writeRaster(mean_Q3, filename=name_Q3, format="GTiff", overwrite = TRUE)
      writeRaster(mean_Q4, filename=name_Q4, format="GTiff", overwrite = TRUE)
    }
  }
}
      
        
 
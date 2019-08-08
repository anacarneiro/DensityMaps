## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Mapping the global distribution of seabird populations
## R script to create year-round maps
## Ana Carneiro and Lizzie Pearmain
## August 2018
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



#####################################################
################ LOADING PACKAGES ###################
#####################################################

library(raster)
library(rgdal)
library(dplyr)
library(tidyr)


#####################################################
######### GENERAL DIRECTIONS AND FILES ##############
#####################################################

## GENERAL DIR
dir <- "C:/Users/ana.carneiro/Documents/GEF analysis/supplements"  ## copy and paste here your working directory

## DIRECTION TO YOUR 5 x 5 grids
dir_5by5_quarter <- "C:/Users/ana.carneiro/Documents/GEF analysis/supplements/scripts_results/10_aggregate_5by5_grid"

## DIRECTION TO YOUR RESULTS
dir_5by5_year <- "C:/Users/ana.carneiro/Documents/GEF analysis/supplements/scripts_results/11_year_combination"



#####################################################
########### AVERAGE YEAR-ROUND DISTRIBUTIONS ########
#####################################################

# Create a metadata
files <- list.files(dir_5by5_quarter, pattern = ".tif$")
df_files <- data.frame(files)

df_files$files <- gsub(".tif", "", files, fixed = TRUE)
df_split1 <- df_files %>% separate(files, c("Species", "IslandGroup","Quarter"), sep = "_")
meta_all <- df_split1  ## renaming file

meta_all$Species <- factor(meta_all$Species)
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
 
    ## create patterns to match raster files 
    pattern_files_Qs <- paste0(tmp_1_species, "_", tmp_2_islandgroups, "_", c("(Q1|Q2|Q3|Q4)"), ".tif$")
    
    ## list files based on the previous patterns
    files_Qs <- list.files(path = dir_5by5_quarter, pattern = pattern_files_Qs, full.names = TRUE)
    
    ## creating year quarter stacks
    stack_Qs <- stack(files_Qs)
    
    ## averaging year quarters
    mean_Qs <- overlay(stack_Qs, fun = mean)
    
    ## rasters names
    name_Qs <- paste0(dir_5by5_year, "/", tmp_1_species, "_", 
                      tmp_2_islandgroups, "_year-round")
      
    ## export results
    writeRaster(mean_Qs, filename=name_Qs, format="GTiff", overwrite = TRUE)
  }
}



## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Mapping the global distribution of seabird populations
## R script to combine rasters according to demography classes
## Ana Carneiro & Lizzie Pearmain
## June 2018
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



#####################################################
################ LOADING PACKAGES ###################
#####################################################

library(raster)
library(dplyr)
library(tidyr)



#####################################################
######### GENERAL DIRECTIONS AND FILES ##############
#####################################################

## GENERAL DIR
dir <- "C:/Users/ana.carneiro/Documents/GEF analysis/supplements"  ## copy and paste here your working directory

## DIRECTION TO YOUR RASTERS (PER DEM CLASS AND YEAR QUARTER)
dir_RASTER <- "C:/Users/ana.carneiro/Documents/GEF analysis/supplements/scripts_results/08_quarter_combination"

## DIRECTION TO YOUR RESULTS
dir_demClasses <- "C:/Users/ana.carneiro/Documents/GEF analysis/supplements/scripts_results/09_sum_demClasses"



#####################################################
########### SUM DEM CLASSES PER YEAR QUARTER ########
#####################################################

# Create a metadata
files <- list.files(dir_RASTER)
df_files <- data.frame(files)

df_files$files <- gsub(".tif", "", files, fixed = TRUE)
df_split1 <- df_files %>% separate(files, c("Species", "IslandGroup","DemClass", "Quarter"), sep = "_")
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
    meta_s2$Quarter <- factor(meta_s2$Quarter)  ## [change quarter to factor so we can loop over it]
    list_3_Quarter <- levels(meta_s2$Quarter) ## what quarter are 

    ## loop over each Quarter -> meta_s3
    for(k in 1:length(levels(meta_s2$Quarter))){
      tmp_3_Quarter <- levels(meta_s2$Quarter)[k]
      print(tmp_3_Quarter)
      meta_s3 <- subset(meta_s2, Quarter == tmp_3_Quarter)
      
      ## create patterns to match raster files (per year quarter)
      pattern_files <- paste0(tmp_1_species, "_", tmp_2_islandgroups, "_", c("(beta|delta|gamma|theta|zeta)"), "_", 
                                 tmp_3_Quarter, ".tif")
      
      ## list files based on the previous patterns
      files <- list.files(path = dir_RASTER, pattern = pattern_files, full.names = TRUE)

      ## creating year quarter stacks
      stack_Quarter <- stack(files)

      ## sum demClasses per quarters
      sum_demClasses <- overlay(stack_Quarter, fun = sum)

      ## rasters names
      name_raster <- paste0(dir_demClasses, "/", tmp_1_species, "_", 
                            tmp_2_islandgroups, "_", tmp_3_Quarter)
      
      
      ## export results
      
      writeRaster(sum_demClasses, filename=name_raster, format="GTiff", overwrite = TRUE)
    }
  }
}

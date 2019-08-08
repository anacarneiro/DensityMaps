## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Mapping the global distribution of seabird populations
## R script to combine, select and rename UDs
## Ana Carneiro
## Aug 2018
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

## THE DIRECTION OF YOUR KERNELS 
dir_kernels_v1 <- "C:/Users/ana.carneiro/Documents/GEF analysis/supplements/scripts_results/03_kernels/breed_stage_v1/"  
dir_kernels_v2 <- "C:/Users/ana.carneiro/Documents/GEF analysis/supplements/scripts_results/03_kernels/breed_stage_v2/"



#####################################################
################# UD COMBINATION ####################
#####################################################

## When tracked from different breeding sites within the same island group:
# When seabirds were tracked from different breeding sites (i.e. colony) within the same island group, 
# UDs were combined by weighting the percentage of the total population involved. 
# Example: Tracking data available for the black-browed albatross for two breeding sites (within the Falklands)
# during brood-guard: 
# New Island and Steeple Jason. New Island has 13,331 and Steeple Jason has 145,964 breeding pairs according
# to ACAP species assessements, which represent 0.08 and 0.92, respectively. This represents, in coding:

# a <- raster(paste0(dir_kernels_v1, "Black-browed Albatross_New Island_GPSorPTT_adult_brood-guard_sum.tif"))      ## reading raster (UD)
# b <- raster(paste0(dir_kernels_v1, "Black-browed Albatross_Steeple Jason_GPSorPTT_adult_brood-guard_sum.tif"))   ## reading raster (UD)
# c <- ((a*0.08)+(b*0.92))   ## weighting according to the percentage
# rast_name_final <- paste0(dir, "/scripts_results/04_combining_selecting_renaming/combining/", "Black-browed Albatross_New Island_Steeple Jason_GPSorPTT_adult_brood-guard_sum")
# writeRaster(c, filename=rast_name_final, format="GTiff", overwrite = TRUE)


## When tracked using both PTT and/or GPS and GLS:
# If a data group had either fewer than 5 individuals or representativeness lower than 70% (with GPS and/or PTT data), 
# GPS and/or PTT data were combined with GLS data to increase sample sizes by weighting the UDs by the proportion of the 
# total number of individuals each sample represented. 
# Example: immature Black-browed Albatross GPS and/or PTT were combined with GLS data

a <- raster(paste0(dir_kernels_v2, "/",  "Black-browed Albatross_New Island_GLS_immature_non-breeding_Q4_Q1_sum.tif"))
b <- raster(paste0(dir_kernels_v2, "/",  "Black-browed Albatross_New Island_GPSorPTT_immature_non-breeding_Q4_Q1_sum.tif"))
c <- ((a*17)+(b*59))/76    
rast_name_final <- paste0(dir, "/scripts_results/04_combining_selecting_renaming/combining/", "Black-browed Albatross_New Island_GPSorPTTorGLS_immature_non-breeding_Q4_Q1_sum")
writeRaster(c, filename=rast_name_final, format="GTiff", overwrite = TRUE)



#####################################################
################## UD SELECTION #####################
#####################################################

# Prior to do the renaming you have to do some manual work. It may be that you will have many more UDs than the ones 
# you will be using for your next analysis. For example, you will have both the original UDs and the UDs derived from 
# the combinations above, and you will be using only the results from the combined ones. We recommend that 
# you create a new folder with the selected UDs containing only the UDs you will be working with from now on.
# Bootstrap step-by-step guide and codes can be found in Lascelles et al. (2016) and Oppel et al. (2018) 

# Our selected kernels for Tristan Albatross are here: 
# "C:/Users/ana.carneiro/Documents/GEF analysis/supplements/scripts_results/04_combining_selecting_renaming/selecting"



#####################################################
################### UD RENAMING #####################
#####################################################

# Breding sites and island/island group using the same definition as ACAP breeding site database, https://data.acap.aq/

# The UDs created for each breeding site (or the combination of breeding sites) were used as an extrapolation to the 
# wider island or island group population (and renamed accordinly).

# To rename your UDs you can follow the next steps:

# 1) List all your colony_name from the metadata that were used to create the UDs (your UDs were named referencing them)
# Here you will have just one species as example but the loop was created already for future use when you may have multiple files 

file_name <- list.files(paste0(dir, "/metadata_files/02_species_to_run_kernels/"), pattern = "*_v1.csv$")  ## there is no difference between the colonies in v1 or v2

all_colonies <- NULL

for(i in 1:length(file_name)){
  
  file <- read.csv(paste0(dir, "/metadata_files/02_species_to_run_kernels/", file_name[i]))
  colony <- unique(file$colony_name)
  all_colonies <- rbind(all_colonies, data.frame(colony))
}


# 2) Remove duplicates
all_colonies_unique <- all_colonies[!duplicated(all_colonies), ]


# 3) Create a dataframe with the colony_names and the corresponding island group 
df_colonies <- data.frame(all_colonies_unique)
colnames(df_colonies)[colnames(df_colonies)== "all_colonies_unique"] <- "colony"
df_colonies$IslandGroup <- df_colonies$colony
# Here you can list as many colonies you wish to rename
levels(df_colonies$IslandGroup)[levels(df_colonies$IslandGroup)=="Gough Island"] <- "Gough"


# 3) Rename UDs according to IslandGroup names

## loop over each dataframe row
for (i in 1:nrow(df_colonies)){
  
  sub_row <- df_colonies[i,]
  print(sub_row)

  # pattern to match all UDs within the row
  pattern_rast <- glob2rx(paste0("*_", sub_row$colony, "_*"))
  files_subset <- list.files(path = paste0(dir, "/scripts_results/04_combining_selecting_renaming/selecting"), pattern = pattern_rast, full.names = TRUE)
  
  if (length(files_subset) > 0){
    
    # loop over each UD
    for (j in 1:length(files_subset)){
      
      # read in file
      rast <- raster(files_subset[j])
      
      # rename according to island group
      rast_name_1 <- gsub(paste0("_", sub_row$colony, "_"), paste0("_", sub_row$IslandGroup, "_"), files_subset[j], fixed = TRUE)
      rast_name_final <- gsub("/selecting/", "/renaming/", rast_name_1, fixed = TRUE)
      
      # export results
      writeRaster(rast, filename=rast_name_final, format="GTiff", overwrite = TRUE)
    }
  }
}
  
    
    
    
    
  
  
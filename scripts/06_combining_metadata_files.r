## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Mapping the global distribution of seabird populations
## R script to combine phenology metadata files, and to link these data with the 
## number of individuals per life-history stage
## Ana Carneiro 
## Aug 2018
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



#####################################################
################ LOADING PACKAGES ###################
#####################################################

library("openxlsx")



#####################################################
######### GENERAL DIRECTIONS AND FILES ##############
#####################################################

## GENERAL DIR
dir <- "C:/Users/ana.carneiro/Documents/GEF analysis/supplements"  ## copy and paste here your working directory

## DIRECTION TO YOUR METADATA FILES (PHENOLOGY)
dir_phenology <- "C:/Users/ana.carneiro/Documents/GEF analysis/supplements/metadata_files/06_metadata_phenology/02_metadata_monthly_datatypes"


## DIRECTION TO YOUR METADATA FILE (DEMOGRAPHY)
dir_demography <- "C:/Users/ana.carneiro/Documents/GEF analysis/supplements/scripts_results/01_demography/"

## SAVING RESULTS
## Create a folder in your computer to save results  
metadata_monthly_equations <- "C:/Users/ana.carneiro/Documents/GEF analysis/supplements/scripts_results/06_combining_metadata_files/"



#####################################################
############# LOADING AND COMBINING FILES ###########
#####################################################

files <- list.files(path = dir_phenology, pattern = "*.xlsx", full.names = T)
all_files <- data.frame()

for (i in 1:length(files)){
  read.file <- read.xlsx(files[i], sheet=1)
  all_files <- rbind(all_files, as.data.frame(read.file))
}



#####################################################
### COMBINE METADATA FOR THE MONTHLY DISTRIBUTION ###
######## WITH NUMBER OF INDIVIDUALS PER CLASS #######
#####################################################

# Read in metadata for the number of individuals per life-history stage
# Results are from demography models

num_individuals <- read.csv(paste0(dir_demography, "num_individuals_per_class_2_transposed.csv"))
num_individuals$TaxonOrder <-  NULL

# Rename column to match phenology file
colnames(num_individuals)[colnames(num_individuals)== "CommonName"] <- "Species"
# Simplify Island Groups to match the phenology
levels(num_individuals$IslandGroup)[levels(num_individuals$IslandGroup)=="New Zealand (Great Barrier Island)"] <- "New Zealand"
levels(num_individuals$IslandGroup)[levels(num_individuals$IslandGroup)=="South Georgia (Islas Georgias del Sur)"] <- "South Georgia"
levels(num_individuals$IslandGroup)[levels(num_individuals$IslandGroup)=="Falkland Islands (Islas Malvinas)"] <- "Falkland Islands"

# Merge datasets
data_merge <- merge(all_files, num_individuals, by=c("Species","IslandGroup","DemClass" ))

# Export output table
write.csv(data_merge, paste0(metadata_monthly_equations, "metadata_monthly_equations.csv"), row.names=FALSE)

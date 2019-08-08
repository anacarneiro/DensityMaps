## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Mapping the global distribution of seabird populations
## R script to assess population proportions in each life-history stage 
## Steffen Oppel (steffen.oppel@rspb.org.uk), January 2018
## Updated by Lizzie Pearmain, April 2019
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



#####################################################
################ LOADING PACKAGES ###################
#####################################################

rm(list=ls())

library(popbio)
library(ggplot2)
library(doParallel)
library(foreach)
library(reshape2)
library(dplyr)



#####################################################
######### GENERAL DIRECTIONS AND FILES ##############
#####################################################

## GENERAL DIR
dir <- "C:/Users/ana.carneiro/Documents/GEF analysis/supplements"  ## copy and paste here your working directory

## TO SAVE RESULTS
## Create a folder in your computer to save results 
results <- "C:/Users/ana.carneiro/Documents/GEF analysis/supplements/scripts_results/01_demography"  



#####################################################
########## READ IN DEMOGRAPHY METADATA ##############
#####################################################

ddem <- read.csv(paste0(dir, "/metadata_files/01_demography/metadata_demography_to_run_simulation.csv"))
tmp <- sapply(ddem, is.factor)
ddem[tmp] <- lapply(ddem[tmp], as.character)  ## Convert fields of class factor into class character



#####################################################
############# START LOOP OVER EACH SPECIES ##########
#####################################################

DEM_RESULTS <- data.frame()

for (i in 1:nrow(ddem)){
  
  # i <- 1 ## testing
  ddem_spp <- ddem[i,]
  
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  #### a. DEFINE THE PARAMETERS FOR THIS SPECIES ----
  pop.size <- as.numeric(as.character(ddem_spp$PopEstimate))              ### population size in breeding pairs
  Sa <- as.numeric(as.character(ddem_spp$AnnualAdultSurvProb))            ### survival of adults
  Sj <- as.numeric(as.character(ddem_spp$AnnualJuveSurvProb))             ### survival of juveniles (up to a.f.b)
  BS <- as.numeric(as.character(ddem_spp$BreedingSuccess))                ### breeding success
  BF <- as.numeric(as.character(ddem_spp$BreedingFrequency))              ### breeding frequency
  afb <- round(as.numeric(as.character(ddem_spp$AgeFirstBreeding)), 0)    ### age at first breeding AS AN INTEGER (i.e. round to 0 decimal places)
  
  #### b. SET UP TABLE OF DEMOGRAPHY PARAMETERS ----
  ### if testing multiple params, use this table to loop over
  simul_in <- expand.grid(pop.size, Sa, Sj, BS, BF, afb)
  dim(simul_in)
  names(simul_in) <- c('pop.size','Sa','Sj','BS','BF','afb')
  # simul_in ## view the table to check
  SIM_OUT <- data.frame() ## create results data frame for simulating stable state for this species
  
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  #### c. DEFINE THE MATRIX FOR THIS SPECIES ----
  ## Model follows Abraham et al. 2016 - see diagram in manuscript.
  ## Age classes are Adult, age 1, age 2, ..., age a.f.b - 1.
  ## Modelling a population of females (= half of the total species population)
  
  dim1 <- afb      ## number of age classes = a.f.b.
  dim2 <- dim1^2   ## dimensions of matrix are a.f.b. * a.f.b.
  
  species.matrix <- (rep.int(0, dim2))                        ## construct Leslie matrix, filled with 0s
  species.matrix[0*dim1 + 1] <- expression(Sa)                ## row 1  , col 1    : adults in year t surviving to be adults in year t+1
  species.matrix[1*dim1 + 1] <- expression(BF*BS*0.5)         ## row 2  , col 1    : production of first year juveniles in year t+1 from adults in year t
  species.matrix[0*dim1 + dim1] <- expression(Sj)             ## row 1  , col dim1 (last column): immatures of age a.f.b.-1 in year t surviving to become adults in year t+1
  
  for (k in 2:(dim1-1)){
    species.matrix[k*dim1 + k] <- expression(Sj)              ## row k+1, col k    : fill in the off-diagonal Sj terms - juveniles/immatures surviving from year t to t+1
  }
  
  # species.matrix ## View the matrix to check (in the form of a 1-dimensional list for now)

  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  #### d. CALCULATING STABLE AGE DISTRIBUTION ----
    
  #### Create Leslie matrix using the species' vital rates
  seabird.vr <- list(Sa=simul_in[1,2], Sj=simul_in[1,3],
                     BS=simul_in[1,4], BF=simul_in[1,5],
                     afb=simul_in[1,6]) ## extract parameters for this simulation from the table of input params
   
  A <- matrix(sapply(species.matrix, eval, seabird.vr, NULL), nrow=sqrt(length(species.matrix)), byrow=TRUE) ## create Leslie matrix
  x <- stable.stage(A) ## find stable stage distribution
  # A ## View matrix to check
  # x ## View stable stage results
  
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  #### e. SUM INTO THE DIFFERENT AGE CLASSES ----
  out <- data.frame(matrix(ncol=9, nrow=1))               ## create data frame to output the age classes
  colnames(out) <- c("TaxonOrder","CommonName","IslandGroup","pop.size","Sa","Sj","BS","BF","afb")
  out[1, 1:3] <- c(ddem_spp$TaxonOrder[1],ddem_spp$CommonName[1],ddem_spp$IslandGroup[1])            ## Fill in info on species / island group
  out[1, 4:9] <- c(as.numeric(simul_in[[1]]),as.numeric(simul_in[[2]]),as.numeric(simul_in[[3]]),
                   as.numeric(simul_in[[4]]),as.numeric(simul_in[[5]]),as.numeric(simul_in[[6]]))    ## Fill in in demography params used
                
  out
  str(out)
  
  ### e.1: keep as proportions:
  out$Prop_AdBreed <- x[1]*BF                             ## Calculate number of adults breeders from number of adults * breeding frequency
  out$Prop_AdNonbreed <- x[1]*(1-BF)                      ## Calculate number of adults breeders from number of adults * (1 - breeding frequency)
  out$Prop_Juvenile <- x[2]                               ## Number of juveniles = number of birds aged 1 year
  out$Prop_Immature <- sum(x[3:(afb)])                    ## Number of immatures = number of birds aged 2 to a.f.b.-1 inclusive (between juvenile and adult)
  # out                                                     ## View to check
  # sum(c(out$Prop_AdBreed, out$Prop_AdNonbreed,
  #       out$Prop_Juvenile, out$Prop_Immature))            ## Check that these proportions sum to 1
  
  ### e.2: include number of birds:
  ## NB. pop.size is in units of breeding pairs = number of breeding females = PropAdBreed
  out$TotalPopNum <- out$pop.size/out$Prop_AdBreed*2      ## total population = Number of breeding females / PropAdBreed. Multiply by 2 to get number of both sexes.
  out$AdBreed <- out$Prop_AdBreed*out$TotalPopNum         ## Calculate number of adult breeders from proportion of pop * total population number
  out$AdNonbreed <- out$Prop_AdNonbreed*out$TotalPopNum   ## Calculate number of adult non-breeders
  out$Juvenile <- out$Prop_Juvenile*out$TotalPopNum       ## Calculate number of juveniles
  out$Immature <- out$Prop_Immature*out$TotalPopNum       ## Calculate number of immatures
  # out                                                     ## View to check that out$AdBreed (number of adult breeders) = out$pop.size * 2
  # sum(c(out$AdBreed, out$AdNonbreed,
  #       out$Juvenile, out$Immature))                      ## Check that these values sum to out$TotalPopNum

  ### e.3. add into DEM_RESULTS data frame
  SIM_OUT <- out
  DEM_RESULTS <- rbind(DEM_RESULTS, as.data.frame(SIM_OUT))

}



#####################################################
############### CHECK RESULTS DATA FRAME ############
#####################################################

View(DEM_RESULTS)
str(DEM_RESULTS)

df.check <- DEM_RESULTS

df.check$check1 <- df.check$Prop_AdBreed+df.check$Prop_AdNonbreed+df.check$Prop_Juvenile+df.check$Prop_Immature
df.check$check2 <- df.check$AdBreed+df.check$AdNonbreed+df.check$Juvenile+df.check$Immature
str(df.check)

round(df.check$check1, 8) == 1                                     ## Check that the proportions sum to 1
round(df.check$check2, 8) == round(df.check$TotalPopNum, 8)     ## Check that the numbers of birds sum to the calculated TotalPopNum
round(df.check$AdBreed, 8) == 2 * round(df.check$pop.size, 8)   ## Check that the number of adult breeders = 2 * pop.size (which was the number of breeding pairs)

rm(df.check)



#####################################################
## CONVERT THE DEMOGRAPHY OUTPUT INTO NUMBER OF IND #
###### PER DEMOGRAPHY CLASS (ALPHA, BETA, ETC) ######
#####################################################

# According to:
# Pools of birds:	    # What they are:		                 # How they're calculated:	
# alpha	              # adult breeders		                 # (from model - as NumberOfAdults * BreedingFrequency)	
# beta	              # adult successful breeders		       # alpha * breeding success	
# gamma	              # adult fail breeders		             # alpha * (1 - breeding success)	
# delta	              # adult non-breeders / sabbatical		 # (from model, as NumberOfAdults * (1 - BreedingFrequency))	
# zeta	              # juveniles		                       # (from model, as number of first year birds)	
# theta	              # immatures		                       # (from model, as number of birds in all age classes from year 2 to the year before recruiting into breeding pop (a.f.b - 1))	

ALL_RESULTS <- DEM_RESULTS 

ALL_RESULTS$alpha <- ALL_RESULTS$AdBreed
ALL_RESULTS$beta <- ALL_RESULTS$AdBreed*ALL_RESULTS$BS        ## Calculate number of successful breeders from number of breeders * breeding success
ALL_RESULTS$gamma <- ALL_RESULTS$AdBreed*(1-ALL_RESULTS$BS)   ## Calculate number of unsuccessful (fail) breeders from number of breeders * (1 - breeding success)
ALL_RESULTS$delta <- ALL_RESULTS$AdNonbreed
ALL_RESULTS$zeta <- ALL_RESULTS$Juvenile
ALL_RESULTS$theta <- ALL_RESULTS$Immature



#####################################################
###### RE-FORMAT THE DATA FOR THE NEXT STEPS OF #####
################### DATA ANALYSIS ###################
#####################################################

### a. current format
ALL_RESULTS_clean <- select(ALL_RESULTS, TaxonOrder=TaxonOrder, CommonName=CommonName,
                            IslandGroup=IslandGroup, alpha=alpha, beta=beta, gamma=gamma,
                            delta=delta, zeta=zeta, theta=theta)
str(ALL_RESULTS_clean)
write.csv(ALL_RESULTS_clean, paste0(results, "/num_individuals_per_class_1_non_transpose.csv"), row.names = FALSE)

### b. new format (sort of transposed)
ALL_RESULTS_clean_transpose <- reshape2::melt(ALL_RESULTS_clean, id=(c("CommonName", "IslandGroup","TaxonOrder")))
str(ALL_RESULTS_clean_transpose)
unique(ALL_RESULTS_clean_transpose$IslandGroup)
colnames(ALL_RESULTS_clean_transpose)[colnames(ALL_RESULTS_clean_transpose)=="variable"] <- "DemClass" ## Change field name to DemClass

ALL_RESULTS_clean_transpose
write.csv(ALL_RESULTS_clean_transpose, paste0(results, "/num_individuals_per_class_2_transposed.csv"), row.names = FALSE)
 


#####################################################
################## PLOT FIGURE 1 ####################
#####################################################

library(tidyverse)

df.ori <- ALL_RESULTS_clean_transpose

df <- df.ori ## save the original data frame 'df.ori' for comparison
head(df)
str(df)
## convert factor columns to character columns
i <- sapply(df, is.factor)
df[i] <- lapply(df[i], as.character)
str(df)

## simplify the island group names for plotting
df$IslandGroup <- gsub("South Georgia \\(Islas Georgias del Sur\\)", "South Georgia", df$IslandGroup)
df$IslandGroup <- gsub("Falkland Islands \\(Islas Malvinas\\)", "Falkland Islands", df$IslandGroup)
df$IslandGroup <- gsub("New Zealand \\(Great Barrier Island\\)", "New Zealand", df$IslandGroup)
df$IslandGroup <- gsub("Amsterdam and St Paul", "Amsterdam & St Paul", df$IslandGroup)
df

head(df)
unique(df$DemClass)

df <- df[-which(df$DemClass == "alpha") , ]
df$DemClass <- gsub("beta", "Adult.b.success.breeders", df$DemClass)
df$DemClass <- gsub("gamma", "Adult.c.fail.breeders", df$DemClass)
df$DemClass <- gsub("delta", "Adult.a.non.breeders", df$DemClass)
df$DemClass <- gsub("zeta", "Juveniles", df$DemClass)
df$DemClass <- gsub("theta", "Immatures", df$DemClass)
head(df)

## make a new df for plotting
df$name_group <- paste0(df$CommonName, ", ", df$IslandGroup)
head(df)
tab <- select(df, SpeciesIsland = name_group,
              DemClass = DemClass, value = value, TaxonOrder = TaxonOrder)
head(tab)

## Make a df to use to order the plot.
df.ord <- select(tab, TaxonOrder=TaxonOrder, SpeciesIsland=SpeciesIsland)
df.ord <- unique(df.ord)
df.ord$TaxonOrder <- as.numeric(df.ord$TaxonOrder)
str(df.ord)


#### make the plot ----
# dat <- format
# head(dat)
# dat$Species <- NULL
# dat$Island.Group <- NULL
# head(dat)
# 
library(reshape)
library(ggplot2)
library(wesanderson)
library(scales)
# datm <- melt(cbind(dat, Common.name = rownames(dat)), id.vars = c('Common.name'))
# head(datm)

datm <- tab
str(datm)
datm$SpeciesIsland <- as.factor(datm$SpeciesIsland)
datm$SpeciesIsland
datm$SpeciesIsland <- factor(datm$SpeciesIsland, levels = datm$SpeciesIsland[order(df.ord$TaxonOrder)])
str(datm)

win.graph(12,8)
ggplot(datm, aes(x = SpeciesIsland , y = value, fill = DemClass)) + 
  geom_bar(position = "fill",stat = "identity") +
  theme_bw() + 
  scale_y_continuous(labels = percent_format()) +
  scale_fill_brewer(palette="Spectral") +
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5)) 

ggsave(paste0(results, "/plot_num_individuals_per_class.pdf"))













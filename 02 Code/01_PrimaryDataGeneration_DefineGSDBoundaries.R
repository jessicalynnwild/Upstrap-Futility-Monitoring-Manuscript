###############################################################
### file: 01_PrimaryDataGeneration_DefineGSDBoundaries      ###
### authors: Jess Wild                                      ###
### creation date:    07/20/2023                            ###
### latest edit date: 07/20/2023                            ###
### description: code to calculate GSD p-value boundaries   ###
###############################################################

###############################################################
### SETUP WORKSPACE                                         ###
###############################################################

#clear environment
rm(list = ls())

#set working directory
setwd('C:/Users/wildje/Desktop/Dissertation/Paper 1')

#import packages
library(tidyverse)
library(rpact) #for group sequential

###############################################################
###  CALCULATE GROUP SEQUENTIAL BOUNDARIES                  ###
###############################################################

#calculate p-value boundaries for pocock GSD, futility only
design_PO_FO <- getDesignGroupSequential(
  sided = 2, alpha = 0.05, beta = 0.2, #statistical test specifications
  informationRates = c(0.25, 0.5, 0.75, 1), #information fraction for each interim look
  typeOfDesign = "asUser", #defines how efficacy monitoring is done
  userAlphaSpending = c(0, 0, 0, 0.05), #defines when efficacy monitoring is done
  typeBetaSpending = "bsP", #defines how futility monitoring is done
  bindingFutility = TRUE
)

#calculate p-value boundaries for obrian-fleming GSD, futility only
design_OBF_FO <- getDesignGroupSequential(
  sided = 2, alpha = 0.05, beta = 0.2,
  informationRates = c(0.25, 0.5, 0.75, 1), 
  typeOfDesign = "asUser",
  userAlphaSpending = c(0, 0, 0, 0.05), 
  typeBetaSpending = "bsOF",
  bindingFutility = TRUE
)

#calculate p-value boundaries for pocock GSD, efficacy only
design_PO_EO <- getDesignGroupSequential(
  sided = 2, alpha = 0.05, beta = 0.2,
  informationRates = c(0.25, 0.5, 0.75, 1), 
  typeOfDesign = "asP")

#calculate p-value boundaries for obrian-fleming GSD, efficacy only  
design_OBF_EO <- getDesignGroupSequential(
  sided = 2, alpha = 0.05, beta = 0.2,
  informationRates = c(0.25, 0.5, 0.75, 1), 
  typeOfDesign = "asOF")

#calculate p-value boundaries for pocock GSD, futility and efficacy 
design_PO_FE <- getDesignGroupSequential(
  sided = 2, alpha = 0.05,
  informationRates = c(0.25, 0.5, 0.75, 1), 
  typeOfDesign = "asP",
  typeBetaSpending = 'bsP', 
  bindingFutility = TRUE
)

#calculate p-value boundaries for obrian-fleming GSD, futility and efficacy 
design_OBF_FE <- getDesignGroupSequential(
  sided = 2, alpha = 0.05,
  informationRates = c(0.25, 0.5, 0.75, 1), 
  typeOfDesign = "asOF",
  typeBetaSpending = 'bsOF', 
  bindingFutility = TRUE
)

###############################################################
###  CONSOLIDATE GROUP SEQUENTIAL BOUNDARIES                ###
###############################################################

#combining all GSD efficacy and futility bounds into a single reference dataset
interim_fraction <- c(0.25,0.50,0.75,1.00)
MonitoringType <- c('FO','EO','FE')
BoundType <- c('OBF','PO')
dat_GSDbounds <- expand.grid(interim_fraction = interim_fraction, 
                             MonitoringType = MonitoringType, 
                             BoundType = BoundType) %>%
  mutate(FutilityCriticalValue = c(c(design_OBF_FO$futilityBounds,design_OBF_FO$criticalValues[4]),
                                   c(design_OBF_EO$futilityBounds,design_OBF_EO$criticalValues[4]),
                                   c(design_OBF_FE$futilityBounds,design_OBF_FE$criticalValues[4]),
                                   c(design_PO_FO$futilityBounds,design_PO_FO$criticalValues[4]),
                                   c(design_PO_EO$futilityBounds,design_PO_EO$criticalValues[4]),
                                   c(design_PO_FE$futilityBounds,design_PO_FE$criticalValues[4])),
         EfficacyCriticalValue = c(design_OBF_FO$criticalValues,
                                   design_OBF_EO$criticalValues,
                                   design_OBF_FE$criticalValues,
                                   design_PO_FO$criticalValues,
                                   design_PO_EO$criticalValues,
                                   design_PO_FE$criticalValues)) %>%
  mutate(FutilityPValue = 2*(1- pnorm(FutilityCriticalValue)),
         EfficacyPValue = 2*(1- pnorm(EfficacyCriticalValue)))

###############################################################
### SAVE OUTPUT DATASETS                                    ###
###############################################################

save(dat_GSDbounds, file = './01 Data/Data Processed/dat_GSDbounds.RData')
#load('./01 Data/Data Processed/dat_GSDbounds.RData')

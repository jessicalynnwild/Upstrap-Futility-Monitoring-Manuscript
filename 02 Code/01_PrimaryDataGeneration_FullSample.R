###############################################################
### file: 01_PrimaryDataGeneration_FullSample               ###
### authors: Jess Wild                                      ###
### creation date:    07/18/2023                            ###
### latest edit date: 07/20/2023                            ###
### description: code to calculate full sample hypothesis   ###
###   rejection rate                                        ###
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

#load input datasets
load('./01 Data/Data Processed/dat_results.RData')
load('./01 Data/Data Processed/dat_setting.RData')

###############################################################
### CALCULATE FULL SAMPLE REJECTION RATE                    ###
###############################################################

#consolidate full sample p value and rejection by sample size and iteration
dat_fullsample_reject <- dat_results %>%
  select(P_Observed_Full, Setting, Iteration) %>%
  unique() %>%
  mutate(Reject_FS = ifelse(as.numeric(P_Observed_Full) < 0.05, yes = 1, no = 0)) %>%
  merge(dat_setting, by='Setting') %>%
  select(n,power,Iteration,P_Observed_Full,Reject_FS) %>%
  unique()

###############################################################
### SAVE OUTPUT DATASETS                                    ###
###############################################################

save(dat_fullsample_reject, file = './01 Data/Data Processed/dat_fullsample_reject.RData')
#load('./01 Data/Data Processed/dat_fullsample_reject.RData')
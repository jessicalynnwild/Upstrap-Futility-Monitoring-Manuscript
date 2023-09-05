###############################################################
### file: 01_PrimaryDataGeneration_GSDFullSample            ###
### authors: Jess Wild                                      ###
### creation date:    07/20/2023                            ###
### latest edit date: 07/21/2023                            ###
### description: code to calculate full sample hypothesis   ###
###   rejection rate for GSD boundaries                     ###
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
load('./01 Data/Data Processed/dat_AB_settings.RData')

###############################################################
### CALCULATE GSD FULL SAMPLE REJECTION RATE                ###
###############################################################

#consolidate full sample p value and rejection by sample size and iteration (using GSD full sample p values)
dat_fullsample_reject_AB <- dat_results %>%
  select(P_Observed_Full, Setting, Iteration) %>%
  unique() %>%
  merge(dat_setting %>% select(Setting,n,power), by = 'Setting') %>%
  merge(dat_AB_settings %>% filter(interim_fraction == 1.00) %>% select(n,MonitoringType,EfficacyP), by = c('n')) %>%
  select(-Setting) %>%
  unique() %>%
  mutate(Reject_FS = ifelse(as.numeric(P_Observed_Full) < EfficacyP, yes = 1, no = 0)) %>%
  select(n,power,MonitoringType,Iteration,P_Observed_Full,Reject_FS)

###############################################################
### SAVE OUTPUT DATASETS                                    ###
###############################################################

save(dat_fullsample_reject_AB, file = './01 Data/Data Processed/dat_fullsample_reject_AB.RData')
#load('./01 Data/Data Processed/dat_fullsample_reject_AB.RData')

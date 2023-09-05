###############################################################
### file: 01_PrimaryDataGeneration_PowerTIE                 ###
### authors: Jess Wild                                      ###
### creation date:    07/18/2023                            ###
### latest edit date: 07/18/2023                            ###
### description: code to calculate power and type I error   ###
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

###############################################################
### Calculate Power/Type I Error                            ###
###############################################################

#calculate power (proportion of iterations per setting with p<0.05) 
dat_pow_err <- dat_results %>%
  select(P_Observed_Full, Setting, Iteration) %>%
  unique() %>%
  mutate(Reject = ifelse(as.numeric(P_Observed_Full) < 0.05, yes = 1, no = 0)) %>%
  group_by(Setting) %>%
  mutate(Power = round(mean(Reject),digits=2)) %>%
  select(Setting, Power) %>%
  unique()

###############################################################
### SAVE OUTPUT DATASETS                                    ###
###############################################################

save(dat_pow_err, file = './01 Data/Data Processed/dat_pow_err.RData')
#load('./01 Data/Data Processed/dat_pow_err.RData')

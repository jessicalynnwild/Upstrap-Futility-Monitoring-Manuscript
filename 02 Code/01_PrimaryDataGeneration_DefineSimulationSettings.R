###############################################################
### file: 01_PrimaryDataGeneration_DeffineSimulationSettings###
### authors: Jess Wild                                      ###
### creation date:   07/18/2023                             ###
### latest edit date:07/18/2023                             ###
### description: code to define overall simulation design   ###
### and conditions                                          ###
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

###############################################################
### DEFINE SIMULATION SETTINGS                              ###
###############################################################

#set up list of simulation scenario variables
n <- c(40,160,600,2000) #total sample size
power <- c(0.05,0.5,0.8,0.95) #statistical power to detect an effect
accrual_type <- 'unif' #static input parameter
outcome_name <- 'outcome_constant_high' #static input parameter
interim_fraction <- c(0.25,0.5,0.75) #define stopping points by information fraction

#generate reference dataset with simulation setting details
dat_setting <- expand.grid(n = n, power = power, accrual_type = accrual_type, 
                           outcome_name = outcome_name, interim_fraction = interim_fraction) %>%
  mutate(Setting = seq(1,48)) %>%
  rowwise() %>%
  mutate(p2 = 0.6, #calculate effect size and relative risk based on power and sample size
         p1 = power.prop.test(n = n/2, p2=0.6, power=power, sig.level=0.05)$p1,
         RR = (power.prop.test(n = n/2, p2=0.6, power=power, sig.level=0.05)$p1)/0.6) %>%
  mutate(p1 = ifelse(power == 0.05, yes = 0.6, no = p1), #manually set null relative risk scenarios with type I error = 0.05
         RR = ifelse(power == 0.05, yes = 1, no = RR)) %>%
  select(Setting,n,interim_fraction,power,p1,p2,RR,outcome_name,accrual_type)

###############################################################
### SAVE OUTPUT DATASETS                                    ###
###############################################################

save(dat_setting, file = './01 Data/Data Processed/dat_setting.RData')
#load('./01 Data/Data Processed/dat_setting.RData')

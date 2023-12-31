###############################################################
### file: 01_PrimaryDataGeneration_GU_Upstrap               ###
### authors: Jess Wild                                      ###
### creation date:    07/20/2023                            ###
### latest edit date: 07/20/2023                            ###
### description: code to apply AU/CU upstrapping boundaries ###
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
load('./01 Data/Data Processed/dat_fullsample_reject_AB.RData')

###############################################################
### APPLY ALPHA/BETA SPENDING CRITERIA TO SIMULATION RESULTS###
###############################################################

#apply alpha and beta spending upstrapped stopping boundaries to simulated results
#NOTE: data processing split into multiple parts by sample size due to RAM limitations
dat_calibrated_AB_results <- rbind(
  #n = 40
  merge(dat_results %>% select(Setting,Iteration,P_Upstrap), 
        dat_setting %>% select(Setting,n,interim_fraction,power) , by = 'Setting') %>%
    filter(n == 40) %>%
    merge(dat_AB_settings %>% filter(interim_fraction != 1.00 & n == 40), by = c('n','interim_fraction')) %>%
    group_by(Setting,Iteration,MonitoringType) %>% 
    mutate(p_lessthan_fut = ifelse(is.na(FutilityP) == F & (as.numeric(P_Upstrap) < FutilityP), yes = 1, no = 0),
           p_lessthan_eff = ifelse(is.na(EfficacyP) == F & (as.numeric(P_Upstrap) < EfficacyP), yes = 1, no = 0)) %>%
    summarise(proportion_sim_fut = mean(p_lessthan_fut),
              proportion_sim_eff = mean(p_lessthan_eff)) %>%
    merge(dat_setting %>% select(Setting,n,interim_fraction), by = 'Setting') %>%
    merge(dat_AB_settings, by = c('n','interim_fraction','MonitoringType')) %>%
    mutate(stop_futility = ifelse(proportion_sim_fut < FutilityProp, yes = 1, no = 0),
           stop_efficacy = ifelse(proportion_sim_eff > EfficacyProp, yes = 1, no = 0)),
  #n = 160
  merge(dat_results %>% select(Setting,Iteration,P_Upstrap), 
        dat_setting %>% select(Setting,n,interim_fraction,power) , by = 'Setting') %>%
    filter(n == 160) %>%
    merge(dat_AB_settings %>% filter(interim_fraction != 1.00 & n == 160), by = c('n','interim_fraction')) %>%
    group_by(Setting,Iteration,MonitoringType) %>% 
    mutate(p_lessthan_fut = ifelse(is.na(FutilityP) == F & (as.numeric(P_Upstrap) < FutilityP), yes = 1, no = 0),
           p_lessthan_eff = ifelse(is.na(EfficacyP) == F & (as.numeric(P_Upstrap) < EfficacyP), yes = 1, no = 0)) %>%
    summarise(proportion_sim_fut = mean(p_lessthan_fut),
              proportion_sim_eff = mean(p_lessthan_eff)) %>%
    merge(dat_setting %>% select(Setting,n,interim_fraction), by = 'Setting') %>%
    merge(dat_AB_settings, by = c('n','interim_fraction','MonitoringType')) %>%
    mutate(stop_futility = ifelse(proportion_sim_fut < FutilityProp, yes = 1, no = 0),
           stop_efficacy = ifelse(proportion_sim_eff > EfficacyProp, yes = 1, no = 0)),
  #n = 600
  merge(dat_results %>% select(Setting,Iteration,P_Upstrap), 
        dat_setting %>% select(Setting,n,interim_fraction,power) , by = 'Setting') %>%
    filter(n == 600) %>%
    merge(dat_AB_settings %>% filter(interim_fraction != 1.00 & n == 600), by = c('n','interim_fraction')) %>%
    group_by(Setting,Iteration,MonitoringType) %>% 
    mutate(p_lessthan_fut = ifelse(is.na(FutilityP) == F & (as.numeric(P_Upstrap) < FutilityP), yes = 1, no = 0),
           p_lessthan_eff = ifelse(is.na(EfficacyP) == F & (as.numeric(P_Upstrap) < EfficacyP), yes = 1, no = 0)) %>%
    summarise(proportion_sim_fut = mean(p_lessthan_fut),
              proportion_sim_eff = mean(p_lessthan_eff)) %>%
    merge(dat_setting %>% select(Setting,n,interim_fraction), by = 'Setting') %>%
    merge(dat_AB_settings, by = c('n','interim_fraction','MonitoringType')) %>%
    mutate(stop_futility = ifelse(proportion_sim_fut < FutilityProp, yes = 1, no = 0),
           stop_efficacy = ifelse(proportion_sim_eff > EfficacyProp, yes = 1, no = 0)),
  #n = 2000
  merge(dat_results %>% select(Setting,Iteration,P_Upstrap), 
        dat_setting %>% select(Setting,n,interim_fraction,power) , by = 'Setting') %>%
    filter(n == 2000) %>%
    merge(dat_AB_settings %>% filter(interim_fraction != 1.00 & n == 2000), by = c('n','interim_fraction')) %>%
    group_by(Setting,Iteration,MonitoringType) %>% 
    mutate(p_lessthan_fut = ifelse(is.na(FutilityP) == F & (as.numeric(P_Upstrap) < FutilityP), yes = 1, no = 0),
           p_lessthan_eff = ifelse(is.na(EfficacyP) == F & (as.numeric(P_Upstrap) < EfficacyP), yes = 1, no = 0)) %>%
    summarise(proportion_sim_fut = mean(p_lessthan_fut),
              proportion_sim_eff = mean(p_lessthan_eff)) %>%
    merge(dat_setting %>% select(Setting,n,interim_fraction), by = 'Setting') %>%
    merge(dat_AB_settings, by = c('n','interim_fraction','MonitoringType')) %>%
    mutate(stop_futility = ifelse(proportion_sim_fut < FutilityProp, yes = 1, no = 0),
           stop_efficacy = ifelse(proportion_sim_eff > EfficacyProp, yes = 1, no = 0))
)

###############################################################
### CALCULATE EVALUATION METRICS                            ###
###############################################################

#use stopping boundary results to get expected sample size and full trial result statistics for upstrapping with alpha beta spending
dat_final_results_AB <- merge(dat_setting %>% select(Setting,n,interim_fraction,power), 
                              dat_calibrated_AB_results %>% select(Setting,Iteration,MonitoringType,stop_futility,stop_efficacy), 
                              by='Setting') %>%
  select(-Setting) %>%
  mutate(interim_fraction = factor(interim_fraction)) %>%
  mutate(stop_futility = ifelse(is.na(stop_futility) == F, yes = stop_futility, no = 0),
         stop_efficacy = ifelse(is.na(stop_efficacy) == F, yes = stop_efficacy, no = 0)) %>%
  reshape(idvar = c("n","Iteration","power","MonitoringType"), 
          timevar = "interim_fraction", 
          direction = "wide") %>%
  select(n,power,MonitoringType,Iteration,stop_futility.0.25,stop_futility.0.5,stop_futility.0.75,stop_efficacy.0.25,stop_efficacy.0.5,stop_efficacy.0.75) %>%
  mutate(ExpectedN_FO = ifelse(stop_futility.0.25==1, yes = 0.25*n, no = 
                                 ifelse(stop_futility.0.5==1, yes = 0.5*n, no = 
                                          ifelse(stop_futility.0.75==1, yes = 0.75*n, no = n))),
         ExpectedN_EO = ifelse(stop_efficacy.0.25==1, yes = 0.25*n, no = 
                                 ifelse(stop_efficacy.0.5==1, yes = 0.5*n, no = 
                                          ifelse(stop_efficacy.0.75==1, yes = 0.75*n, no = n))),
         ExpectedN_FE = ifelse(stop_efficacy.0.25==1|stop_futility.0.25==1, yes = 0.25*n, no = 
                                 ifelse(stop_efficacy.0.5==1|stop_futility.0.5==1, yes = 0.5*n, no = 
                                          ifelse(stop_efficacy.0.75==1|stop_futility.0.75==1, yes = 0.75*n, no = n))),
         Decision_FO = ifelse(ExpectedN_FO<n,yes = 1, no = 0),
         Decision_EO = ifelse(ExpectedN_EO<n,yes = 1, no = 0),
         Decision_FE = ifelse(ExpectedN_FE<n,yes = 1, no = 0),
         Decision_FE_F = ifelse(Decision_FE==1&ExpectedN_FO<ExpectedN_EO,yes = 1, no = 0),
         Decision_FE_E = ifelse(Decision_FE==1&ExpectedN_EO<ExpectedN_FO,yes = 1, no = 0),
         Decision_FO_25 = ifelse(ExpectedN_FO==0.25*n,yes = 1, no = 0),
         Decision_EO_25 = ifelse(ExpectedN_EO==0.25*n,yes = 1, no = 0),
         Decision_FE_25 = ifelse(ExpectedN_FE==0.25*n,yes = 1, no = 0),
         Decision_FO_50 = ifelse(ExpectedN_FO==0.5*n,yes = 1, no = 0),
         Decision_EO_50 = ifelse(ExpectedN_EO==0.5*n,yes = 1, no = 0),
         Decision_FE_50 = ifelse(ExpectedN_FE==0.5*n,yes = 1, no = 0),
         Decision_FO_75 = ifelse(ExpectedN_FO==0.75*n,yes = 1, no = 0),
         Decision_EO_75 = ifelse(ExpectedN_EO==0.75*n,yes = 1, no = 0),
         Decision_FE_75 = ifelse(ExpectedN_FE==0.75*n,yes = 1, no = 0)) %>%
  merge(dat_fullsample_reject_AB, by=c('n','power','MonitoringType','Iteration')) %>%
  mutate(RejectionRate_FO_FS = ifelse(Decision_FO==0&Reject_FS==1,yes=1,no=0),
         RejectionRate_EO_FS = ifelse(Decision_EO==1|Reject_FS==1,yes=1,no=0),
         RejectionRate_FE_FS = ifelse(Decision_FE_E==1|(Decision_FE_F==0&Reject_FS==1),yes=1,no=0)) %>%
  select(n,power,MonitoringType,Iteration,ExpectedN_FO,ExpectedN_EO,ExpectedN_FE,Decision_FO,Decision_EO,Decision_FE,
         Decision_FE_F,Decision_FE_E,Decision_FO_25,Decision_EO_25,Decision_FE_25,Decision_FO_50,
         Decision_EO_50,Decision_FE_50,Decision_FO_75,Decision_EO_75,Decision_FE_75,RejectionRate_FO_FS,
         RejectionRate_EO_FS,RejectionRate_FE_FS) %>%
  group_by(n,power,MonitoringType) %>%
  summarise(ExpectedN_FO_mean = mean(ExpectedN_FO),
            ExpectedN_FO_sd = sd(ExpectedN_FO),
            ExpectedN_EO_mean = mean(ExpectedN_EO),
            ExpectedN_EO_sd = sd(ExpectedN_EO),
            ExpectedN_FE_mean = mean(ExpectedN_FE),
            ExpectedN_FE_sd = sd(ExpectedN_FE),
            DecisionFO = mean(Decision_FO),
            DecisionEO = mean(Decision_EO),
            DecisionFE = mean(Decision_FE),
            DecisionFE_F = mean(Decision_FE_F),
            DecisionFE_E = mean(Decision_FE_E),
            DecisionFO_25 = mean(Decision_FO_25),
            DecisionEO_25 = mean(Decision_EO_25),
            DecisionFE_25 = mean(Decision_FE_25),
            DecisionFO_50 = mean(Decision_FO_50),
            DecisionEO_50 = mean(Decision_EO_50),
            DecisionFE_50 = mean(Decision_FE_50),
            DecisionFO_75 = mean(Decision_FO_75),
            DecisionEO_75 = mean(Decision_EO_75),
            DecisionFE_75 = mean(Decision_FE_75),
            RejectionRate_FO_FS = mean(RejectionRate_FO_FS),
            RejectionRate_EO_FS = mean(RejectionRate_EO_FS),
            RejectionRate_FE_FS = mean(RejectionRate_FE_FS)) %>%
  mutate(MonitoringType = factor(MonitoringType)) %>%
  pivot_wider(
    id_cols = c('n','power'),
    names_from = c('MonitoringType'),
    values_from = c('ExpectedN_FO_mean','ExpectedN_FO_sd','ExpectedN_EO_mean','ExpectedN_EO_sd',
                    'ExpectedN_FE_mean','ExpectedN_FE_sd','DecisionFO','DecisionEO','DecisionFE',
                    'DecisionFO_25','DecisionFO_50','DecisionFO_75','DecisionEO_25','DecisionEO_50',
                    'DecisionEO_75','DecisionFE_25','DecisionFE_50','DecisionFE_75',
                    'DecisionFE_F','DecisionFE_E','RejectionRate_FO_FS','RejectionRate_EO_FS',
                    'RejectionRate_FE_FS')
  ) %>%
  select(n,power,ExpectedN_FO_mean_FO,ExpectedN_FO_sd_FO,
         ExpectedN_EO_mean_EO,ExpectedN_EO_sd_EO,
         ExpectedN_FE_mean_FE,ExpectedN_FE_sd_FE,
         DecisionFO_FO,DecisionEO_EO,DecisionFE_FE,DecisionFE_F_FE,DecisionFE_E_FE,
         DecisionFO_25_FO,DecisionEO_25_EO,DecisionFE_25_FE,
         DecisionFO_50_FO,DecisionEO_50_EO,DecisionFE_50_FE,
         DecisionFO_75_FO,DecisionEO_75_EO,DecisionFE_75_FE,
         RejectionRate_FO_FS_FO,RejectionRate_EO_FS_EO,RejectionRate_FE_FS_FE) %>%
  rename('ExpectedN_FO_mean' = ExpectedN_FO_mean_FO,
         'ExpectedN_FO_sd' = ExpectedN_FO_sd_FO,
         'ExpectedN_EO_mean' = ExpectedN_EO_mean_EO,
         'ExpectedN_EO_sd' = ExpectedN_EO_sd_EO,
         'ExpectedN_FE_mean' = ExpectedN_FE_mean_FE,
         'ExpectedN_FE_sd' = ExpectedN_FE_sd_FE,
         'DecisionFO' = DecisionFO_FO,
         'DecisionEO' = DecisionEO_EO,
         'DecisionFE' = DecisionFE_FE,
         'DecisionFE_F' = DecisionFE_F_FE,
         'DecisionFE_E' = DecisionFE_E_FE,
         'DecisionFO_25' = DecisionFO_25_FO,
         'DecisionEO_25' = DecisionEO_25_EO,
         'DecisionFE_25' = DecisionFE_25_FE,
         'DecisionFO_50' = DecisionFO_50_FO,
         'DecisionEO_50' = DecisionEO_50_EO,
         'DecisionFE_50' = DecisionFE_50_FE,
         'DecisionFO_75' = DecisionFO_75_FO,
         'DecisionEO_75' = DecisionEO_75_EO,
         'DecisionFE_75' = DecisionFE_75_FE,
         'RejectionRate_FO_FS' = RejectionRate_FO_FS_FO,
         'RejectionRate_EO_FS' = RejectionRate_EO_FS_EO,
         'RejectionRate_FE_FS' = RejectionRate_FE_FS_FE)

###############################################################
### SAVE OUTPUT DATASETS                                    ###
###############################################################

save(dat_calibrated_AB_results, file = './01 Data/Data Processed/dat_calibrated_AB_results.RData')
#load('./01 Data/Data Processed/dat_calibrated_AB_results.RData')
save(dat_final_results_AB, file = './01 Data/Data Processed/dat_final_results_AB.RData')
#load('./01 Data/Data Processed/dat_final_results_AB.RData')

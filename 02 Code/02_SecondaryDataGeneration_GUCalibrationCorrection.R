###############################################################
### file: 02_SecondaryDataGeneration_GUCalibrationCorrection###
### authors: Jess Wild                                      ###
### creation date:    07/25/2023                            ###
### latest edit date: 07/26/2023                            ###
### description: code to apply multiple testing correction  ###
### to GU upstrapping calibration                           ###
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
load('./01 Data/Data Processed/dat_setting.RData')
load('./01 Data/Data Processed/dat_results.RData')
load('./01 Data/Data Processed/dat_AB_settings.RData')
load('./01 Data/Data Processed/dat_fullsample_reject_AB.RData')
load('./01 Data/Data Processed/dat_final_output.RData')

###############################################################
### DEFINE CORRECTED PROPORTION THRESHOLDS                  ###
###############################################################

#define calibration thresholds for the calibration corrections
dat_AB_setting_calibration_corrected <- dat_AB_settings %>%
  mutate(FutilityProp_multCorr1 = FutilityProp/(4+1),
         FutilityProp_multCorr2 = FutilityProp/(4),
         FutilityProp_multCorr3 = FutilityProp/(4-1),
         FutilityProp_multCorr4 = FutilityProp/(4-1.5),
         FutilityProp_multCorr5 = FutilityProp/(4-2),
         FutilityProp_expCorr1 = FutilityProp^0.95,
         FutilityProp_expCorr2 = FutilityProp^0.70,
         EfficacyProp_multCorr1 = EfficacyProp+(1-EfficacyProp)/(4+1),
         EfficacyProp_multCorr2 = EfficacyProp+(1-EfficacyProp)/(4),
         EfficacyProp_multCorr3 = EfficacyProp+(1-EfficacyProp)/(4-1),
         EfficacyProp_multCorr4 = EfficacyProp+(1-EfficacyProp)/(4-1.5),
         EfficacyProp_multCorr5 = EfficacyProp+(1-EfficacyProp)/(4-2),
         EfficacyProp_expCorr1 = EfficacyProp^(1-0.95),
         EfficacyProp_expCorr2 = EfficacyProp^(1-0.70))

###############################################################
### APPLY M1 CALIBRATION CORRECTIONS TO GU RESULTS          ###
###############################################################

#apply M1 calibration correction stopping boundaries to GU results
dat_ABcalibrated_corrected_M1_upstraped_results <- merge(dat_results %>% select(Setting,Iteration,P_Upstrap), 
        dat_setting %>% select(Setting,n,interim_fraction,power) , by = 'Setting') %>%
    merge(dat_AB_setting_calibration_corrected %>% filter(interim_fraction != 1.00), by = c('n','interim_fraction')) %>%
    group_by(Setting,Iteration,MonitoringType) %>% 
    mutate(p_lessthan_fut = ifelse(is.na(FutilityP) == F & (as.numeric(P_Upstrap) < FutilityP), yes = 1, no = 0),
           p_lessthan_eff = ifelse(is.na(EfficacyP) == F & (as.numeric(P_Upstrap) < EfficacyP), yes = 1, no = 0)) %>%
    summarise(proportion_sim_fut = mean(p_lessthan_fut),
              proportion_sim_eff = mean(p_lessthan_eff)) %>%
    merge(dat_setting %>% select(Setting,n,interim_fraction), by = 'Setting') %>%
    merge(dat_AB_setting_calibration_corrected, by = c('n','interim_fraction','MonitoringType')) %>%
    mutate(stop_futility = ifelse(proportion_sim_fut < FutilityProp_multCorr1, yes = 1, no = 0),
           stop_efficacy = ifelse(proportion_sim_eff > EfficacyProp_multCorr1, yes = 1, no = 0))

#use stopping boundary results to get expected sample size and full trial result statistics
dat_ABcalibrated_corrected_M1_upstraped_final_results <- dat_ABcalibrated_corrected_M1_upstraped_results %>%
    select(Setting,Iteration,MonitoringType,stop_futility,stop_efficacy) %>%
    merge(dat_setting %>% select(n,Setting,power,interim_fraction), by='Setting') %>%
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
    select(n,power,MonitoringType,Iteration,ExpectedN_FO,ExpectedN_EO,ExpectedN_FE,Decision_FO,
           Decision_EO,Decision_FE,Decision_FE_F,Decision_FE_E,Decision_FO_25,Decision_EO_25,Decision_FE_25,Decision_FO_50,
           Decision_EO_50,Decision_FE_50,Decision_FO_75,Decision_EO_75,Decision_FE_75,
           RejectionRate_FO_FS,RejectionRate_EO_FS,RejectionRate_FE_FS) %>%
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
                      'DecisionFE_F','DecisionFE_E','DecisionFO_25','DecisionEO_25','DecisionFE_25',
                      'DecisionFO_50','DecisionEO_50','DecisionFE_50','DecisionFO_75','DecisionEO_75',
                      'DecisionFE_75','RejectionRate_FO_FS','RejectionRate_EO_FS',
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
           'RejectionRate_FE_FS' = RejectionRate_FE_FS_FE)#,

###############################################################
### APPLY M2 CALIBRATION CORRECTIONS TO GU RESULTS          ###
###############################################################

#apply M2 calibration correction stopping boundaries to GU results
dat_ABcalibrated_corrected_M2_upstraped_results <- merge(dat_results %>% select(Setting,Iteration,P_Upstrap), 
        dat_setting %>% select(Setting,n,interim_fraction,power) , by = 'Setting') %>%
    merge(dat_AB_setting_calibration_corrected %>% filter(interim_fraction != 1.00), by = c('n','interim_fraction')) %>%
    group_by(Setting,Iteration,MonitoringType) %>% 
    mutate(p_lessthan_fut = ifelse(is.na(FutilityP) == F & (as.numeric(P_Upstrap) < FutilityP), yes = 1, no = 0),
           p_lessthan_eff = ifelse(is.na(EfficacyP) == F & (as.numeric(P_Upstrap) < EfficacyP), yes = 1, no = 0)) %>%
    summarise(proportion_sim_fut = mean(p_lessthan_fut),
              proportion_sim_eff = mean(p_lessthan_eff)) %>%
    merge(dat_setting %>% select(Setting,n,interim_fraction), by = 'Setting') %>%
    merge(dat_AB_setting_calibration_corrected, by = c('n','interim_fraction','MonitoringType')) %>%
    mutate(stop_futility = ifelse(proportion_sim_fut < FutilityProp_multCorr2, yes = 1, no = 0),
           stop_efficacy = ifelse(proportion_sim_eff > EfficacyProp_multCorr2, yes = 1, no = 0)) %>%
    select(Setting,Iteration,MonitoringType,stop_futility,stop_efficacy) %>%
    merge(dat_setting %>% select(n,Setting,power,interim_fraction), by='Setting') %>%
    select(-Setting) %>%
    mutate(interim_fraction = factor(interim_fraction)) %>%
    mutate(stop_futility = ifelse(is.na(stop_futility) == F, yes = stop_futility, no = 0),
           stop_efficacy = ifelse(is.na(stop_efficacy) == F, yes = stop_efficacy, no = 0))

#use stopping boundary results to get expected sample size and full trial result statistics
dat_ABcalibrated_corrected_M2_upstraped_final_results <- dat_ABcalibrated_corrected_M2_upstraped_results %>%
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
    select(n,power,MonitoringType,Iteration,ExpectedN_FO,ExpectedN_EO,ExpectedN_FE,Decision_FO,
           Decision_EO,Decision_FE,Decision_FE_F,Decision_FE_E,Decision_FO_25,Decision_EO_25,Decision_FE_25,Decision_FO_50,
           Decision_EO_50,Decision_FE_50,Decision_FO_75,Decision_EO_75,Decision_FE_75,
           RejectionRate_FO_FS,RejectionRate_EO_FS,RejectionRate_FE_FS) %>%
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
                      'DecisionFE_F','DecisionFE_E','DecisionFO_25','DecisionEO_25','DecisionFE_25',
                      'DecisionFO_50','DecisionEO_50','DecisionFE_50','DecisionFO_75','DecisionEO_75',
                      'DecisionFE_75','RejectionRate_FO_FS','RejectionRate_EO_FS',
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
### APPLY M3 CALIBRATION CORRECTIONS TO GU RESULTS          ###
###############################################################

#apply M3 calibration correction stopping boundaries to GU results
dat_ABcalibrated_corrected_M3_upstraped_results <- merge(dat_results %>% select(Setting,Iteration,P_Upstrap), 
        dat_setting %>% select(Setting,n,interim_fraction,power) , by = 'Setting') %>%
    merge(dat_AB_setting_calibration_corrected %>% filter(interim_fraction != 1.00), by = c('n','interim_fraction')) %>%
    group_by(Setting,Iteration,MonitoringType) %>% 
    mutate(p_lessthan_fut = ifelse(is.na(FutilityP) == F & (as.numeric(P_Upstrap) < FutilityP), yes = 1, no = 0),
           p_lessthan_eff = ifelse(is.na(EfficacyP) == F & (as.numeric(P_Upstrap) < EfficacyP), yes = 1, no = 0)) %>%
    summarise(proportion_sim_fut = mean(p_lessthan_fut),
              proportion_sim_eff = mean(p_lessthan_eff)) %>%
    merge(dat_setting %>% select(Setting,n,interim_fraction), by = 'Setting') %>%
    merge(dat_AB_setting_calibration_corrected, by = c('n','interim_fraction','MonitoringType')) %>%
    mutate(stop_futility = ifelse(proportion_sim_fut < FutilityProp_multCorr3, yes = 1, no = 0),
           stop_efficacy = ifelse(proportion_sim_eff > EfficacyProp_multCorr3, yes = 1, no = 0))

#use stopping boundary results to get expected sample size and full trial result statistics
dat_ABcalibrated_corrected_M3_upstraped_final_results <- dat_ABcalibrated_corrected_M3_upstraped_results %>%
    select(Setting,Iteration,MonitoringType,stop_futility,stop_efficacy) %>%
    merge(dat_setting %>% select(n,Setting,power,interim_fraction), by='Setting') %>%
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
    select(n,power,MonitoringType,Iteration,ExpectedN_FO,ExpectedN_EO,ExpectedN_FE,Decision_FO,
           Decision_EO,Decision_FE,Decision_FE_F,Decision_FE_E,Decision_FO_25,Decision_EO_25,Decision_FE_25,Decision_FO_50,
           Decision_EO_50,Decision_FE_50,Decision_FO_75,Decision_EO_75,Decision_FE_75,
           RejectionRate_FO_FS,RejectionRate_EO_FS,RejectionRate_FE_FS) %>%
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
                      'DecisionFE_F','DecisionFE_E','DecisionFO_25','DecisionEO_25','DecisionFE_25',
                      'DecisionFO_50','DecisionEO_50','DecisionFE_50','DecisionFO_75','DecisionEO_75',
                      'DecisionFE_75','RejectionRate_FO_FS','RejectionRate_EO_FS',
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
### APPLY M4 CALIBRATION CORRECTIONS TO GU RESULTS          ###
###############################################################

#apply M14 calibration correction stopping boundaries to GU results
dat_ABcalibrated_corrected_M4_upstraped_results <- merge(dat_results %>% select(Setting,Iteration,P_Upstrap), 
        dat_setting %>% select(Setting,n,interim_fraction,power) , by = 'Setting') %>%
    merge(dat_AB_setting_calibration_corrected %>% filter(interim_fraction != 1.00), by = c('n','interim_fraction')) %>%
    group_by(Setting,Iteration,MonitoringType) %>% 
    mutate(p_lessthan_fut = ifelse(is.na(FutilityP) == F & (as.numeric(P_Upstrap) < FutilityP), yes = 1, no = 0),
           p_lessthan_eff = ifelse(is.na(EfficacyP) == F & (as.numeric(P_Upstrap) < EfficacyP), yes = 1, no = 0)) %>%
    summarise(proportion_sim_fut = mean(p_lessthan_fut),
              proportion_sim_eff = mean(p_lessthan_eff)) %>%
    merge(dat_setting %>% select(Setting,n,interim_fraction), by = 'Setting') %>%
    merge(dat_AB_setting_calibration_corrected, by = c('n','interim_fraction','MonitoringType')) %>%
    mutate(stop_futility = ifelse(proportion_sim_fut < FutilityProp_multCorr4, yes = 1, no = 0),
           stop_efficacy = ifelse(proportion_sim_eff > EfficacyProp_multCorr4, yes = 1, no = 0))

#use stopping boundary results to get expected sample size and full trial result statistics
dat_ABcalibrated_corrected_M4_upstraped_final_results <- dat_ABcalibrated_corrected_M4_upstraped_results%>%
    select(Setting,Iteration,MonitoringType,stop_futility,stop_efficacy) %>%
    merge(dat_setting %>% select(n,Setting,power,interim_fraction), by='Setting') %>%
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
    select(n,power,MonitoringType,Iteration,ExpectedN_FO,ExpectedN_EO,ExpectedN_FE,Decision_FO,
           Decision_EO,Decision_FE,Decision_FE_F,Decision_FE_E,Decision_FO_25,Decision_EO_25,Decision_FE_25,Decision_FO_50,
           Decision_EO_50,Decision_FE_50,Decision_FO_75,Decision_EO_75,Decision_FE_75,
           RejectionRate_FO_FS,RejectionRate_EO_FS,RejectionRate_FE_FS) %>%
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
                      'DecisionFE_F','DecisionFE_E','DecisionFO_25','DecisionEO_25','DecisionFE_25',
                      'DecisionFO_50','DecisionEO_50','DecisionFE_50','DecisionFO_75','DecisionEO_75',
                      'DecisionFE_75','RejectionRate_FO_FS','RejectionRate_EO_FS',
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
### APPLY M5 CALIBRATION CORRECTIONS TO GU RESULTS          ###
###############################################################

#apply M5 calibration correction stopping boundaries to GU results
dat_ABcalibrated_corrected_M5_upstraped_results <- merge(dat_results %>% select(Setting,Iteration,P_Upstrap), 
        dat_setting %>% select(Setting,n,interim_fraction,power) , by = 'Setting') %>%
    merge(dat_AB_setting_calibration_corrected %>% filter(interim_fraction != 1.00), by = c('n','interim_fraction')) %>%
    group_by(Setting,Iteration,MonitoringType) %>% 
    mutate(p_lessthan_fut = ifelse(is.na(FutilityP) == F & (as.numeric(P_Upstrap) < FutilityP), yes = 1, no = 0),
           p_lessthan_eff = ifelse(is.na(EfficacyP) == F & (as.numeric(P_Upstrap) < EfficacyP), yes = 1, no = 0)) %>%
    summarise(proportion_sim_fut = mean(p_lessthan_fut),
              proportion_sim_eff = mean(p_lessthan_eff)) %>%
    merge(dat_setting %>% select(Setting,n,interim_fraction), by = 'Setting') %>%
    merge(dat_AB_setting_calibration_corrected, by = c('n','interim_fraction','MonitoringType')) %>%
    mutate(stop_futility = ifelse(proportion_sim_fut < FutilityProp_multCorr5, yes = 1, no = 0),
           stop_efficacy = ifelse(proportion_sim_eff > EfficacyProp_multCorr5, yes = 1, no = 0))

#use stopping boundary results to get expected sample size and full trial result statistics
dat_ABcalibrated_corrected_M5_upstraped_final_results <- dat_ABcalibrated_corrected_M5_upstraped_results %>%
  select(Setting,Iteration,MonitoringType,stop_futility,stop_efficacy) %>%
  merge(dat_setting %>% select(n,Setting,power,interim_fraction), by='Setting') %>%
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
    select(n,power,MonitoringType,Iteration,ExpectedN_FO,ExpectedN_EO,ExpectedN_FE,Decision_FO,
           Decision_EO,Decision_FE,Decision_FE_F,Decision_FE_E,Decision_FO_25,Decision_EO_25,Decision_FE_25,Decision_FO_50,
           Decision_EO_50,Decision_FE_50,Decision_FO_75,Decision_EO_75,Decision_FE_75,
           RejectionRate_FO_FS,RejectionRate_EO_FS,RejectionRate_FE_FS) %>%
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
                      'DecisionFE_F','DecisionFE_E','DecisionFO_25','DecisionEO_25','DecisionFE_25',
                      'DecisionFO_50','DecisionEO_50','DecisionFE_50','DecisionFO_75','DecisionEO_75',
                      'DecisionFE_75','RejectionRate_FO_FS','RejectionRate_EO_FS',
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
### APPLY E1 CALIBRATION CORRECTIONS TO GU RESULTS          ###
###############################################################

#apply E1 calibration correction stopping boundaries to GU results
dat_ABcalibrated_corrected_E1_upstraped_results <- merge(dat_results %>% select(Setting,Iteration,P_Upstrap), 
        dat_setting %>% select(Setting,n,interim_fraction,power) , by = 'Setting') %>%
    merge(dat_AB_setting_calibration_corrected %>% filter(interim_fraction != 1.00), by = c('n','interim_fraction')) %>%
    group_by(Setting,Iteration,MonitoringType) %>% 
    mutate(p_lessthan_fut = ifelse(is.na(FutilityP) == F & (as.numeric(P_Upstrap) < FutilityP), yes = 1, no = 0),
           p_lessthan_eff = ifelse(is.na(EfficacyP) == F & (as.numeric(P_Upstrap) < EfficacyP), yes = 1, no = 0)) %>%
    summarise(proportion_sim_fut = mean(p_lessthan_fut),
              proportion_sim_eff = mean(p_lessthan_eff)) %>%
    merge(dat_setting %>% select(Setting,n,interim_fraction), by = 'Setting') %>%
    merge(dat_AB_setting_calibration_corrected, by = c('n','interim_fraction','MonitoringType')) %>%
    mutate(stop_futility = ifelse(proportion_sim_fut < FutilityProp_expCorr1, yes = 1, no = 0),
           stop_efficacy = ifelse(proportion_sim_eff > EfficacyProp_expCorr1, yes = 1, no = 0))
    
#use stopping boundary results to get expected sample size and full trial result statistics
dat_ABcalibrated_corrected_E1_upstraped_final_results <- dat_ABcalibrated_corrected_E1_upstraped_results %>%
  select(Setting,Iteration,MonitoringType,stop_futility,stop_efficacy) %>%
  merge(dat_setting %>% select(n,Setting,power,interim_fraction), by='Setting') %>%
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
    select(n,power,MonitoringType,Iteration,ExpectedN_FO,ExpectedN_EO,ExpectedN_FE,Decision_FO,
           Decision_EO,Decision_FE,Decision_FE_F,Decision_FE_E,Decision_FO_25,Decision_EO_25,Decision_FE_25,Decision_FO_50,
           Decision_EO_50,Decision_FE_50,Decision_FO_75,Decision_EO_75,Decision_FE_75,
           RejectionRate_FO_FS,RejectionRate_EO_FS,RejectionRate_FE_FS) %>%
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
                      'DecisionFE_F','DecisionFE_E','DecisionFO_25','DecisionEO_25','DecisionFE_25',
                      'DecisionFO_50','DecisionEO_50','DecisionFE_50','DecisionFO_75','DecisionEO_75',
                      'DecisionFE_75','RejectionRate_FO_FS','RejectionRate_EO_FS',
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
### APPLY E2 CALIBRATION CORRECTIONS TO GU RESULTS          ###
###############################################################

#apply E2 calibration correction stopping boundaries to GU results
dat_ABcalibrated_corrected_E2_upstraped_results <- merge(dat_results %>% select(Setting,Iteration,P_Upstrap), 
        dat_setting %>% select(Setting,n,interim_fraction,power) , by = 'Setting') %>%
    merge(dat_AB_setting_calibration_corrected %>% filter(interim_fraction != 1.00), by = c('n','interim_fraction')) %>%
    group_by(Setting,Iteration,MonitoringType) %>% 
    mutate(p_lessthan_fut = ifelse(is.na(FutilityP) == F & (as.numeric(P_Upstrap) < FutilityP), yes = 1, no = 0),
           p_lessthan_eff = ifelse(is.na(EfficacyP) == F & (as.numeric(P_Upstrap) < EfficacyP), yes = 1, no = 0)) %>%
    summarise(proportion_sim_fut = mean(p_lessthan_fut),
              proportion_sim_eff = mean(p_lessthan_eff)) %>%
    merge(dat_setting %>% select(Setting,n,interim_fraction), by = 'Setting') %>%
    merge(dat_AB_setting_calibration_corrected, by = c('n','interim_fraction','MonitoringType')) %>%
    mutate(stop_futility = ifelse(proportion_sim_fut < FutilityProp_expCorr2, yes = 1, no = 0),
           stop_efficacy = ifelse(proportion_sim_eff > EfficacyProp_expCorr2, yes = 1, no = 0))

#use stopping boundary results to get expected sample size and full trial result statistics
dat_ABcalibrated_corrected_E2_upstraped_final_results <- dat_ABcalibrated_corrected_E2_upstraped_results %>%
  select(Setting,Iteration,MonitoringType,stop_futility,stop_efficacy) %>%
  merge(dat_setting %>% select(n,Setting,power,interim_fraction), by='Setting') %>%
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
    select(n,power,MonitoringType,Iteration,ExpectedN_FO,ExpectedN_EO,ExpectedN_FE,Decision_FO,
           Decision_EO,Decision_FE,Decision_FE_F,Decision_FE_E,Decision_FO_25,Decision_EO_25,Decision_FE_25,Decision_FO_50,
           Decision_EO_50,Decision_FE_50,Decision_FO_75,Decision_EO_75,Decision_FE_75,
           RejectionRate_FO_FS,RejectionRate_EO_FS,RejectionRate_FE_FS) %>%
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
                      'DecisionFE_F','DecisionFE_E','DecisionFO_25','DecisionEO_25','DecisionFE_25',
                      'DecisionFO_50','DecisionEO_50','DecisionFE_50','DecisionFO_75','DecisionEO_75',
                      'DecisionFE_75','RejectionRate_FO_FS','RejectionRate_EO_FS',
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
### CONSOLIDATE AND STANDARDIZE DATA FORMATTING             ###
###############################################################

#combine data from all GU calibration corrections
dat_final_calibration_correction_GU <- rbind(
  dat_ABcalibrated_corrected_M1_upstraped_final_results %>% mutate(Method = 'AB Upstrap M1'),
  dat_ABcalibrated_corrected_M2_upstraped_final_results %>% mutate(Method = 'AB Upstrap M2'),
  dat_ABcalibrated_corrected_M3_upstraped_final_results %>% mutate(Method = 'AB Upstrap M3'),
  dat_ABcalibrated_corrected_M4_upstraped_final_results %>% mutate(Method = 'AB Upstrap M4'),
  dat_ABcalibrated_corrected_M5_upstraped_final_results %>% mutate(Method = 'AB Upstrap M5'),
  dat_ABcalibrated_corrected_E1_upstraped_final_results %>% mutate(Method = 'AB Upstrap E1'),
  dat_ABcalibrated_corrected_E2_upstraped_final_results %>% mutate(Method = 'AB Upstrap E2'),
  dat_final_output %>% filter(Method %in% c('OBrien Fleming GSD','Alpha Beta Spending'))
)

###############################################################
### SAVE OUTPUT DATASETS                                    ###
###############################################################

save(dat_AB_setting_calibration_corrected, file = './01 Data/Data Processed/dat_AB_setting_calibration_corrected.RData')
#load('./01 Data/Data Processed/dat_AB_setting_calibration_corrected.RData')
save(dat_ABcalibrated_corrected_M1_upstraped_results, file = './01 Data/Data Processed/dat_ABcalibrated_corrected_M1_upstraped_results.RData')
#load('./01 Data/Data Processed/dat_ABcalibrated_corrected_M1_upstraped_results.RData')
save(dat_ABcalibrated_corrected_M1_upstraped_final_results, file = './01 Data/Data Processed/dat_ABcalibrated_corrected_M1_upstraped_final_results.RData')
#load('./01 Data/Data Processed/dat_ABcalibrated_corrected_M1_upstraped_final_results.RData')
save(dat_ABcalibrated_corrected_M2_upstraped_results, file = './01 Data/Data Processed/dat_ABcalibrated_corrected_M2_upstraped_results.RData')
#load('./01 Data/Data Processed/dat_ABcalibrated_corrected_M2_upstraped_results.RData')
save(dat_ABcalibrated_corrected_M2_upstraped_final_results, file = './01 Data/Data Processed/dat_ABcalibrated_corrected_M2_upstraped_final_results.RData')
#load('./01 Data/Data Processed/dat_ABcalibrated_corrected_M2_upstraped_final_results.RData')
save(dat_ABcalibrated_corrected_M3_upstraped_results, file = './01 Data/Data Processed/dat_ABcalibrated_corrected_M3_upstraped_results.RData')
#load('./01 Data/Data Processed/dat_ABcalibrated_corrected_M3_upstraped_results.RData')
save(dat_ABcalibrated_corrected_M3_upstraped_final_results, file = './01 Data/Data Processed/dat_ABcalibrated_corrected_M3_upstraped_final_results.RData')
#load('./01 Data/Data Processed/dat_ABcalibrated_corrected_M3_upstraped_final_results.RData')
save(dat_ABcalibrated_corrected_M4_upstraped_results, file = './01 Data/Data Processed/dat_ABcalibrated_corrected_M4_upstraped_results.RData')
#load('./01 Data/Data Processed/dat_ABcalibrated_corrected_M4_upstraped_results.RData')
save(dat_ABcalibrated_corrected_M4_upstraped_final_results, file = './01 Data/Data Processed/dat_ABcalibrated_corrected_M4_upstraped_final_results.RData')
#load('./01 Data/Data Processed/dat_ABcalibrated_corrected_M4_upstraped_final_results.RData')
save(dat_ABcalibrated_corrected_M5_upstraped_results, file = './01 Data/Data Processed/dat_ABcalibrated_corrected_M5_upstraped_results.RData')
#load('./01 Data/Data Processed/dat_ABcalibrated_corrected_M5_upstraped_results.RData')
save(dat_ABcalibrated_corrected_M5_upstraped_final_results, file = './01 Data/Data Processed/dat_ABcalibrated_corrected_M5_upstraped_final_results.RData')
#load('./01 Data/Data Processed/dat_ABcalibrated_corrected_M5_upstraped_final_results.RData')
save(dat_ABcalibrated_corrected_E1_upstraped_results, file = './01 Data/Data Processed/dat_ABcalibrated_corrected_E1_upstraped_results.RData')
#load('./01 Data/Data Processed/dat_ABcalibrated_corrected_E1_upstraped_results.RData')
save(dat_ABcalibrated_corrected_E1_upstraped_final_results, file = './01 Data/Data Processed/dat_ABcalibrated_corrected_E1_upstraped_final_results.RData')
#load('./01 Data/Data Processed/dat_ABcalibrated_corrected_E1_upstraped_final_results.RData')
save(dat_ABcalibrated_corrected_E2_upstraped_results, file = './01 Data/Data Processed/dat_ABcalibrated_corrected_E2_upstraped_results.RData')
#load('./01 Data/Data Processed/dat_ABcalibrated_corrected_E2_upstraped_results.RData')
save(dat_ABcalibrated_corrected_E2_upstraped_final_results, file = './01 Data/Data Processed/dat_ABcalibrated_corrected_E2_upstraped_final_results.RData')
#load('./01 Data/Data Processed/dat_ABcalibrated_corrected_E2_upstraped_final_results.RData')
save(dat_final_calibration_correction_GU, file = './01 Data/Data Processed/dat_final_calibration_correction_GU.RData')
#load('./01 Data/Data Processed/dat_final_calibration_correction_GU.RData')

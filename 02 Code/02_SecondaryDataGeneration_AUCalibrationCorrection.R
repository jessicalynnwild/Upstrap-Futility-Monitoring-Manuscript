###############################################################
### file: 02_SecondaryDataGeneration_AUCalibrationCorrection###
### authors: Jess Wild                                      ###
### creation date:    07/25/2023                            ###
### latest edit date: 07/25/2023                            ###
### description: code to apply multiple testing correction  ###
### to AU upstrapping calibration                           ###
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
load('./01 Data/Data Processed/dat_fullsample_reject.RData')
load('./01 Data/Data Processed/dat_setting_calibration.RData')
load('./01 Data/Data Processed/dat_final_output.RData')

###############################################################
### APPLY M1 CALIBRATION CORRECTIONS TO AU RESULTS          ###
###############################################################

#apply M1 calibration correction stopping boundaries to AU results
dat_arbitrary_corrected_M1_upstraped_results <- dat_results %>%
  group_by(Setting,Iteration) %>% 
  mutate(p_lessthan_fut = ifelse(as.numeric(P_Upstrap) < 0.05, yes = 1, no = 0),
         p_lessthan_eff = ifelse(as.numeric(P_Upstrap) < 0.05, yes = 1, no = 0)) %>%
  summarise(proportion_sim_fut = mean(p_lessthan_fut),
            proportion_sim_eff = mean(p_lessthan_eff)) %>%
  merge(dat_setting_calibration %>% select(Setting,Prop_futility,Prop_efficacy), by = 'Setting') %>%
  mutate(stop_futility = ifelse(proportion_sim_fut < 0.05/(4+1), yes = 1, no = 0),
         stop_efficacy = ifelse(proportion_sim_eff > 0.80+(1-0.80)/(4+1), yes = 1, no = 0))

#use stopping boundary results to get expected sample size and full trial result statistics
dat_arbitrary_corrected_M1_upstraped_final_results <- dat_arbitrary_corrected_M1_upstraped_results %>%
  select(Setting,Iteration, stop_futility, stop_efficacy) %>%
  merge(dat_setting %>% select(Setting,n,interim_fraction,power), by='Setting') %>%
  select(-Setting) %>%
  mutate(interim_fraction = factor(interim_fraction)) %>%
  reshape(idvar = c("n","Iteration","power"), 
          timevar = "interim_fraction", 
          direction = "wide") %>%
  select(n,power,Iteration,stop_futility.0.25,stop_futility.0.5,stop_futility.0.75,stop_efficacy.0.25,stop_efficacy.0.5,stop_efficacy.0.75) %>%
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
  merge(dat_fullsample_reject, by=c('n','power','Iteration')) %>%
  mutate(RejectionRate_FO_FS = ifelse(Decision_FO==0&Reject_FS==1,yes=1,no=0),
         RejectionRate_EO_FS = ifelse(Decision_EO==1|Reject_FS==1,yes=1,no=0),
         RejectionRate_FE_FS = ifelse(Decision_FE_E==1|(Decision_FE_F==0&Reject_FS==1),yes=1,no=0)) %>%
  select(n,power,Iteration,ExpectedN_FO,ExpectedN_EO,ExpectedN_FE,Decision_FO,Decision_EO,Decision_FE,
         Decision_FE_F,Decision_FE_E,Decision_FO_25,Decision_EO_25,Decision_FE_25,Decision_FO_50,
         Decision_EO_50,Decision_FE_50,Decision_FO_75,Decision_EO_75,Decision_FE_75,RejectionRate_FO_FS,
         RejectionRate_EO_FS,RejectionRate_FE_FS) %>%
  group_by(n,power) %>%
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
            RejectionRate_FE_FS = mean(RejectionRate_FE_FS))

###############################################################
### APPLY M2 CALIBRATION CORRECTIONS TO AU RESULTS          ###
###############################################################

#apply M2 calibration correction stopping boundaries to AU results
dat_arbitrary_corrected_M2_upstraped_results <- dat_results %>%
  group_by(Setting,Iteration) %>% 
  mutate(p_lessthan_fut = ifelse(as.numeric(P_Upstrap) < 0.05, yes = 1, no = 0),
         p_lessthan_eff = ifelse(as.numeric(P_Upstrap) < 0.05, yes = 1, no = 0)) %>%
  summarise(proportion_sim_fut = mean(p_lessthan_fut),
            proportion_sim_eff = mean(p_lessthan_eff)) %>%
  merge(dat_setting_calibration %>% select(Setting,Prop_futility,Prop_efficacy), by = 'Setting') %>%
  mutate(stop_futility = ifelse(proportion_sim_fut < 0.05/(4), yes = 1, no = 0),
         stop_efficacy = ifelse(proportion_sim_eff > 0.80+(1-0.80)/(4), yes = 1, no = 0))

#use stopping boundary results to get expected sample size and full trial result statistics
dat_arbitrary_corrected_M2_upstraped_final_results <- dat_arbitrary_corrected_M2_upstraped_results %>%
  select(Setting,Iteration, stop_futility, stop_efficacy) %>%
  merge(dat_setting %>% select(Setting,n,interim_fraction,power), by='Setting') %>%
  select(-Setting) %>%
  mutate(interim_fraction = factor(interim_fraction)) %>%
  reshape(idvar = c("n","Iteration","power"), 
          timevar = "interim_fraction", 
          direction = "wide") %>%
  select(n,power,Iteration,stop_futility.0.25,stop_futility.0.5,stop_futility.0.75,stop_efficacy.0.25,stop_efficacy.0.5,stop_efficacy.0.75) %>%
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
  merge(dat_fullsample_reject, by=c('n','power','Iteration')) %>%
  mutate(RejectionRate_FO_FS = ifelse(Decision_FO==0&Reject_FS==1,yes=1,no=0),
         RejectionRate_EO_FS = ifelse(Decision_EO==1|Reject_FS==1,yes=1,no=0),
         RejectionRate_FE_FS = ifelse(Decision_FE_E==1|(Decision_FE_F==0&Reject_FS==1),yes=1,no=0)) %>%
  select(n,power,Iteration,ExpectedN_FO,ExpectedN_EO,ExpectedN_FE,Decision_FO,Decision_EO,Decision_FE,
         Decision_FE_F,Decision_FE_E,Decision_FO_25,Decision_EO_25,Decision_FE_25,Decision_FO_50,
         Decision_EO_50,Decision_FE_50,Decision_FO_75,Decision_EO_75,Decision_FE_75,RejectionRate_FO_FS,
         RejectionRate_EO_FS,RejectionRate_FE_FS) %>%
  group_by(n,power) %>%
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
            RejectionRate_FE_FS = mean(RejectionRate_FE_FS))

###############################################################
### APPLY M3 CALIBRATION CORRECTIONS TO AU RESULTS          ###
###############################################################

#apply M3 calibration correction stopping boundaries to AU results
dat_arbitrary_corrected_M3_upstraped_results <- dat_results %>%
  group_by(Setting,Iteration) %>% 
  mutate(p_lessthan_fut = ifelse(as.numeric(P_Upstrap) < 0.05, yes = 1, no = 0),
         p_lessthan_eff = ifelse(as.numeric(P_Upstrap) < 0.05, yes = 1, no = 0)) %>%
  summarise(proportion_sim_fut = mean(p_lessthan_fut),
            proportion_sim_eff = mean(p_lessthan_eff)) %>%
  merge(dat_setting_calibration %>% select(Setting,Prop_futility,Prop_efficacy), by = 'Setting') %>%
  mutate(stop_futility = ifelse(proportion_sim_fut < 0.05/(4-1), yes = 1, no = 0),
         stop_efficacy = ifelse(proportion_sim_eff > 0.80+(1-0.80)/(4-1), yes = 1, no = 0))

#use stopping boundary results to get expected sample size and full trial result statistics
dat_arbitrary_corrected_M3_upstraped_final_results <- dat_arbitrary_corrected_M3_upstraped_results %>%
  select(Setting,Iteration, stop_futility, stop_efficacy) %>%
  merge(dat_setting %>% select(Setting,n,interim_fraction,power), by='Setting') %>%
  select(-Setting) %>%
  mutate(interim_fraction = factor(interim_fraction)) %>%
  reshape(idvar = c("n","Iteration","power"), 
          timevar = "interim_fraction", 
          direction = "wide") %>%
  select(n,power,Iteration,stop_futility.0.25,stop_futility.0.5,stop_futility.0.75,stop_efficacy.0.25,stop_efficacy.0.5,stop_efficacy.0.75) %>%
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
  merge(dat_fullsample_reject, by=c('n','power','Iteration')) %>%
  mutate(RejectionRate_FO_FS = ifelse(Decision_FO==0&Reject_FS==1,yes=1,no=0),
         RejectionRate_EO_FS = ifelse(Decision_EO==1|Reject_FS==1,yes=1,no=0),
         RejectionRate_FE_FS = ifelse(Decision_FE_E==1|(Decision_FE_F==0&Reject_FS==1),yes=1,no=0)) %>%
  select(n,power,Iteration,ExpectedN_FO,ExpectedN_EO,ExpectedN_FE,Decision_FO,Decision_EO,Decision_FE,
         Decision_FE_F,Decision_FE_E,Decision_FO_25,Decision_EO_25,Decision_FE_25,Decision_FO_50,
         Decision_EO_50,Decision_FE_50,Decision_FO_75,Decision_EO_75,Decision_FE_75,RejectionRate_FO_FS,
         RejectionRate_EO_FS,RejectionRate_FE_FS) %>%
  group_by(n,power) %>%
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
            RejectionRate_FE_FS = mean(RejectionRate_FE_FS))

###############################################################
### APPLY M4 CALIBRATION CORRECTIONS TO AU RESULTS          ###
###############################################################

#apply M4 calibration correction stopping boundaries to AU results
dat_arbitrary_corrected_M4_upstraped_results <- dat_results %>%
  group_by(Setting,Iteration) %>% 
  mutate(p_lessthan_fut = ifelse(as.numeric(P_Upstrap) < 0.05, yes = 1, no = 0),
         p_lessthan_eff = ifelse(as.numeric(P_Upstrap) < 0.05, yes = 1, no = 0)) %>%
  summarise(proportion_sim_fut = mean(p_lessthan_fut),
            proportion_sim_eff = mean(p_lessthan_eff)) %>%
  merge(dat_setting_calibration %>% select(Setting,Prop_futility,Prop_efficacy), by = 'Setting') %>%
  mutate(stop_futility = ifelse(proportion_sim_fut < 0.05/(4-1.5), yes = 1, no = 0),
         stop_efficacy = ifelse(proportion_sim_eff > 0.80+(1-0.80)/(4-1.5), yes = 1, no = 0))

#use stopping boundary results to get expected sample size and full trial result statistics
dat_arbitrary_corrected_M4_upstraped_final_results <- dat_arbitrary_corrected_M4_upstraped_results %>%
  select(Setting,Iteration, stop_futility, stop_efficacy) %>%
  merge(dat_setting %>% select(Setting,n,interim_fraction,power), by='Setting') %>%
  select(-Setting) %>%
  mutate(interim_fraction = factor(interim_fraction)) %>%
  reshape(idvar = c("n","Iteration","power"), 
          timevar = "interim_fraction", 
          direction = "wide") %>%
  select(n,power,Iteration,stop_futility.0.25,stop_futility.0.5,stop_futility.0.75,stop_efficacy.0.25,stop_efficacy.0.5,stop_efficacy.0.75) %>%
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
  merge(dat_fullsample_reject, by=c('n','power','Iteration')) %>%
  mutate(RejectionRate_FO_FS = ifelse(Decision_FO==0&Reject_FS==1,yes=1,no=0),
         RejectionRate_EO_FS = ifelse(Decision_EO==1|Reject_FS==1,yes=1,no=0),
         RejectionRate_FE_FS = ifelse(Decision_FE_E==1|(Decision_FE_F==0&Reject_FS==1),yes=1,no=0)) %>%
  select(n,power,Iteration,ExpectedN_FO,ExpectedN_EO,ExpectedN_FE,Decision_FO,Decision_EO,Decision_FE,
         Decision_FE_F,Decision_FE_E,Decision_FO_25,Decision_EO_25,Decision_FE_25,Decision_FO_50,
         Decision_EO_50,Decision_FE_50,Decision_FO_75,Decision_EO_75,Decision_FE_75,RejectionRate_FO_FS,
         RejectionRate_EO_FS,RejectionRate_FE_FS) %>%
  group_by(n,power) %>%
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
            RejectionRate_FE_FS = mean(RejectionRate_FE_FS))

###############################################################
### APPLY M5 CALIBRATION CORRECTIONS TO AU RESULTS          ###
###############################################################

#apply M5 calibration correction stopping boundaries to AU results
dat_arbitrary_corrected_M5_upstraped_results <- dat_results %>%
  group_by(Setting,Iteration) %>% 
  mutate(p_lessthan_fut = ifelse(as.numeric(P_Upstrap) < 0.05, yes = 1, no = 0),
         p_lessthan_eff = ifelse(as.numeric(P_Upstrap) < 0.05, yes = 1, no = 0)) %>%
  summarise(proportion_sim_fut = mean(p_lessthan_fut),
            proportion_sim_eff = mean(p_lessthan_eff)) %>%
  merge(dat_setting_calibration %>% select(Setting,Prop_futility,Prop_efficacy), by = 'Setting') %>%
  mutate(stop_futility = ifelse(proportion_sim_fut < 0.05/(4-2), yes = 1, no = 0),
         stop_efficacy = ifelse(proportion_sim_eff > 0.80+(1-0.80)/(4-2), yes = 1, no = 0))

#use stopping boundary results to get expected sample size and full trial result statistics
dat_arbitrary_corrected_M5_upstraped_final_results <- dat_arbitrary_corrected_M5_upstraped_results%>%
  select(Setting,Iteration, stop_futility, stop_efficacy) %>%
  merge(dat_setting %>% select(Setting,n,interim_fraction,power), by='Setting') %>%
  select(-Setting) %>%
  mutate(interim_fraction = factor(interim_fraction)) %>%
  reshape(idvar = c("n","Iteration","power"), 
          timevar = "interim_fraction", 
          direction = "wide") %>%
  select(n,power,Iteration,stop_futility.0.25,stop_futility.0.5,stop_futility.0.75,stop_efficacy.0.25,stop_efficacy.0.5,stop_efficacy.0.75) %>%
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
  merge(dat_fullsample_reject, by=c('n','power','Iteration')) %>%
  mutate(RejectionRate_FO_FS = ifelse(Decision_FO==0&Reject_FS==1,yes=1,no=0),
         RejectionRate_EO_FS = ifelse(Decision_EO==1|Reject_FS==1,yes=1,no=0),
         RejectionRate_FE_FS = ifelse(Decision_FE_E==1|(Decision_FE_F==0&Reject_FS==1),yes=1,no=0)) %>%
  select(n,power,Iteration,ExpectedN_FO,ExpectedN_EO,ExpectedN_FE,Decision_FO,Decision_EO,Decision_FE,
         Decision_FE_F,Decision_FE_E,Decision_FO_25,Decision_EO_25,Decision_FE_25,Decision_FO_50,
         Decision_EO_50,Decision_FE_50,Decision_FO_75,Decision_EO_75,Decision_FE_75,RejectionRate_FO_FS,
         RejectionRate_EO_FS,RejectionRate_FE_FS) %>%
  group_by(n,power) %>%
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
            RejectionRate_FE_FS = mean(RejectionRate_FE_FS))

###############################################################
### APPLY E1 CALIBRATION CORRECTIONS TO AU RESULTS          ###
###############################################################

#apply E1 calibration correction stopping boundaries to AU results
dat_arbitrary_corrected_E1_upstraped_results <- dat_results %>%
  group_by(Setting,Iteration) %>% 
  mutate(p_lessthan_fut = ifelse(as.numeric(P_Upstrap) < 0.05, yes = 1, no = 0),
         p_lessthan_eff = ifelse(as.numeric(P_Upstrap) < 0.05, yes = 1, no = 0)) %>%
  summarise(proportion_sim_fut = mean(p_lessthan_fut),
            proportion_sim_eff = mean(p_lessthan_eff)) %>%
  merge(dat_setting_calibration %>% select(Setting,Prop_futility,Prop_efficacy), by = 'Setting') %>%
  mutate(stop_futility = ifelse(proportion_sim_fut < 0.05^0.95, yes = 1, no = 0),
         stop_efficacy = ifelse(proportion_sim_eff > 0.80^(1-0.95), yes = 1, no = 0))

#use stopping boundary results to get expected sample size and full trial result statistics
dat_arbitrary_corrected_E1_upstraped_final_results <- dat_arbitrary_corrected_E1_upstraped_results %>%
  select(Setting,Iteration, stop_futility, stop_efficacy) %>%
  merge(dat_setting %>% select(Setting,n,interim_fraction,power), by='Setting') %>%
  select(-Setting) %>%
  mutate(interim_fraction = factor(interim_fraction)) %>%
  reshape(idvar = c("n","Iteration","power"), 
          timevar = "interim_fraction", 
          direction = "wide") %>%
  select(n,power,Iteration,stop_futility.0.25,stop_futility.0.5,stop_futility.0.75,stop_efficacy.0.25,stop_efficacy.0.5,stop_efficacy.0.75) %>%
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
  merge(dat_fullsample_reject, by=c('n','power','Iteration')) %>%
  mutate(RejectionRate_FO_FS = ifelse(Decision_FO==0&Reject_FS==1,yes=1,no=0),
         RejectionRate_EO_FS = ifelse(Decision_EO==1|Reject_FS==1,yes=1,no=0),
         RejectionRate_FE_FS = ifelse(Decision_FE_E==1|(Decision_FE_F==0&Reject_FS==1),yes=1,no=0)) %>%
  select(n,power,Iteration,ExpectedN_FO,ExpectedN_EO,ExpectedN_FE,Decision_FO,Decision_EO,Decision_FE,
         Decision_FE_F,Decision_FE_E,Decision_FO_25,Decision_EO_25,Decision_FE_25,Decision_FO_50,
         Decision_EO_50,Decision_FE_50,Decision_FO_75,Decision_EO_75,Decision_FE_75,RejectionRate_FO_FS,
         RejectionRate_EO_FS,RejectionRate_FE_FS) %>%
  group_by(n,power) %>%
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
            RejectionRate_FE_FS = mean(RejectionRate_FE_FS))

###############################################################
### APPLY E2 CALIBRATION CORRECTIONS TO AU RESULTS          ###
###############################################################

#apply E2 calibration correction stopping boundaries to AU results
dat_arbitrary_corrected_E2_upstraped_results <- dat_results %>%
  group_by(Setting,Iteration) %>% 
  mutate(p_lessthan_fut = ifelse(as.numeric(P_Upstrap) < 0.05, yes = 1, no = 0),
         p_lessthan_eff = ifelse(as.numeric(P_Upstrap) < 0.05, yes = 1, no = 0)) %>%
  summarise(proportion_sim_fut = mean(p_lessthan_fut),
            proportion_sim_eff = mean(p_lessthan_eff)) %>%
  merge(dat_setting_calibration %>% select(Setting,Prop_futility,Prop_efficacy), by = 'Setting') %>%
  mutate(stop_futility = ifelse(proportion_sim_fut < 0.05^0.70, yes = 1, no = 0),
         stop_efficacy = ifelse(proportion_sim_eff > 0.80^(1-0.70), yes = 1, no = 0)) 

#use stopping boundary results to get expected sample size and full trial result statistics
dat_arbitrary_corrected_E2_upstraped_final_results <- dat_arbitrary_corrected_E2_upstraped_results %>%
  select(Setting,Iteration, stop_futility, stop_efficacy) %>%
  merge(dat_setting %>% select(Setting,n,interim_fraction,power), by='Setting') %>%
  select(-Setting) %>%
  mutate(interim_fraction = factor(interim_fraction)) %>%
  reshape(idvar = c("n","Iteration","power"), 
          timevar = "interim_fraction", 
          direction = "wide") %>%
  select(n,power,Iteration,stop_futility.0.25,stop_futility.0.5,stop_futility.0.75,stop_efficacy.0.25,stop_efficacy.0.5,stop_efficacy.0.75) %>%
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
  merge(dat_fullsample_reject, by=c('n','power','Iteration')) %>%
  mutate(RejectionRate_FO_FS = ifelse(Decision_FO==0&Reject_FS==1,yes=1,no=0),
         RejectionRate_EO_FS = ifelse(Decision_EO==1|Reject_FS==1,yes=1,no=0),
         RejectionRate_FE_FS = ifelse(Decision_FE_E==1|(Decision_FE_F==0&Reject_FS==1),yes=1,no=0)) %>%
  select(n,power,Iteration,ExpectedN_FO,ExpectedN_EO,ExpectedN_FE,Decision_FO,Decision_EO,Decision_FE,
         Decision_FE_F,Decision_FE_E,Decision_FO_25,Decision_EO_25,Decision_FE_25,Decision_FO_50,
         Decision_EO_50,Decision_FE_50,Decision_FO_75,Decision_EO_75,Decision_FE_75,RejectionRate_FO_FS,
         RejectionRate_EO_FS,RejectionRate_FE_FS) %>%
  group_by(n,power) %>%
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
            RejectionRate_FE_FS = mean(RejectionRate_FE_FS))

###############################################################
### CONSOLIDATE AND STANDARDIZE DATA FORMATTING             ###
###############################################################

#combine data from all AU calibration corrections
dat_final_calibration_correction_AU <- rbind(
  dat_arbitrary_corrected_M1_upstraped_final_results %>% mutate(Method = 'Arbitrary Upstrap M1'),
  dat_arbitrary_corrected_M2_upstraped_final_results %>% mutate(Method = 'Arbitrary Upstrap M2'),
  dat_arbitrary_corrected_M3_upstraped_final_results %>% mutate(Method = 'Arbitrary Upstrap M3'),
  dat_arbitrary_corrected_M4_upstraped_final_results %>% mutate(Method = 'Arbitrary Upstrap M4'),
  dat_arbitrary_corrected_M5_upstraped_final_results %>% mutate(Method = 'Arbitrary Upstrap M5'),
  dat_arbitrary_corrected_E1_upstraped_final_results %>% mutate(Method = 'Arbitrary Upstrap E1'),
  dat_arbitrary_corrected_E2_upstraped_final_results %>% mutate(Method = 'Arbitrary Upstrap E2'),
  dat_final_output %>% filter(Method %in% c('OBrien Fleming GSD','Arbitrary Upstrap'))
)

###############################################################
### SAVE OUTPUT DATASETS                                    ###
###############################################################

save(dat_arbitrary_corrected_M1_upstraped_results, file = './01 Data/Data Processed/dat_arbitrary_corrected_M1_upstraped_results.RData')
#load('./01 Data/Data Processed/dat_arbitrary_corrected_M1_upstraped_results.RData')
save(dat_arbitrary_corrected_M1_upstraped_final_results, file = './01 Data/Data Processed/dat_arbitrary_corrected_M1_upstraped_final_results.RData')
#load('./01 Data/Data Processed/dat_arbitrary_corrected_M1_upstraped_final_results.RData')
save(dat_arbitrary_corrected_M2_upstraped_results, file = './01 Data/Data Processed/dat_arbitrary_corrected_M2_upstraped_results.RData')
#load('./01 Data/Data Processed/dat_arbitrary_corrected_M2_upstraped_results.RData')
save(dat_arbitrary_corrected_M2_upstraped_final_results, file = './01 Data/Data Processed/dat_arbitrary_corrected_M2_upstraped_final_results.RData')
#load('./01 Data/Data Processed/dat_arbitrary_corrected_M2_upstraped_final_results.RData')
save(dat_arbitrary_corrected_M3_upstraped_results, file = './01 Data/Data Processed/dat_arbitrary_corrected_M3_upstraped_results.RData')
#load('./01 Data/Data Processed/dat_arbitrary_corrected_M3_upstraped_results.RData')
save(dat_arbitrary_corrected_M3_upstraped_final_results, file = './01 Data/Data Processed/dat_arbitrary_corrected_M3_upstraped_final_results.RData')
#load('./01 Data/Data Processed/dat_arbitrary_corrected_M3_upstraped_final_results.RData')
save(dat_arbitrary_corrected_M4_upstraped_results, file = './01 Data/Data Processed/dat_arbitrary_corrected_M4_upstraped_results.RData')
#load('./01 Data/Data Processed/dat_arbitrary_corrected_M4_upstraped_results.RData')
save(dat_arbitrary_corrected_M4_upstraped_final_results, file = './01 Data/Data Processed/dat_arbitrary_corrected_M4_upstraped_final_results.RData')
#load('./01 Data/Data Processed/dat_arbitrary_corrected_M4_upstraped_final_results.RData')
save(dat_arbitrary_corrected_M5_upstraped_results, file = './01 Data/Data Processed/dat_arbitrary_corrected_M5_upstraped_results.RData')
#load('./01 Data/Data Processed/dat_arbitrary_corrected_M5_upstraped_results.RData')
save(dat_arbitrary_corrected_M5_upstraped_final_results, file = './01 Data/Data Processed/dat_arbitrary_corrected_M5_upstraped_final_results.RData')
#load('./01 Data/Data Processed/dat_arbitrary_corrected_M5_upstraped_final_results.RData')
save(dat_arbitrary_corrected_E1_upstraped_results, file = './01 Data/Data Processed/dat_arbitrary_corrected_E1_upstraped_results.RData')
#load('./01 Data/Data Processed/dat_arbitrary_corrected_E1_upstraped_results.RData')
save(dat_arbitrary_corrected_E1_upstraped_final_results, file = './01 Data/Data Processed/dat_arbitrary_corrected_E1_upstraped_final_results.RData')
#load('./01 Data/Data Processed/dat_arbitrary_corrected_E1_upstraped_final_results.RData')
save(dat_arbitrary_corrected_E2_upstraped_results, file = './01 Data/Data Processed/dat_arbitrary_corrected_E2_upstraped_results.RData')
#load('./01 Data/Data Processed/dat_arbitrary_corrected_E2_upstraped_results.RData')
save(dat_arbitrary_corrected_E2_upstraped_final_results, file = './01 Data/Data Processed/dat_arbitrary_corrected_E2_upstraped_final_results.RData')
#load('./01 Data/Data Processed/dat_arbitrary_corrected_E2_upstraped_final_results.RData')
save(dat_final_calibration_correction_AU, file = './01 Data/Data Processed/dat_final_calibration_correction_AU.RData')
#load('./01 Data/Data Processed/dat_final_calibration_correction_AU.RData')

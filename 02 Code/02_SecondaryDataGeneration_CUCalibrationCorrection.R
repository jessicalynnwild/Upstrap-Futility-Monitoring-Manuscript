###############################################################
### file: 02_SecondaryDataGeneration_CUCalibrationCorrection###
### authors: Jess Wild                                      ###
### creation date:    07/25/2023                            ###
### latest edit date: 07/25/2023                            ###
### description: code to apply multiple testing correction  ###
### to CU upstrapping calibration                           ###
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
### DEFINE CORRECTED PROPORTION THRESHOLDS                  ###
###############################################################

#define calibration thresholds for the calibration corrections
dat_setting_calibration_corrected <- dat_setting_calibration %>%
  mutate(Prop_futility_multCorr1 = Prop_futility/(4+1),
         Prop_futility_multCorr2 = Prop_futility/(4),
         Prop_futility_multCorr3 = Prop_futility/(4-1),
         Prop_futility_multCorr4 = Prop_futility/(4-1.5),
         Prop_futility_multCorr5 = Prop_futility/(4-2),
         Prop_futility_expCorr1 = Prop_futility^0.95,
         Prop_futility_expCorr2 = Prop_futility^0.70,
         Prop_efficacy_multCorr1 = Prop_efficacy+(1-Prop_efficacy)/(4+1),
         Prop_efficacy_multCorr2 = Prop_efficacy+(1-Prop_efficacy)/(4),
         Prop_efficacy_multCorr3 = Prop_efficacy+(1-Prop_efficacy)/(4-1),
         Prop_efficacy_multCorr4 = Prop_efficacy+(1-Prop_efficacy)/(4-1.5),
         Prop_efficacy_multCorr5 = Prop_efficacy+(1-Prop_efficacy)/(4-2),
         Prop_efficacy_expCorr1 = Prop_efficacy^(1-0.95),
         Prop_efficacy_expCorr2 = Prop_efficacy^(1-0.70))

###############################################################
### APPLY M1 CALIBRATION CORRECTIONS TO CU RESULTS          ###
###############################################################

#apply M1 calibration correction stopping boundaries to CU results
dat_calibrated_corrected_M1_upstraped_results <- merge(dat_results, dat_setting_calibration_corrected %>% 
                                                         select(Setting,P_futility,Prop_efficacy_multCorr1,P_efficacy,Prop_futility_multCorr1) %>%
                                                         rename(Prop_efficacy = Prop_efficacy_multCorr1,
                                                                Prop_futility = Prop_futility_multCorr1), 
                                                       by = 'Setting') %>%
  group_by(Setting,Iteration) %>% 
  mutate(p_lessthan_fut = ifelse(as.numeric(P_Upstrap) < P_futility, yes = 1, no = 0),
         p_lessthan_eff = ifelse(as.numeric(P_Upstrap) < P_efficacy, yes = 1, no = 0)) %>%
  summarise(proportion_sim_fut = mean(p_lessthan_fut),
            proportion_sim_eff = mean(p_lessthan_eff)) %>%
  merge(dat_setting_calibration_corrected %>% 
          select(Setting,P_futility,Prop_efficacy_multCorr1,P_efficacy,Prop_futility_multCorr1) %>%
          rename(Prop_efficacy = Prop_efficacy_multCorr1,
                 Prop_futility = Prop_futility_multCorr1), 
        by = 'Setting') %>%
  mutate(stop_futility = ifelse(proportion_sim_fut < Prop_futility, yes = 1, no = 0),
         stop_efficacy = ifelse(proportion_sim_eff > Prop_efficacy, yes = 1, no = 0))

#use stopping boundary results to get expected sample size and full trial result statistics
dat_calibrated_corrected_M1_upstraped_final_results <- dat_calibrated_corrected_M1_upstraped_results %>%
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
### APPLY M2 CALIBRATION CORRECTIONS TO CU RESULTS          ###
###############################################################

#apply M2 calibration correction stopping boundaries to CU results
dat_calibrated_corrected_M2_upstraped_results <- merge(dat_results, dat_setting_calibration_corrected %>% 
                                                         select(Setting,P_futility,Prop_efficacy_multCorr2,P_efficacy,Prop_futility_multCorr2) %>%
                                                         rename(Prop_efficacy = Prop_efficacy_multCorr2,
                                                                Prop_futility = Prop_futility_multCorr2), 
                                                       by = 'Setting') %>%
  group_by(Setting,Iteration) %>% 
  mutate(p_lessthan_fut = ifelse(as.numeric(P_Upstrap) < P_futility, yes = 1, no = 0),
         p_lessthan_eff = ifelse(as.numeric(P_Upstrap) < P_efficacy, yes = 1, no = 0)) %>%
  summarise(proportion_sim_fut = mean(p_lessthan_fut),
            proportion_sim_eff = mean(p_lessthan_eff)) %>%
  merge(dat_setting_calibration_corrected %>% 
          select(Setting,P_futility,Prop_efficacy_multCorr2,P_efficacy,Prop_futility_multCorr2) %>%
          rename(Prop_efficacy = Prop_efficacy_multCorr2,
                 Prop_futility = Prop_futility_multCorr2), 
        by = 'Setting') %>%
  mutate(stop_futility = ifelse(proportion_sim_fut < Prop_futility, yes = 1, no = 0),
         stop_efficacy = ifelse(proportion_sim_eff > Prop_efficacy, yes = 1, no = 0)) 

#use stopping boundary results to get expected sample size and full trial result statistics
dat_calibrated_corrected_M2_upstraped_final_results <- dat_calibrated_corrected_M2_upstraped_results %>%
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
### APPLY M3 CALIBRATION CORRECTIONS TO CU RESULTS          ###
###############################################################

#apply M3 calibration correction stopping boundaries to CU results
dat_calibrated_corrected_M3_upstraped_results <- merge(dat_results, dat_setting_calibration_corrected %>% 
                                                         select(Setting,P_futility,Prop_efficacy_multCorr3,P_efficacy,Prop_futility_multCorr3) %>%
                                                         rename(Prop_efficacy = Prop_efficacy_multCorr3,
                                                                Prop_futility = Prop_futility_multCorr3), 
                                                       by = 'Setting') %>%
  group_by(Setting,Iteration) %>% 
  mutate(p_lessthan_fut = ifelse(as.numeric(P_Upstrap) < P_futility, yes = 1, no = 0),
         p_lessthan_eff = ifelse(as.numeric(P_Upstrap) < P_efficacy, yes = 1, no = 0)) %>%
  summarise(proportion_sim_fut = mean(p_lessthan_fut),
            proportion_sim_eff = mean(p_lessthan_eff)) %>%
  merge(dat_setting_calibration_corrected %>% 
          select(Setting,P_futility,Prop_efficacy_multCorr3,P_efficacy,Prop_futility_multCorr3) %>%
          rename(Prop_efficacy = Prop_efficacy_multCorr3,
                 Prop_futility = Prop_futility_multCorr3), 
        by = 'Setting') %>%
  mutate(stop_futility = ifelse(proportion_sim_fut < Prop_futility, yes = 1, no = 0),
         stop_efficacy = ifelse(proportion_sim_eff > Prop_efficacy, yes = 1, no = 0))

#use stopping boundary results to get expected sample size and full trial result statistics
dat_calibrated_corrected_M3_upstraped_final_results <- dat_calibrated_corrected_M3_upstraped_results %>%
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
### APPLY M4 CALIBRATION CORRECTIONS TO CU RESULTS          ###
###############################################################

#apply M4 calibration correction stopping boundaries to CU results
dat_calibrated_corrected_M4_upstraped_results <- merge(dat_results, dat_setting_calibration_corrected %>% 
                                                         select(Setting,P_futility,Prop_efficacy_multCorr4,P_efficacy,Prop_futility_multCorr4) %>%
                                                         rename(Prop_efficacy = Prop_efficacy_multCorr4,
                                                                Prop_futility = Prop_futility_multCorr4), 
                                                       by = 'Setting') %>%
  group_by(Setting,Iteration) %>% 
  mutate(p_lessthan_fut = ifelse(as.numeric(P_Upstrap) < P_futility, yes = 1, no = 0),
         p_lessthan_eff = ifelse(as.numeric(P_Upstrap) < P_efficacy, yes = 1, no = 0)) %>%
  summarise(proportion_sim_fut = mean(p_lessthan_fut),
            proportion_sim_eff = mean(p_lessthan_eff)) %>%
  merge(dat_setting_calibration_corrected %>% 
          select(Setting,P_futility,Prop_efficacy_multCorr4,P_efficacy,Prop_futility_multCorr4) %>%
          rename(Prop_efficacy = Prop_efficacy_multCorr4,
                 Prop_futility = Prop_futility_multCorr4), 
        by = 'Setting') %>%
  mutate(stop_futility = ifelse(proportion_sim_fut < Prop_futility, yes = 1, no = 0),
         stop_efficacy = ifelse(proportion_sim_eff > Prop_efficacy, yes = 1, no = 0)) 

#use stopping boundary results to get expected sample size and full trial result statistics
dat_calibrated_corrected_M4_upstraped_final_results <- dat_calibrated_corrected_M4_upstraped_results%>%
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
### APPLY M5 CALIBRATION CORRECTIONS TO CU RESULTS          ###
###############################################################

#apply M5 calibration correction stopping boundaries to CU results
dat_calibrated_corrected_M5_upstraped_results <- merge(dat_results, dat_setting_calibration_corrected %>% 
                                                         select(Setting,P_futility,Prop_efficacy_multCorr5,P_efficacy,Prop_futility_multCorr5) %>%
                                                         rename(Prop_efficacy = Prop_efficacy_multCorr5,
                                                                Prop_futility = Prop_futility_multCorr5), 
                                                       by = 'Setting') %>%
  group_by(Setting,Iteration) %>% 
  mutate(p_lessthan_fut = ifelse(as.numeric(P_Upstrap) < P_futility, yes = 1, no = 0),
         p_lessthan_eff = ifelse(as.numeric(P_Upstrap) < P_efficacy, yes = 1, no = 0)) %>%
  summarise(proportion_sim_fut = mean(p_lessthan_fut),
            proportion_sim_eff = mean(p_lessthan_eff)) %>%
  merge(dat_setting_calibration_corrected %>% 
          select(Setting,P_futility,Prop_efficacy_multCorr5,P_efficacy,Prop_futility_multCorr5) %>%
          rename(Prop_efficacy = Prop_efficacy_multCorr5,
                 Prop_futility = Prop_futility_multCorr5), 
        by = 'Setting') %>%
  mutate(stop_futility = ifelse(proportion_sim_fut < Prop_futility, yes = 1, no = 0),
         stop_efficacy = ifelse(proportion_sim_eff > Prop_efficacy, yes = 1, no = 0))

#use stopping boundary results to get expected sample size and full trial result statistics
dat_calibrated_corrected_M5_upstraped_final_results <- dat_calibrated_corrected_M5_upstraped_results %>%
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
### APPLY E1 CALIBRATION CORRECTIONS TO CU RESULTS          ###
###############################################################

#apply E1 calibration correction stopping boundaries to CU results
dat_calibrated_corrected_E1_upstraped_results <- merge(dat_results, dat_setting_calibration_corrected %>% 
                                                         select(Setting,P_futility,Prop_efficacy_expCorr1,P_efficacy,Prop_futility_expCorr1) %>%
                                                         rename(Prop_efficacy = Prop_efficacy_expCorr1,
                                                                Prop_futility = Prop_futility_expCorr1), 
                                                       by = 'Setting') %>%
  group_by(Setting,Iteration) %>% 
  mutate(p_lessthan_fut = ifelse(as.numeric(P_Upstrap) < P_futility, yes = 1, no = 0),
         p_lessthan_eff = ifelse(as.numeric(P_Upstrap) < P_efficacy, yes = 1, no = 0)) %>%
  summarise(proportion_sim_fut = mean(p_lessthan_fut),
            proportion_sim_eff = mean(p_lessthan_eff)) %>%
  merge(dat_setting_calibration_corrected %>% 
          select(Setting,P_futility,Prop_efficacy_expCorr1,P_efficacy,Prop_futility_expCorr1) %>%
          rename(Prop_efficacy = Prop_efficacy_expCorr1,
                 Prop_futility = Prop_futility_expCorr1), 
        by = 'Setting') %>%
  mutate(stop_futility = ifelse(proportion_sim_fut < Prop_futility, yes = 1, no = 0),
         stop_efficacy = ifelse(proportion_sim_eff > Prop_efficacy, yes = 1, no = 0))

#use stopping boundary results to get expected sample size and full trial result statistics
dat_calibrated_corrected_E1_upstraped_final_results <- dat_calibrated_corrected_E1_upstraped_results %>%
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
### APPLY E2 CALIBRATION CORRECTIONS TO CU RESULTS          ###
###############################################################

#apply E2 calibration correction stopping boundaries to CU results
dat_calibrated_corrected_E2_upstraped_results <- merge(dat_results, dat_setting_calibration_corrected %>% 
                                                         select(Setting,P_futility,Prop_efficacy_expCorr2,P_efficacy,Prop_futility_expCorr2) %>%
                                                         rename(Prop_efficacy = Prop_efficacy_expCorr2,
                                                                Prop_futility = Prop_futility_expCorr2), 
                                                       by = 'Setting') %>%
  group_by(Setting,Iteration) %>% 
  mutate(p_lessthan_fut = ifelse(as.numeric(P_Upstrap) < P_futility, yes = 1, no = 0),
         p_lessthan_eff = ifelse(as.numeric(P_Upstrap) < P_efficacy, yes = 1, no = 0)) %>%
  summarise(proportion_sim_fut = mean(p_lessthan_fut),
            proportion_sim_eff = mean(p_lessthan_eff)) %>%
  merge(dat_setting_calibration_corrected %>% 
          select(Setting,P_futility,Prop_efficacy_expCorr2,P_efficacy,Prop_futility_expCorr2) %>%
          rename(Prop_efficacy = Prop_efficacy_expCorr2,
                 Prop_futility = Prop_futility_expCorr2), 
        by = 'Setting') %>%
  mutate(stop_futility = ifelse(proportion_sim_fut < Prop_futility, yes = 1, no = 0),
         stop_efficacy = ifelse(proportion_sim_eff > Prop_efficacy, yes = 1, no = 0))

#use stopping boundary results to get expected sample size and full trial result statistics
dat_calibrated_corrected_E2_upstraped_final_results <- dat_calibrated_corrected_E2_upstraped_results %>%
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

#combine data from all CU calibration corrections
dat_final_calibration_correction_CU <- rbind(
  dat_calibrated_corrected_M1_upstraped_final_results %>% mutate(Method = 'Calibrated Upstrap M1'),
  dat_calibrated_corrected_M2_upstraped_final_results %>% mutate(Method = 'Calibrated Upstrap M2'),
  dat_calibrated_corrected_M3_upstraped_final_results %>% mutate(Method = 'Calibrated Upstrap M3'),
  dat_calibrated_corrected_M4_upstraped_final_results %>% mutate(Method = 'Calibrated Upstrap M4'),
  dat_calibrated_corrected_M5_upstraped_final_results %>% mutate(Method = 'Calibrated Upstrap M5'),
  dat_calibrated_corrected_E1_upstraped_final_results %>% mutate(Method = 'Calibrated Upstrap E1'),
  dat_calibrated_corrected_E2_upstraped_final_results %>% mutate(Method = 'Calibrated Upstrap E2'),
  dat_final_output %>% filter(Method %in% c('OBrien Fleming GSD','Calibrated Upstrap'))
)

###############################################################
### SAVE OUTPUT DATASETS                                    ###
###############################################################

save(dat_setting_calibration_corrected, file = './01 Data/Data Processed/dat_setting_calibration_corrected.RData')
#load('./01 Data/Data Processed/dat_setting_calibration_corrected.RData')
save(dat_calibrated_corrected_M1_upstraped_results, file = './01 Data/Data Processed/dat_calibrated_corrected_M1_upstraped_results.RData')
#load('./01 Data/Data Processed/dat_calibrated_corrected_M1_upstraped_results.RData')
save(dat_calibrated_corrected_M1_upstraped_final_results, file = './01 Data/Data Processed/dat_calibrated_corrected_M1_upstraped_final_results.RData')
#load('./01 Data/Data Processed/dat_calibrated_corrected_M1_upstraped_final_results.RData')
save(dat_calibrated_corrected_M2_upstraped_results, file = './01 Data/Data Processed/dat_calibrated_corrected_M2_upstraped_results.RData')
#load('./01 Data/Data Processed/dat_calibrated_corrected_M2_upstraped_results.RData')
save(dat_calibrated_corrected_M2_upstraped_final_results, file = './01 Data/Data Processed/dat_calibrated_corrected_M2_upstraped_final_results.RData')
#load('./01 Data/Data Processed/dat_calibrated_corrected_M2_upstraped_final_results.RData')
save(dat_calibrated_corrected_M3_upstraped_results, file = './01 Data/Data Processed/dat_calibrated_corrected_M3_upstraped_results.RData')
#load('./01 Data/Data Processed/dat_calibrated_corrected_M3_upstraped_results.RData')
save(dat_calibrated_corrected_M3_upstraped_final_results, file = './01 Data/Data Processed/dat_calibrated_corrected_M3_upstraped_final_results.RData')
#load('./01 Data/Data Processed/dat_calibrated_corrected_M3_upstraped_final_results.RData')
save(dat_calibrated_corrected_M4_upstraped_results, file = './01 Data/Data Processed/dat_calibrated_corrected_M4_upstraped_results.RData')
#load('./01 Data/Data Processed/dat_calibrated_corrected_M4_upstraped_results.RData')
save(dat_calibrated_corrected_M4_upstraped_final_results, file = './01 Data/Data Processed/dat_calibrated_corrected_M4_upstraped_final_results.RData')
#load('./01 Data/Data Processed/dat_calibrated_corrected_M4_upstraped_final_results.RData')
save(dat_calibrated_corrected_M5_upstraped_results, file = './01 Data/Data Processed/dat_calibrated_corrected_M5_upstraped_results.RData')
#load('./01 Data/Data Processed/dat_calibrated_corrected_M5_upstraped_results.RData')
save(dat_calibrated_corrected_M5_upstraped_final_results, file = './01 Data/Data Processed/dat_calibrated_corrected_M5_upstraped_final_results.RData')
#load('./01 Data/Data Processed/dat_calibrated_corrected_M5_upstraped_final_results.RData')
save(dat_calibrated_corrected_E1_upstraped_results, file = './01 Data/Data Processed/dat_calibrated_corrected_E1_upstraped_results.RData')
#load('./01 Data/Data Processed/dat_calibrated_corrected_E1_upstraped_results.RData')
save(dat_calibrated_corrected_E1_upstraped_final_results, file = './01 Data/Data Processed/dat_calibrated_corrected_E1_upstraped_final_results.RData')
#load('./01 Data/Data Processed/dat_calibrated_corrected_E1_upstraped_final_results.RData')
save(dat_calibrated_corrected_E2_upstraped_results, file = './01 Data/Data Processed/dat_calibrated_corrected_E2_upstraped_results.RData')
#load('./01 Data/Data Processed/dat_calibrated_corrected_E2_upstraped_results.RData')
save(dat_calibrated_corrected_E2_upstraped_final_results, file = './01 Data/Data Processed/dat_calibrated_corrected_E2_upstraped_final_results.RData')
#load('./01 Data/Data Processed/dat_calibrated_corrected_E2_upstraped_final_results.RData')
save(dat_final_calibration_correction_CU, file = './01 Data/Data Processed/dat_final_calibration_correction_CU.RData')
#load('./01 Data/Data Processed/dat_final_calibration_correction_CU.RData')

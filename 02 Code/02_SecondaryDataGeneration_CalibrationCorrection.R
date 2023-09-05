###############################################################
### file: 02_SecondaryDataGeneration_GUCalibrationCorrection###
### authors: Jess Wild                                      ###
### creation date:    07/25/2023                            ###
### latest edit date: 07/25/2023                            ###
### description: code to apply multiple testing correction  ###
### to GU upstrapping calibration                           ###
###############################################################

#NOTE: code needs to be edited to split into the dat_results, dat_final_rests format instead of doing both in one step:
###############################################################
### APPLY UPSTRAPPING CRITERIA TO SIMULATION RESULTS        ###
###############################################################
###############################################################
### CALCULATE EVALUATION METRICS                            ###
###############################################################

###############################################################
### SETUP WORKSPACE                                         ###
###############################################################

#set working directory
setwd('C:/Users/wildje/Desktop/Dissertation/Paper 1')

#import packages
library(tidyverse)

#load input datasets
load('./01 Data/Data Processed/dat_setting.RData')
load('./01 Data/Data Processed/dat_results.RData')
load('./01 Data/Data Processed/dat_pow_err.RData')
load('./01 Data/Data Processed/dat_fullsample_reject.RData')
load('./01 Data/Data Processed/dat_setting_calibration.RData')
load('./01 Data/Data Processed/dat_AB_settings.RData')
load('./01 Data/Data Processed/dat_fullsample_reject_AB.RData')
load('./01 Data/Data Processed/dat_final_output.RData')
load('./01 Data/Data Processed/dat_final_calibration_no25.RData')

###############################################################
### APPLY CALIBRATION CORRECTIONS TO AU SIMULATION RESULTS  ###
###############################################################

#apply calibration corrections (multiplicative, exponential) to AU results
dat_arbitrary_corrected_M1_upstraped_results <- dat_results %>%
  group_by(Setting,Iteration) %>% 
  mutate(p_lessthan_fut = ifelse(as.numeric(P_Upstrap) < 0.05, yes = 1, no = 0),
         p_lessthan_eff = ifelse(as.numeric(P_Upstrap) < 0.05, yes = 1, no = 0)) %>%
  summarise(proportion_sim_fut = mean(p_lessthan_fut),
            proportion_sim_eff = mean(p_lessthan_eff)) %>%
  merge(dat_setting_calibration %>% select(Setting,Prop_futility,Prop_efficacy), by = 'Setting') %>%
  mutate(stop_futility = ifelse(proportion_sim_fut < 0.05/(4+1), yes = 1, no = 0),
         stop_efficacy = ifelse(proportion_sim_eff > 0.80+(1-0.80)/(4+1), yes = 1, no = 0)) %>%
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
dat_arbitrary_corrected_M2_upstraped_results <- dat_results %>%
  group_by(Setting,Iteration) %>% 
  mutate(p_lessthan_fut = ifelse(as.numeric(P_Upstrap) < 0.05, yes = 1, no = 0),
         p_lessthan_eff = ifelse(as.numeric(P_Upstrap) < 0.05, yes = 1, no = 0)) %>%
  summarise(proportion_sim_fut = mean(p_lessthan_fut),
            proportion_sim_eff = mean(p_lessthan_eff)) %>%
  merge(dat_setting_calibration %>% select(Setting,Prop_futility,Prop_efficacy), by = 'Setting') %>%
  mutate(stop_futility = ifelse(proportion_sim_fut < 0.05/(4), yes = 1, no = 0),
         stop_efficacy = ifelse(proportion_sim_eff > 0.80+(1-0.80)/(4), yes = 1, no = 0)) %>%
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
dat_arbitrary_corrected_M3_upstraped_results <- dat_results %>%
  group_by(Setting,Iteration) %>% 
  mutate(p_lessthan_fut = ifelse(as.numeric(P_Upstrap) < 0.05, yes = 1, no = 0),
         p_lessthan_eff = ifelse(as.numeric(P_Upstrap) < 0.05, yes = 1, no = 0)) %>%
  summarise(proportion_sim_fut = mean(p_lessthan_fut),
            proportion_sim_eff = mean(p_lessthan_eff)) %>%
  merge(dat_setting_calibration %>% select(Setting,Prop_futility,Prop_efficacy), by = 'Setting') %>%
  mutate(stop_futility = ifelse(proportion_sim_fut < 0.05/(4-1), yes = 1, no = 0),
         stop_efficacy = ifelse(proportion_sim_eff > 0.80+(1-0.80)/(4-1), yes = 1, no = 0)) %>%
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
dat_arbitrary_corrected_M4_upstraped_results <- dat_results %>%
  group_by(Setting,Iteration) %>% 
  mutate(p_lessthan_fut = ifelse(as.numeric(P_Upstrap) < 0.05, yes = 1, no = 0),
         p_lessthan_eff = ifelse(as.numeric(P_Upstrap) < 0.05, yes = 1, no = 0)) %>%
  summarise(proportion_sim_fut = mean(p_lessthan_fut),
            proportion_sim_eff = mean(p_lessthan_eff)) %>%
  merge(dat_setting_calibration %>% select(Setting,Prop_futility,Prop_efficacy), by = 'Setting') %>%
  mutate(stop_futility = ifelse(proportion_sim_fut < 0.05/(4-1.5), yes = 1, no = 0),
         stop_efficacy = ifelse(proportion_sim_eff > 0.80+(1-0.80)/(4-1.5), yes = 1, no = 0)) %>%
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
dat_arbitrary_corrected_M5_upstraped_results <- dat_results %>%
  group_by(Setting,Iteration) %>% 
  mutate(p_lessthan_fut = ifelse(as.numeric(P_Upstrap) < 0.05, yes = 1, no = 0),
         p_lessthan_eff = ifelse(as.numeric(P_Upstrap) < 0.05, yes = 1, no = 0)) %>%
  summarise(proportion_sim_fut = mean(p_lessthan_fut),
            proportion_sim_eff = mean(p_lessthan_eff)) %>%
  merge(dat_setting_calibration %>% select(Setting,Prop_futility,Prop_efficacy), by = 'Setting') %>%
  mutate(stop_futility = ifelse(proportion_sim_fut < 0.05/(4-2), yes = 1, no = 0),
         stop_efficacy = ifelse(proportion_sim_eff > 0.80+(1-0.80)/(4-2), yes = 1, no = 0)) %>%
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
dat_arbitrary_corrected_E1_upstraped_results <- dat_results %>%
  group_by(Setting,Iteration) %>% 
  mutate(p_lessthan_fut = ifelse(as.numeric(P_Upstrap) < 0.05, yes = 1, no = 0),
         p_lessthan_eff = ifelse(as.numeric(P_Upstrap) < 0.05, yes = 1, no = 0)) %>%
  summarise(proportion_sim_fut = mean(p_lessthan_fut),
            proportion_sim_eff = mean(p_lessthan_eff)) %>%
  merge(dat_setting_calibration %>% select(Setting,Prop_futility,Prop_efficacy), by = 'Setting') %>%
  mutate(stop_futility = ifelse(proportion_sim_fut < 0.05^0.95, yes = 1, no = 0),
         stop_efficacy = ifelse(proportion_sim_eff > 0.80^(1-0.95), yes = 1, no = 0)) %>%
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
dat_arbitrary_corrected_E2_upstraped_results <- dat_results %>%
  group_by(Setting,Iteration) %>% 
  mutate(p_lessthan_fut = ifelse(as.numeric(P_Upstrap) < 0.05, yes = 1, no = 0),
         p_lessthan_eff = ifelse(as.numeric(P_Upstrap) < 0.05, yes = 1, no = 0)) %>%
  summarise(proportion_sim_fut = mean(p_lessthan_fut),
            proportion_sim_eff = mean(p_lessthan_eff)) %>%
  merge(dat_setting_calibration %>% select(Setting,Prop_futility,Prop_efficacy), by = 'Setting') %>%
  mutate(stop_futility = ifelse(proportion_sim_fut < 0.05^0.70, yes = 1, no = 0),
         stop_efficacy = ifelse(proportion_sim_eff > 0.80^(1-0.70), yes = 1, no = 0)) %>%
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
### APPLY CALIBRATION CORRECTIONS TO CU SIMULATION RESULTS  ###
###############################################################

#apply calibration corrections (multiplicative, exponential) to CU results
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

#apply (corrected) calibrated upstrapping stopping boundaries to simulated results
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
         stop_efficacy = ifelse(proportion_sim_eff > Prop_efficacy, yes = 1, no = 0)) %>%
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
         stop_efficacy = ifelse(proportion_sim_eff > Prop_efficacy, yes = 1, no = 0)) %>%
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
         stop_efficacy = ifelse(proportion_sim_eff > Prop_efficacy, yes = 1, no = 0)) %>%
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
         stop_efficacy = ifelse(proportion_sim_eff > Prop_efficacy, yes = 1, no = 0)) %>%
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
         stop_efficacy = ifelse(proportion_sim_eff > Prop_efficacy, yes = 1, no = 0)) %>%
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
         stop_efficacy = ifelse(proportion_sim_eff > Prop_efficacy, yes = 1, no = 0)) %>%
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
         stop_efficacy = ifelse(proportion_sim_eff > Prop_efficacy, yes = 1, no = 0)) %>%
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
### APPLY CALIBRATION CORRECTIONS TO GU SIMULATION RESULTS  ###
###############################################################

#apply calibration corrections (multiplicative, exponential) to GU results
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

dat_ABcalibrated_corrected_M1_upstraped_results <- rbind(
  #n = 40
  merge(dat_results %>% select(Setting,Iteration,P_Upstrap), 
        dat_setting %>% select(Setting,n,interim_fraction,power) , by = 'Setting') %>%
    filter(n == 40) %>%
    merge(dat_AB_setting_calibration_corrected %>% filter(interim_fraction != 1.00 & n == 40), by = c('n','interim_fraction')) %>%
    group_by(Setting,Iteration,MonitoringType) %>% 
    mutate(p_lessthan_fut = ifelse(is.na(FutilityP) == F & (as.numeric(P_Upstrap) < FutilityP), yes = 1, no = 0),
           p_lessthan_eff = ifelse(is.na(EfficacyP) == F & (as.numeric(P_Upstrap) < EfficacyP), yes = 1, no = 0)) %>%
    summarise(proportion_sim_fut = mean(p_lessthan_fut),
              proportion_sim_eff = mean(p_lessthan_eff)) %>%
    merge(dat_setting %>% select(Setting,n,interim_fraction), by = 'Setting') %>%
    merge(dat_AB_setting_calibration_corrected, by = c('n','interim_fraction','MonitoringType')) %>%
    mutate(stop_futility = ifelse(proportion_sim_fut < FutilityProp_multCorr1, yes = 1, no = 0),
           stop_efficacy = ifelse(proportion_sim_eff > EfficacyProp_multCorr1, yes = 1, no = 0)) %>%
    select(Setting,Iteration,MonitoringType,stop_futility,stop_efficacy) %>%
    merge(dat_setting %>% filter(n == 40) %>% select(n,Setting,power,interim_fraction), by='Setting') %>%
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
    merge(dat_fullsample_reject_AB %>% filter(n == 40), by=c('n','power','MonitoringType','Iteration')) %>%
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
           'RejectionRate_FE_FS' = RejectionRate_FE_FS_FE),
  #n = 160
  merge(dat_results %>% select(Setting,Iteration,P_Upstrap), 
        dat_setting %>% select(Setting,n,interim_fraction,power) , by = 'Setting') %>%
    filter(n == 160) %>%
    merge(dat_AB_setting_calibration_corrected %>% filter(interim_fraction != 1.00 & n == 160), by = c('n','interim_fraction')) %>%
    group_by(Setting,Iteration,MonitoringType) %>% 
    mutate(p_lessthan_fut = ifelse(is.na(FutilityP) == F & (as.numeric(P_Upstrap) < FutilityP), yes = 1, no = 0),
           p_lessthan_eff = ifelse(is.na(EfficacyP) == F & (as.numeric(P_Upstrap) < EfficacyP), yes = 1, no = 0)) %>%
    summarise(proportion_sim_fut = mean(p_lessthan_fut),
              proportion_sim_eff = mean(p_lessthan_eff)) %>%
    merge(dat_setting %>% select(Setting,n,interim_fraction), by = 'Setting') %>%
    merge(dat_AB_setting_calibration_corrected, by = c('n','interim_fraction','MonitoringType')) %>%
    mutate(stop_futility = ifelse(proportion_sim_fut < FutilityProp_multCorr1, yes = 1, no = 0),
           stop_efficacy = ifelse(proportion_sim_eff > EfficacyProp_multCorr1, yes = 1, no = 0)) %>%
    select(Setting,Iteration,MonitoringType,stop_futility,stop_efficacy) %>%
    merge(dat_setting %>% filter(n == 160) %>% select(n,Setting,power,interim_fraction), by='Setting') %>%
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
    merge(dat_fullsample_reject_AB %>% filter(n == 160), by=c('n','power','MonitoringType','Iteration')) %>%
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
           'RejectionRate_FE_FS' = RejectionRate_FE_FS_FE),
  #n = 600
  merge(dat_results %>% select(Setting,Iteration,P_Upstrap), 
        dat_setting %>% select(Setting,n,interim_fraction,power) , by = 'Setting') %>%
    filter(n == 600) %>%
    merge(dat_AB_setting_calibration_corrected %>% filter(interim_fraction != 1.00 & n == 600), by = c('n','interim_fraction')) %>%
    group_by(Setting,Iteration,MonitoringType) %>% 
    mutate(p_lessthan_fut = ifelse(is.na(FutilityP) == F & (as.numeric(P_Upstrap) < FutilityP), yes = 1, no = 0),
           p_lessthan_eff = ifelse(is.na(EfficacyP) == F & (as.numeric(P_Upstrap) < EfficacyP), yes = 1, no = 0)) %>%
    summarise(proportion_sim_fut = mean(p_lessthan_fut),
              proportion_sim_eff = mean(p_lessthan_eff)) %>%
    merge(dat_setting %>% select(Setting,n,interim_fraction), by = 'Setting') %>%
    merge(dat_AB_setting_calibration_corrected, by = c('n','interim_fraction','MonitoringType')) %>%
    mutate(stop_futility = ifelse(proportion_sim_fut < FutilityProp_multCorr1, yes = 1, no = 0),
           stop_efficacy = ifelse(proportion_sim_eff > EfficacyProp_multCorr1, yes = 1, no = 0)) %>%
    select(Setting,Iteration,MonitoringType,stop_futility,stop_efficacy) %>%
    merge(dat_setting %>% filter(n == 600) %>% select(n,Setting,power,interim_fraction), by='Setting') %>%
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
    merge(dat_fullsample_reject_AB %>% filter(n == 600), by=c('n','power','MonitoringType','Iteration')) %>%
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
           'RejectionRate_FE_FS' = RejectionRate_FE_FS_FE),
  #n = 2000
  merge(dat_results %>% select(Setting,Iteration,P_Upstrap), 
        dat_setting %>% select(Setting,n,interim_fraction,power) , by = 'Setting') %>%
    filter(n == 2000) %>%
    merge(dat_AB_setting_calibration_corrected %>% filter(interim_fraction != 1.00 & n == 2000), by = c('n','interim_fraction')) %>%
    group_by(Setting,Iteration,MonitoringType) %>% 
    mutate(p_lessthan_fut = ifelse(is.na(FutilityP) == F & (as.numeric(P_Upstrap) < FutilityP), yes = 1, no = 0),
           p_lessthan_eff = ifelse(is.na(EfficacyP) == F & (as.numeric(P_Upstrap) < EfficacyP), yes = 1, no = 0)) %>%
    summarise(proportion_sim_fut = mean(p_lessthan_fut),
              proportion_sim_eff = mean(p_lessthan_eff)) %>%
    merge(dat_setting %>% select(Setting,n,interim_fraction), by = 'Setting') %>%
    merge(dat_AB_setting_calibration_corrected, by = c('n','interim_fraction','MonitoringType')) %>%
    mutate(stop_futility = ifelse(proportion_sim_fut < FutilityProp_multCorr1, yes = 1, no = 0),
           stop_efficacy = ifelse(proportion_sim_eff > EfficacyProp_multCorr1, yes = 1, no = 0)) %>%
    select(Setting,Iteration,MonitoringType,stop_futility,stop_efficacy) %>%
    merge(dat_setting %>% filter(n == 2000) %>% select(n,Setting,power,interim_fraction), by='Setting') %>%
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
    merge(dat_fullsample_reject_AB %>% filter(n == 2000), by=c('n','power','MonitoringType','Iteration')) %>%
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
)
dat_ABcalibrated_corrected_M2_upstraped_results <- rbind(
  #n = 40
  merge(dat_results %>% select(Setting,Iteration,P_Upstrap), 
        dat_setting %>% select(Setting,n,interim_fraction,power) , by = 'Setting') %>%
    filter(n == 40) %>%
    merge(dat_AB_setting_calibration_corrected %>% filter(interim_fraction != 1.00 & n == 40), by = c('n','interim_fraction')) %>%
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
    merge(dat_setting %>% filter(n == 40) %>% select(n,Setting,power,interim_fraction), by='Setting') %>%
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
    merge(dat_fullsample_reject_AB %>% filter(n == 40), by=c('n','power','MonitoringType','Iteration')) %>%
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
           'RejectionRate_FE_FS' = RejectionRate_FE_FS_FE),
  #n = 160
  merge(dat_results %>% select(Setting,Iteration,P_Upstrap), 
        dat_setting %>% select(Setting,n,interim_fraction,power) , by = 'Setting') %>%
    filter(n == 160) %>%
    merge(dat_AB_setting_calibration_corrected %>% filter(interim_fraction != 1.00 & n == 160), by = c('n','interim_fraction')) %>%
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
    merge(dat_setting %>% filter(n == 160) %>% select(n,Setting,power,interim_fraction), by='Setting') %>%
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
           Decision_FE_E = ifelse(Decision_FE==1&ExpectedN_EO<ExpectedN_FO,yes = 1, no = 0)) %>%
    merge(dat_fullsample_reject_AB %>% filter(n == 160), by=c('n','power','MonitoringType','Iteration')) %>%
    mutate(RejectionRate_FO_FS = ifelse(Decision_FO==0&Reject_FS==1,yes=1,no=0),
           RejectionRate_EO_FS = ifelse(Decision_EO==1|Reject_FS==1,yes=1,no=0),
           RejectionRate_FE_FS = ifelse(Decision_FE_E==1|(Decision_FE_F==0&Reject_FS==1),yes=1,no=0)) %>%
    select(n,power,MonitoringType,Iteration,ExpectedN_FO,ExpectedN_EO,ExpectedN_FE,Decision_FO,Decision_EO,Decision_FE,Decision_FE_F,Decision_FE_E,RejectionRate_FO_FS,RejectionRate_EO_FS,RejectionRate_FE_FS) %>%
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
              RejectionRate_FO_FS = mean(RejectionRate_FO_FS),
              RejectionRate_EO_FS = mean(RejectionRate_EO_FS),
              RejectionRate_FE_FS = mean(RejectionRate_FE_FS)) %>%
    mutate(MonitoringType = factor(MonitoringType)) %>%
    pivot_wider(
      id_cols = c('n','power'),
      names_from = c('MonitoringType'),
      values_from = c('ExpectedN_FO_mean','ExpectedN_FO_sd','ExpectedN_EO_mean','ExpectedN_EO_sd',
                      'ExpectedN_FE_mean','ExpectedN_FE_sd','DecisionFO','DecisionEO','DecisionFE',
                      'DecisionFE_F','DecisionFE_E','RejectionRate_FO_FS','RejectionRate_EO_FS',
                      'RejectionRate_FE_FS')
    ) %>%
    select(n,power,ExpectedN_FO_mean_FO,ExpectedN_FO_sd_FO,
           ExpectedN_EO_mean_EO,ExpectedN_EO_sd_EO,
           ExpectedN_FE_mean_FE,ExpectedN_FE_sd_FE,
           DecisionFO_FO,DecisionEO_EO,DecisionFE_FE,DecisionFE_F_FE,DecisionFE_E_FE,
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
           'RejectionRate_FO_FS' = RejectionRate_FO_FS_FO,
           'RejectionRate_EO_FS' = RejectionRate_EO_FS_EO,
           'RejectionRate_FE_FS' = RejectionRate_FE_FS_FE),
  #n = 600
  merge(dat_results %>% select(Setting,Iteration,P_Upstrap), 
        dat_setting %>% select(Setting,n,interim_fraction,power) , by = 'Setting') %>%
    filter(n == 600) %>%
    merge(dat_AB_setting_calibration_corrected %>% filter(interim_fraction != 1.00 & n == 600), by = c('n','interim_fraction')) %>%
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
    merge(dat_setting %>% filter(n == 600) %>% select(n,Setting,power,interim_fraction), by='Setting') %>%
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
    merge(dat_fullsample_reject_AB %>% filter(n == 600), by=c('n','power','MonitoringType','Iteration')) %>%
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
           'RejectionRate_FE_FS' = RejectionRate_FE_FS_FE),
  #n = 2000
  merge(dat_results %>% select(Setting,Iteration,P_Upstrap), 
        dat_setting %>% select(Setting,n,interim_fraction,power) , by = 'Setting') %>%
    filter(n == 2000) %>%
    merge(dat_AB_setting_calibration_corrected %>% filter(interim_fraction != 1.00 & n == 2000), by = c('n','interim_fraction')) %>%
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
    merge(dat_setting %>% filter(n == 2000) %>% select(n,Setting,power,interim_fraction), by='Setting') %>%
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
    merge(dat_fullsample_reject_AB %>% filter(n == 2000), by=c('n','power','MonitoringType','Iteration')) %>%
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
)
dat_ABcalibrated_corrected_M3_upstraped_results <- rbind(
  #n = 40
  merge(dat_results %>% select(Setting,Iteration,P_Upstrap), 
        dat_setting %>% select(Setting,n,interim_fraction,power) , by = 'Setting') %>%
    filter(n == 40) %>%
    merge(dat_AB_setting_calibration_corrected %>% filter(interim_fraction != 1.00 & n == 40), by = c('n','interim_fraction')) %>%
    group_by(Setting,Iteration,MonitoringType) %>% 
    mutate(p_lessthan_fut = ifelse(is.na(FutilityP) == F & (as.numeric(P_Upstrap) < FutilityP), yes = 1, no = 0),
           p_lessthan_eff = ifelse(is.na(EfficacyP) == F & (as.numeric(P_Upstrap) < EfficacyP), yes = 1, no = 0)) %>%
    summarise(proportion_sim_fut = mean(p_lessthan_fut),
              proportion_sim_eff = mean(p_lessthan_eff)) %>%
    merge(dat_setting %>% select(Setting,n,interim_fraction), by = 'Setting') %>%
    merge(dat_AB_setting_calibration_corrected, by = c('n','interim_fraction','MonitoringType')) %>%
    mutate(stop_futility = ifelse(proportion_sim_fut < FutilityProp_multCorr3, yes = 1, no = 0),
           stop_efficacy = ifelse(proportion_sim_eff > EfficacyProp_multCorr3, yes = 1, no = 0)) %>%
    select(Setting,Iteration,MonitoringType,stop_futility,stop_efficacy) %>%
    merge(dat_setting %>% filter(n == 40) %>% select(n,Setting,power,interim_fraction), by='Setting') %>%
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
    merge(dat_fullsample_reject_AB %>% filter(n == 40), by=c('n','power','MonitoringType','Iteration')) %>%
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
           'RejectionRate_FE_FS' = RejectionRate_FE_FS_FE),
  #n = 160
  merge(dat_results %>% select(Setting,Iteration,P_Upstrap), 
        dat_setting %>% select(Setting,n,interim_fraction,power) , by = 'Setting') %>%
    filter(n == 160) %>%
    merge(dat_AB_setting_calibration_corrected %>% filter(interim_fraction != 1.00 & n == 160), by = c('n','interim_fraction')) %>%
    group_by(Setting,Iteration,MonitoringType) %>% 
    mutate(p_lessthan_fut = ifelse(is.na(FutilityP) == F & (as.numeric(P_Upstrap) < FutilityP), yes = 1, no = 0),
           p_lessthan_eff = ifelse(is.na(EfficacyP) == F & (as.numeric(P_Upstrap) < EfficacyP), yes = 1, no = 0)) %>%
    summarise(proportion_sim_fut = mean(p_lessthan_fut),
              proportion_sim_eff = mean(p_lessthan_eff)) %>%
    merge(dat_setting %>% select(Setting,n,interim_fraction), by = 'Setting') %>%
    merge(dat_AB_setting_calibration_corrected, by = c('n','interim_fraction','MonitoringType')) %>%
    mutate(stop_futility = ifelse(proportion_sim_fut < FutilityProp_multCorr3, yes = 1, no = 0),
           stop_efficacy = ifelse(proportion_sim_eff > EfficacyProp_multCorr3, yes = 1, no = 0)) %>%
    select(Setting,Iteration,MonitoringType,stop_futility,stop_efficacy) %>%
    merge(dat_setting %>% filter(n == 160) %>% select(n,Setting,power,interim_fraction), by='Setting') %>%
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
           Decision_FE_E = ifelse(Decision_FE==1&ExpectedN_EO<ExpectedN_FO,yes = 1, no = 0)) %>%
    merge(dat_fullsample_reject_AB %>% filter(n == 160), by=c('n','power','MonitoringType','Iteration')) %>%
    mutate(RejectionRate_FO_FS = ifelse(Decision_FO==0&Reject_FS==1,yes=1,no=0),
           RejectionRate_EO_FS = ifelse(Decision_EO==1|Reject_FS==1,yes=1,no=0),
           RejectionRate_FE_FS = ifelse(Decision_FE_E==1|(Decision_FE_F==0&Reject_FS==1),yes=1,no=0)) %>%
    select(n,power,MonitoringType,Iteration,ExpectedN_FO,ExpectedN_EO,ExpectedN_FE,Decision_FO,Decision_EO,Decision_FE,Decision_FE_F,Decision_FE_E,RejectionRate_FO_FS,RejectionRate_EO_FS,RejectionRate_FE_FS) %>%
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
              RejectionRate_FO_FS = mean(RejectionRate_FO_FS),
              RejectionRate_EO_FS = mean(RejectionRate_EO_FS),
              RejectionRate_FE_FS = mean(RejectionRate_FE_FS)) %>%
    mutate(MonitoringType = factor(MonitoringType)) %>%
    pivot_wider(
      id_cols = c('n','power'),
      names_from = c('MonitoringType'),
      values_from = c('ExpectedN_FO_mean','ExpectedN_FO_sd','ExpectedN_EO_mean','ExpectedN_EO_sd',
                      'ExpectedN_FE_mean','ExpectedN_FE_sd','DecisionFO','DecisionEO','DecisionFE',
                      'DecisionFE_F','DecisionFE_E','RejectionRate_FO_FS','RejectionRate_EO_FS',
                      'RejectionRate_FE_FS')
    ) %>%
    select(n,power,ExpectedN_FO_mean_FO,ExpectedN_FO_sd_FO,
           ExpectedN_EO_mean_EO,ExpectedN_EO_sd_EO,
           ExpectedN_FE_mean_FE,ExpectedN_FE_sd_FE,
           DecisionFO_FO,DecisionEO_EO,DecisionFE_FE,DecisionFE_F_FE,DecisionFE_E_FE,
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
           'RejectionRate_FO_FS' = RejectionRate_FO_FS_FO,
           'RejectionRate_EO_FS' = RejectionRate_EO_FS_EO,
           'RejectionRate_FE_FS' = RejectionRate_FE_FS_FE),
  #n = 600
  merge(dat_results %>% select(Setting,Iteration,P_Upstrap), 
        dat_setting %>% select(Setting,n,interim_fraction,power) , by = 'Setting') %>%
    filter(n == 600) %>%
    merge(dat_AB_setting_calibration_corrected %>% filter(interim_fraction != 1.00 & n == 600), by = c('n','interim_fraction')) %>%
    group_by(Setting,Iteration,MonitoringType) %>% 
    mutate(p_lessthan_fut = ifelse(is.na(FutilityP) == F & (as.numeric(P_Upstrap) < FutilityP), yes = 1, no = 0),
           p_lessthan_eff = ifelse(is.na(EfficacyP) == F & (as.numeric(P_Upstrap) < EfficacyP), yes = 1, no = 0)) %>%
    summarise(proportion_sim_fut = mean(p_lessthan_fut),
              proportion_sim_eff = mean(p_lessthan_eff)) %>%
    merge(dat_setting %>% select(Setting,n,interim_fraction), by = 'Setting') %>%
    merge(dat_AB_setting_calibration_corrected, by = c('n','interim_fraction','MonitoringType')) %>%
    mutate(stop_futility = ifelse(proportion_sim_fut < FutilityProp_multCorr3, yes = 1, no = 0),
           stop_efficacy = ifelse(proportion_sim_eff > EfficacyProp_multCorr3, yes = 1, no = 0)) %>%
    select(Setting,Iteration,MonitoringType,stop_futility,stop_efficacy) %>%
    merge(dat_setting %>% filter(n == 600) %>% select(n,Setting,power,interim_fraction), by='Setting') %>%
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
    merge(dat_fullsample_reject_AB %>% filter(n == 600), by=c('n','power','MonitoringType','Iteration')) %>%
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
           'RejectionRate_FE_FS' = RejectionRate_FE_FS_FE),
  #n = 2000
  merge(dat_results %>% select(Setting,Iteration,P_Upstrap), 
        dat_setting %>% select(Setting,n,interim_fraction,power) , by = 'Setting') %>%
    filter(n == 2000) %>%
    merge(dat_AB_setting_calibration_corrected %>% filter(interim_fraction != 1.00 & n == 2000), by = c('n','interim_fraction')) %>%
    group_by(Setting,Iteration,MonitoringType) %>% 
    mutate(p_lessthan_fut = ifelse(is.na(FutilityP) == F & (as.numeric(P_Upstrap) < FutilityP), yes = 1, no = 0),
           p_lessthan_eff = ifelse(is.na(EfficacyP) == F & (as.numeric(P_Upstrap) < EfficacyP), yes = 1, no = 0)) %>%
    summarise(proportion_sim_fut = mean(p_lessthan_fut),
              proportion_sim_eff = mean(p_lessthan_eff)) %>%
    merge(dat_setting %>% select(Setting,n,interim_fraction), by = 'Setting') %>%
    merge(dat_AB_setting_calibration_corrected, by = c('n','interim_fraction','MonitoringType')) %>%
    mutate(stop_futility = ifelse(proportion_sim_fut < FutilityProp_multCorr3, yes = 1, no = 0),
           stop_efficacy = ifelse(proportion_sim_eff > EfficacyProp_multCorr3, yes = 1, no = 0)) %>%
    select(Setting,Iteration,MonitoringType,stop_futility,stop_efficacy) %>%
    merge(dat_setting %>% filter(n == 2000) %>% select(n,Setting,power,interim_fraction), by='Setting') %>%
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
    merge(dat_fullsample_reject_AB %>% filter(n == 2000), by=c('n','power','MonitoringType','Iteration')) %>%
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
)
dat_ABcalibrated_corrected_M4_upstraped_results <- rbind(
  #n = 40
  merge(dat_results %>% select(Setting,Iteration,P_Upstrap), 
        dat_setting %>% select(Setting,n,interim_fraction,power) , by = 'Setting') %>%
    filter(n == 40) %>%
    merge(dat_AB_setting_calibration_corrected %>% filter(interim_fraction != 1.00 & n == 40), by = c('n','interim_fraction')) %>%
    group_by(Setting,Iteration,MonitoringType) %>% 
    mutate(p_lessthan_fut = ifelse(is.na(FutilityP) == F & (as.numeric(P_Upstrap) < FutilityP), yes = 1, no = 0),
           p_lessthan_eff = ifelse(is.na(EfficacyP) == F & (as.numeric(P_Upstrap) < EfficacyP), yes = 1, no = 0)) %>%
    summarise(proportion_sim_fut = mean(p_lessthan_fut),
              proportion_sim_eff = mean(p_lessthan_eff)) %>%
    merge(dat_setting %>% select(Setting,n,interim_fraction), by = 'Setting') %>%
    merge(dat_AB_setting_calibration_corrected, by = c('n','interim_fraction','MonitoringType')) %>%
    mutate(stop_futility = ifelse(proportion_sim_fut < FutilityProp_multCorr4, yes = 1, no = 0),
           stop_efficacy = ifelse(proportion_sim_eff > EfficacyProp_multCorr4, yes = 1, no = 0)) %>%
    select(Setting,Iteration,MonitoringType,stop_futility,stop_efficacy) %>%
    merge(dat_setting %>% filter(n == 40) %>% select(n,Setting,power,interim_fraction), by='Setting') %>%
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
    merge(dat_fullsample_reject_AB %>% filter(n == 40), by=c('n','power','MonitoringType','Iteration')) %>%
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
           'RejectionRate_FE_FS' = RejectionRate_FE_FS_FE),
  #n = 160
  merge(dat_results %>% select(Setting,Iteration,P_Upstrap), 
        dat_setting %>% select(Setting,n,interim_fraction,power) , by = 'Setting') %>%
    filter(n == 160) %>%
    merge(dat_AB_setting_calibration_corrected %>% filter(interim_fraction != 1.00 & n == 160), by = c('n','interim_fraction')) %>%
    group_by(Setting,Iteration,MonitoringType) %>% 
    mutate(p_lessthan_fut = ifelse(is.na(FutilityP) == F & (as.numeric(P_Upstrap) < FutilityP), yes = 1, no = 0),
           p_lessthan_eff = ifelse(is.na(EfficacyP) == F & (as.numeric(P_Upstrap) < EfficacyP), yes = 1, no = 0)) %>%
    summarise(proportion_sim_fut = mean(p_lessthan_fut),
              proportion_sim_eff = mean(p_lessthan_eff)) %>%
    merge(dat_setting %>% select(Setting,n,interim_fraction), by = 'Setting') %>%
    merge(dat_AB_setting_calibration_corrected, by = c('n','interim_fraction','MonitoringType')) %>%
    mutate(stop_futility = ifelse(proportion_sim_fut < FutilityProp_multCorr4, yes = 1, no = 0),
           stop_efficacy = ifelse(proportion_sim_eff > EfficacyProp_multCorr4, yes = 1, no = 0)) %>%
    select(Setting,Iteration,MonitoringType,stop_futility,stop_efficacy) %>%
    merge(dat_setting %>% filter(n == 160) %>% select(n,Setting,power,interim_fraction), by='Setting') %>%
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
           Decision_FE_E = ifelse(Decision_FE==1&ExpectedN_EO<ExpectedN_FO,yes = 1, no = 0)) %>%
    merge(dat_fullsample_reject_AB %>% filter(n == 160), by=c('n','power','MonitoringType','Iteration')) %>%
    mutate(RejectionRate_FO_FS = ifelse(Decision_FO==0&Reject_FS==1,yes=1,no=0),
           RejectionRate_EO_FS = ifelse(Decision_EO==1|Reject_FS==1,yes=1,no=0),
           RejectionRate_FE_FS = ifelse(Decision_FE_E==1|(Decision_FE_F==0&Reject_FS==1),yes=1,no=0)) %>%
    select(n,power,MonitoringType,Iteration,ExpectedN_FO,ExpectedN_EO,ExpectedN_FE,Decision_FO,Decision_EO,Decision_FE,Decision_FE_F,Decision_FE_E,RejectionRate_FO_FS,RejectionRate_EO_FS,RejectionRate_FE_FS) %>%
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
              RejectionRate_FO_FS = mean(RejectionRate_FO_FS),
              RejectionRate_EO_FS = mean(RejectionRate_EO_FS),
              RejectionRate_FE_FS = mean(RejectionRate_FE_FS)) %>%
    mutate(MonitoringType = factor(MonitoringType)) %>%
    pivot_wider(
      id_cols = c('n','power'),
      names_from = c('MonitoringType'),
      values_from = c('ExpectedN_FO_mean','ExpectedN_FO_sd','ExpectedN_EO_mean','ExpectedN_EO_sd',
                      'ExpectedN_FE_mean','ExpectedN_FE_sd','DecisionFO','DecisionEO','DecisionFE',
                      'DecisionFE_F','DecisionFE_E','RejectionRate_FO_FS','RejectionRate_EO_FS',
                      'RejectionRate_FE_FS')
    ) %>%
    select(n,power,ExpectedN_FO_mean_FO,ExpectedN_FO_sd_FO,
           ExpectedN_EO_mean_EO,ExpectedN_EO_sd_EO,
           ExpectedN_FE_mean_FE,ExpectedN_FE_sd_FE,
           DecisionFO_FO,DecisionEO_EO,DecisionFE_FE,DecisionFE_F_FE,DecisionFE_E_FE,
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
           'RejectionRate_FO_FS' = RejectionRate_FO_FS_FO,
           'RejectionRate_EO_FS' = RejectionRate_EO_FS_EO,
           'RejectionRate_FE_FS' = RejectionRate_FE_FS_FE),
  #n = 600
  merge(dat_results %>% select(Setting,Iteration,P_Upstrap), 
        dat_setting %>% select(Setting,n,interim_fraction,power) , by = 'Setting') %>%
    filter(n == 600) %>%
    merge(dat_AB_setting_calibration_corrected %>% filter(interim_fraction != 1.00 & n == 600), by = c('n','interim_fraction')) %>%
    group_by(Setting,Iteration,MonitoringType) %>% 
    mutate(p_lessthan_fut = ifelse(is.na(FutilityP) == F & (as.numeric(P_Upstrap) < FutilityP), yes = 1, no = 0),
           p_lessthan_eff = ifelse(is.na(EfficacyP) == F & (as.numeric(P_Upstrap) < EfficacyP), yes = 1, no = 0)) %>%
    summarise(proportion_sim_fut = mean(p_lessthan_fut),
              proportion_sim_eff = mean(p_lessthan_eff)) %>%
    merge(dat_setting %>% select(Setting,n,interim_fraction), by = 'Setting') %>%
    merge(dat_AB_setting_calibration_corrected, by = c('n','interim_fraction','MonitoringType')) %>%
    mutate(stop_futility = ifelse(proportion_sim_fut < FutilityProp_multCorr4, yes = 1, no = 0),
           stop_efficacy = ifelse(proportion_sim_eff > EfficacyProp_multCorr4, yes = 1, no = 0)) %>%
    select(Setting,Iteration,MonitoringType,stop_futility,stop_efficacy) %>%
    merge(dat_setting %>% filter(n == 600) %>% select(n,Setting,power,interim_fraction), by='Setting') %>%
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
    merge(dat_fullsample_reject_AB %>% filter(n == 600), by=c('n','power','MonitoringType','Iteration')) %>%
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
           'RejectionRate_FE_FS' = RejectionRate_FE_FS_FE),
  #n = 2000
  merge(dat_results %>% select(Setting,Iteration,P_Upstrap), 
        dat_setting %>% select(Setting,n,interim_fraction,power) , by = 'Setting') %>%
    filter(n == 2000) %>%
    merge(dat_AB_setting_calibration_corrected %>% filter(interim_fraction != 1.00 & n == 2000), by = c('n','interim_fraction')) %>%
    group_by(Setting,Iteration,MonitoringType) %>% 
    mutate(p_lessthan_fut = ifelse(is.na(FutilityP) == F & (as.numeric(P_Upstrap) < FutilityP), yes = 1, no = 0),
           p_lessthan_eff = ifelse(is.na(EfficacyP) == F & (as.numeric(P_Upstrap) < EfficacyP), yes = 1, no = 0)) %>%
    summarise(proportion_sim_fut = mean(p_lessthan_fut),
              proportion_sim_eff = mean(p_lessthan_eff)) %>%
    merge(dat_setting %>% select(Setting,n,interim_fraction), by = 'Setting') %>%
    merge(dat_AB_setting_calibration_corrected, by = c('n','interim_fraction','MonitoringType')) %>%
    mutate(stop_futility = ifelse(proportion_sim_fut < FutilityProp_multCorr4, yes = 1, no = 0),
           stop_efficacy = ifelse(proportion_sim_eff > EfficacyProp_multCorr4, yes = 1, no = 0)) %>%
    select(Setting,Iteration,MonitoringType,stop_futility,stop_efficacy) %>%
    merge(dat_setting %>% filter(n == 2000) %>% select(n,Setting,power,interim_fraction), by='Setting') %>%
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
    merge(dat_fullsample_reject_AB %>% filter(n == 2000), by=c('n','power','MonitoringType','Iteration')) %>%
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
)
dat_ABcalibrated_corrected_M5_upstraped_results <- rbind(
  #n = 40
  merge(dat_results %>% select(Setting,Iteration,P_Upstrap), 
        dat_setting %>% select(Setting,n,interim_fraction,power) , by = 'Setting') %>%
    filter(n == 40) %>%
    merge(dat_AB_setting_calibration_corrected %>% filter(interim_fraction != 1.00 & n == 40), by = c('n','interim_fraction')) %>%
    group_by(Setting,Iteration,MonitoringType) %>% 
    mutate(p_lessthan_fut = ifelse(is.na(FutilityP) == F & (as.numeric(P_Upstrap) < FutilityP), yes = 1, no = 0),
           p_lessthan_eff = ifelse(is.na(EfficacyP) == F & (as.numeric(P_Upstrap) < EfficacyP), yes = 1, no = 0)) %>%
    summarise(proportion_sim_fut = mean(p_lessthan_fut),
              proportion_sim_eff = mean(p_lessthan_eff)) %>%
    merge(dat_setting %>% select(Setting,n,interim_fraction), by = 'Setting') %>%
    merge(dat_AB_setting_calibration_corrected, by = c('n','interim_fraction','MonitoringType')) %>%
    mutate(stop_futility = ifelse(proportion_sim_fut < FutilityProp_multCorr5, yes = 1, no = 0),
           stop_efficacy = ifelse(proportion_sim_eff > EfficacyProp_multCorr5, yes = 1, no = 0)) %>%
    select(Setting,Iteration,MonitoringType,stop_futility,stop_efficacy) %>%
    merge(dat_setting %>% filter(n == 40) %>% select(n,Setting,power,interim_fraction), by='Setting') %>%
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
    merge(dat_fullsample_reject_AB %>% filter(n == 40), by=c('n','power','MonitoringType','Iteration')) %>%
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
           'RejectionRate_FE_FS' = RejectionRate_FE_FS_FE),
  #n = 160
  merge(dat_results %>% select(Setting,Iteration,P_Upstrap), 
        dat_setting %>% select(Setting,n,interim_fraction,power) , by = 'Setting') %>%
    filter(n == 160) %>%
    merge(dat_AB_setting_calibration_corrected %>% filter(interim_fraction != 1.00 & n == 160), by = c('n','interim_fraction')) %>%
    group_by(Setting,Iteration,MonitoringType) %>% 
    mutate(p_lessthan_fut = ifelse(is.na(FutilityP) == F & (as.numeric(P_Upstrap) < FutilityP), yes = 1, no = 0),
           p_lessthan_eff = ifelse(is.na(EfficacyP) == F & (as.numeric(P_Upstrap) < EfficacyP), yes = 1, no = 0)) %>%
    summarise(proportion_sim_fut = mean(p_lessthan_fut),
              proportion_sim_eff = mean(p_lessthan_eff)) %>%
    merge(dat_setting %>% select(Setting,n,interim_fraction), by = 'Setting') %>%
    merge(dat_AB_setting_calibration_corrected, by = c('n','interim_fraction','MonitoringType')) %>%
    mutate(stop_futility = ifelse(proportion_sim_fut < FutilityProp_multCorr5, yes = 1, no = 0),
           stop_efficacy = ifelse(proportion_sim_eff > EfficacyProp_multCorr5, yes = 1, no = 0)) %>%
    select(Setting,Iteration,MonitoringType,stop_futility,stop_efficacy) %>%
    merge(dat_setting %>% filter(n == 160) %>% select(n,Setting,power,interim_fraction), by='Setting') %>%
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
           Decision_FE_E = ifelse(Decision_FE==1&ExpectedN_EO<ExpectedN_FO,yes = 1, no = 0)) %>%
    merge(dat_fullsample_reject_AB %>% filter(n == 160), by=c('n','power','MonitoringType','Iteration')) %>%
    mutate(RejectionRate_FO_FS = ifelse(Decision_FO==0&Reject_FS==1,yes=1,no=0),
           RejectionRate_EO_FS = ifelse(Decision_EO==1|Reject_FS==1,yes=1,no=0),
           RejectionRate_FE_FS = ifelse(Decision_FE_E==1|(Decision_FE_F==0&Reject_FS==1),yes=1,no=0)) %>%
    select(n,power,MonitoringType,Iteration,ExpectedN_FO,ExpectedN_EO,ExpectedN_FE,Decision_FO,Decision_EO,Decision_FE,Decision_FE_F,Decision_FE_E,RejectionRate_FO_FS,RejectionRate_EO_FS,RejectionRate_FE_FS) %>%
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
              RejectionRate_FO_FS = mean(RejectionRate_FO_FS),
              RejectionRate_EO_FS = mean(RejectionRate_EO_FS),
              RejectionRate_FE_FS = mean(RejectionRate_FE_FS)) %>%
    mutate(MonitoringType = factor(MonitoringType)) %>%
    pivot_wider(
      id_cols = c('n','power'),
      names_from = c('MonitoringType'),
      values_from = c('ExpectedN_FO_mean','ExpectedN_FO_sd','ExpectedN_EO_mean','ExpectedN_EO_sd',
                      'ExpectedN_FE_mean','ExpectedN_FE_sd','DecisionFO','DecisionEO','DecisionFE',
                      'DecisionFE_F','DecisionFE_E','RejectionRate_FO_FS','RejectionRate_EO_FS',
                      'RejectionRate_FE_FS')
    ) %>%
    select(n,power,ExpectedN_FO_mean_FO,ExpectedN_FO_sd_FO,
           ExpectedN_EO_mean_EO,ExpectedN_EO_sd_EO,
           ExpectedN_FE_mean_FE,ExpectedN_FE_sd_FE,
           DecisionFO_FO,DecisionEO_EO,DecisionFE_FE,DecisionFE_F_FE,DecisionFE_E_FE,
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
           'RejectionRate_FO_FS' = RejectionRate_FO_FS_FO,
           'RejectionRate_EO_FS' = RejectionRate_EO_FS_EO,
           'RejectionRate_FE_FS' = RejectionRate_FE_FS_FE),
  #n = 600
  merge(dat_results %>% select(Setting,Iteration,P_Upstrap), 
        dat_setting %>% select(Setting,n,interim_fraction,power) , by = 'Setting') %>%
    filter(n == 600) %>%
    merge(dat_AB_setting_calibration_corrected %>% filter(interim_fraction != 1.00 & n == 600), by = c('n','interim_fraction')) %>%
    group_by(Setting,Iteration,MonitoringType) %>% 
    mutate(p_lessthan_fut = ifelse(is.na(FutilityP) == F & (as.numeric(P_Upstrap) < FutilityP), yes = 1, no = 0),
           p_lessthan_eff = ifelse(is.na(EfficacyP) == F & (as.numeric(P_Upstrap) < EfficacyP), yes = 1, no = 0)) %>%
    summarise(proportion_sim_fut = mean(p_lessthan_fut),
              proportion_sim_eff = mean(p_lessthan_eff)) %>%
    merge(dat_setting %>% select(Setting,n,interim_fraction), by = 'Setting') %>%
    merge(dat_AB_setting_calibration_corrected, by = c('n','interim_fraction','MonitoringType')) %>%
    mutate(stop_futility = ifelse(proportion_sim_fut < FutilityProp_multCorr5, yes = 1, no = 0),
           stop_efficacy = ifelse(proportion_sim_eff > EfficacyProp_multCorr5, yes = 1, no = 0)) %>%
    select(Setting,Iteration,MonitoringType,stop_futility,stop_efficacy) %>%
    merge(dat_setting %>% filter(n == 600) %>% select(n,Setting,power,interim_fraction), by='Setting') %>%
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
    merge(dat_fullsample_reject_AB %>% filter(n == 600), by=c('n','power','MonitoringType','Iteration')) %>%
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
           'RejectionRate_FE_FS' = RejectionRate_FE_FS_FE),
  #n = 2000
  merge(dat_results %>% select(Setting,Iteration,P_Upstrap), 
        dat_setting %>% select(Setting,n,interim_fraction,power) , by = 'Setting') %>%
    filter(n == 2000) %>%
    merge(dat_AB_setting_calibration_corrected %>% filter(interim_fraction != 1.00 & n == 2000), by = c('n','interim_fraction')) %>%
    group_by(Setting,Iteration,MonitoringType) %>% 
    mutate(p_lessthan_fut = ifelse(is.na(FutilityP) == F & (as.numeric(P_Upstrap) < FutilityP), yes = 1, no = 0),
           p_lessthan_eff = ifelse(is.na(EfficacyP) == F & (as.numeric(P_Upstrap) < EfficacyP), yes = 1, no = 0)) %>%
    summarise(proportion_sim_fut = mean(p_lessthan_fut),
              proportion_sim_eff = mean(p_lessthan_eff)) %>%
    merge(dat_setting %>% select(Setting,n,interim_fraction), by = 'Setting') %>%
    merge(dat_AB_setting_calibration_corrected, by = c('n','interim_fraction','MonitoringType')) %>%
    mutate(stop_futility = ifelse(proportion_sim_fut < FutilityProp_multCorr5, yes = 1, no = 0),
           stop_efficacy = ifelse(proportion_sim_eff > EfficacyProp_multCorr5, yes = 1, no = 0)) %>%
    select(Setting,Iteration,MonitoringType,stop_futility,stop_efficacy) %>%
    merge(dat_setting %>% filter(n == 2000) %>% select(n,Setting,power,interim_fraction), by='Setting') %>%
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
    merge(dat_fullsample_reject_AB %>% filter(n == 2000), by=c('n','power','MonitoringType','Iteration')) %>%
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
)
dat_ABcalibrated_corrected_E1_upstraped_results <- rbind(
  #n = 40
  merge(dat_results %>% select(Setting,Iteration,P_Upstrap), 
        dat_setting %>% select(Setting,n,interim_fraction,power) , by = 'Setting') %>%
    filter(n == 40) %>%
    merge(dat_AB_setting_calibration_corrected %>% filter(interim_fraction != 1.00 & n == 40), by = c('n','interim_fraction')) %>%
    group_by(Setting,Iteration,MonitoringType) %>% 
    mutate(p_lessthan_fut = ifelse(is.na(FutilityP) == F & (as.numeric(P_Upstrap) < FutilityP), yes = 1, no = 0),
           p_lessthan_eff = ifelse(is.na(EfficacyP) == F & (as.numeric(P_Upstrap) < EfficacyP), yes = 1, no = 0)) %>%
    summarise(proportion_sim_fut = mean(p_lessthan_fut),
              proportion_sim_eff = mean(p_lessthan_eff)) %>%
    merge(dat_setting %>% select(Setting,n,interim_fraction), by = 'Setting') %>%
    merge(dat_AB_setting_calibration_corrected, by = c('n','interim_fraction','MonitoringType')) %>%
    mutate(stop_futility = ifelse(proportion_sim_fut < FutilityProp_expCorr1, yes = 1, no = 0),
           stop_efficacy = ifelse(proportion_sim_eff > EfficacyProp_expCorr1, yes = 1, no = 0)) %>%
    select(Setting,Iteration,MonitoringType,stop_futility,stop_efficacy) %>%
    merge(dat_setting %>% filter(n == 40) %>% select(n,Setting,power,interim_fraction), by='Setting') %>%
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
    merge(dat_fullsample_reject_AB %>% filter(n == 40), by=c('n','power','MonitoringType','Iteration')) %>%
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
           'RejectionRate_FE_FS' = RejectionRate_FE_FS_FE),
  #n = 160
  merge(dat_results %>% select(Setting,Iteration,P_Upstrap), 
        dat_setting %>% select(Setting,n,interim_fraction,power) , by = 'Setting') %>%
    filter(n == 160) %>%
    merge(dat_AB_setting_calibration_corrected %>% filter(interim_fraction != 1.00 & n == 160), by = c('n','interim_fraction')) %>%
    group_by(Setting,Iteration,MonitoringType) %>% 
    mutate(p_lessthan_fut = ifelse(is.na(FutilityP) == F & (as.numeric(P_Upstrap) < FutilityP), yes = 1, no = 0),
           p_lessthan_eff = ifelse(is.na(EfficacyP) == F & (as.numeric(P_Upstrap) < EfficacyP), yes = 1, no = 0)) %>%
    summarise(proportion_sim_fut = mean(p_lessthan_fut),
              proportion_sim_eff = mean(p_lessthan_eff)) %>%
    merge(dat_setting %>% select(Setting,n,interim_fraction), by = 'Setting') %>%
    merge(dat_AB_setting_calibration_corrected, by = c('n','interim_fraction','MonitoringType')) %>%
    mutate(stop_futility = ifelse(proportion_sim_fut < FutilityProp_expCorr1, yes = 1, no = 0),
           stop_efficacy = ifelse(proportion_sim_eff > EfficacyProp_expCorr1, yes = 1, no = 0)) %>%
    select(Setting,Iteration,MonitoringType,stop_futility,stop_efficacy) %>%
    merge(dat_setting %>% filter(n == 160) %>% select(n,Setting,power,interim_fraction), by='Setting') %>%
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
    merge(dat_fullsample_reject_AB %>% filter(n == 160), by=c('n','power','MonitoringType','Iteration')) %>%
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
           'RejectionRate_FE_FS' = RejectionRate_FE_FS_FE),
  #n = 600
  merge(dat_results %>% select(Setting,Iteration,P_Upstrap), 
        dat_setting %>% select(Setting,n,interim_fraction,power) , by = 'Setting') %>%
    filter(n == 600) %>%
    merge(dat_AB_setting_calibration_corrected %>% filter(interim_fraction != 1.00 & n == 600), by = c('n','interim_fraction')) %>%
    group_by(Setting,Iteration,MonitoringType) %>% 
    mutate(p_lessthan_fut = ifelse(is.na(FutilityP) == F & (as.numeric(P_Upstrap) < FutilityP), yes = 1, no = 0),
           p_lessthan_eff = ifelse(is.na(EfficacyP) == F & (as.numeric(P_Upstrap) < EfficacyP), yes = 1, no = 0)) %>%
    summarise(proportion_sim_fut = mean(p_lessthan_fut),
              proportion_sim_eff = mean(p_lessthan_eff)) %>%
    merge(dat_setting %>% select(Setting,n,interim_fraction), by = 'Setting') %>%
    merge(dat_AB_setting_calibration_corrected, by = c('n','interim_fraction','MonitoringType')) %>%
    mutate(stop_futility = ifelse(proportion_sim_fut < FutilityProp_expCorr1, yes = 1, no = 0),
           stop_efficacy = ifelse(proportion_sim_eff > EfficacyProp_expCorr1, yes = 1, no = 0)) %>%
    select(Setting,Iteration,MonitoringType,stop_futility,stop_efficacy) %>%
    merge(dat_setting %>% filter(n == 600) %>% select(n,Setting,power,interim_fraction), by='Setting') %>%
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
    merge(dat_fullsample_reject_AB %>% filter(n == 600), by=c('n','power','MonitoringType','Iteration')) %>%
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
           'RejectionRate_FE_FS' = RejectionRate_FE_FS_FE),
  #n = 2000
  merge(dat_results %>% select(Setting,Iteration,P_Upstrap), 
        dat_setting %>% select(Setting,n,interim_fraction,power) , by = 'Setting') %>%
    filter(n == 2000) %>%
    merge(dat_AB_setting_calibration_corrected %>% filter(interim_fraction != 1.00 & n == 2000), by = c('n','interim_fraction')) %>%
    group_by(Setting,Iteration,MonitoringType) %>% 
    mutate(p_lessthan_fut = ifelse(is.na(FutilityP) == F & (as.numeric(P_Upstrap) < FutilityP), yes = 1, no = 0),
           p_lessthan_eff = ifelse(is.na(EfficacyP) == F & (as.numeric(P_Upstrap) < EfficacyP), yes = 1, no = 0)) %>%
    summarise(proportion_sim_fut = mean(p_lessthan_fut),
              proportion_sim_eff = mean(p_lessthan_eff)) %>%
    merge(dat_setting %>% select(Setting,n,interim_fraction), by = 'Setting') %>%
    merge(dat_AB_setting_calibration_corrected, by = c('n','interim_fraction','MonitoringType')) %>%
    mutate(stop_futility = ifelse(proportion_sim_fut < FutilityProp_expCorr1, yes = 1, no = 0),
           stop_efficacy = ifelse(proportion_sim_eff > EfficacyProp_expCorr1, yes = 1, no = 0)) %>%
    select(Setting,Iteration,MonitoringType,stop_futility,stop_efficacy) %>%
    merge(dat_setting %>% filter(n == 2000) %>% select(n,Setting,power,interim_fraction), by='Setting') %>%
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
    merge(dat_fullsample_reject_AB %>% filter(n == 2000), by=c('n','power','MonitoringType','Iteration')) %>%
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
)
dat_ABcalibrated_corrected_E2_upstraped_results <- rbind(
  #n = 40
  merge(dat_results %>% select(Setting,Iteration,P_Upstrap), 
        dat_setting %>% select(Setting,n,interim_fraction,power) , by = 'Setting') %>%
    filter(n == 40) %>%
    merge(dat_AB_setting_calibration_corrected %>% filter(interim_fraction != 1.00 & n == 40), by = c('n','interim_fraction')) %>%
    group_by(Setting,Iteration,MonitoringType) %>% 
    mutate(p_lessthan_fut = ifelse(is.na(FutilityP) == F & (as.numeric(P_Upstrap) < FutilityP), yes = 1, no = 0),
           p_lessthan_eff = ifelse(is.na(EfficacyP) == F & (as.numeric(P_Upstrap) < EfficacyP), yes = 1, no = 0)) %>%
    summarise(proportion_sim_fut = mean(p_lessthan_fut),
              proportion_sim_eff = mean(p_lessthan_eff)) %>%
    merge(dat_setting %>% select(Setting,n,interim_fraction), by = 'Setting') %>%
    merge(dat_AB_setting_calibration_corrected, by = c('n','interim_fraction','MonitoringType')) %>%
    mutate(stop_futility = ifelse(proportion_sim_fut < FutilityProp_expCorr2, yes = 1, no = 0),
           stop_efficacy = ifelse(proportion_sim_eff > EfficacyProp_expCorr2, yes = 1, no = 0)) %>%
    select(Setting,Iteration,MonitoringType,stop_futility,stop_efficacy) %>%
    merge(dat_setting %>% filter(n == 40) %>% select(n,Setting,power,interim_fraction), by='Setting') %>%
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
    merge(dat_fullsample_reject_AB %>% filter(n == 40), by=c('n','power','MonitoringType','Iteration')) %>%
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
           'RejectionRate_FE_FS' = RejectionRate_FE_FS_FE),
  #n = 160
  merge(dat_results %>% select(Setting,Iteration,P_Upstrap), 
        dat_setting %>% select(Setting,n,interim_fraction,power) , by = 'Setting') %>%
    filter(n == 160) %>%
    merge(dat_AB_setting_calibration_corrected %>% filter(interim_fraction != 1.00 & n == 160), by = c('n','interim_fraction')) %>%
    group_by(Setting,Iteration,MonitoringType) %>% 
    mutate(p_lessthan_fut = ifelse(is.na(FutilityP) == F & (as.numeric(P_Upstrap) < FutilityP), yes = 1, no = 0),
           p_lessthan_eff = ifelse(is.na(EfficacyP) == F & (as.numeric(P_Upstrap) < EfficacyP), yes = 1, no = 0)) %>%
    summarise(proportion_sim_fut = mean(p_lessthan_fut),
              proportion_sim_eff = mean(p_lessthan_eff)) %>%
    merge(dat_setting %>% select(Setting,n,interim_fraction), by = 'Setting') %>%
    merge(dat_AB_setting_calibration_corrected, by = c('n','interim_fraction','MonitoringType')) %>%
    mutate(stop_futility = ifelse(proportion_sim_fut < FutilityProp_expCorr2, yes = 1, no = 0),
           stop_efficacy = ifelse(proportion_sim_eff > EfficacyProp_expCorr2, yes = 1, no = 0)) %>%
    select(Setting,Iteration,MonitoringType,stop_futility,stop_efficacy) %>%
    merge(dat_setting %>% filter(n == 160) %>% select(n,Setting,power,interim_fraction), by='Setting') %>%
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
    merge(dat_fullsample_reject_AB %>% filter(n == 160), by=c('n','power','MonitoringType','Iteration')) %>%
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
           'RejectionRate_FE_FS' = RejectionRate_FE_FS_FE),
  #n = 600
  merge(dat_results %>% select(Setting,Iteration,P_Upstrap), 
        dat_setting %>% select(Setting,n,interim_fraction,power) , by = 'Setting') %>%
    filter(n == 600) %>%
    merge(dat_AB_setting_calibration_corrected %>% filter(interim_fraction != 1.00 & n == 600), by = c('n','interim_fraction')) %>%
    group_by(Setting,Iteration,MonitoringType) %>% 
    mutate(p_lessthan_fut = ifelse(is.na(FutilityP) == F & (as.numeric(P_Upstrap) < FutilityP), yes = 1, no = 0),
           p_lessthan_eff = ifelse(is.na(EfficacyP) == F & (as.numeric(P_Upstrap) < EfficacyP), yes = 1, no = 0)) %>%
    summarise(proportion_sim_fut = mean(p_lessthan_fut),
              proportion_sim_eff = mean(p_lessthan_eff)) %>%
    merge(dat_setting %>% select(Setting,n,interim_fraction), by = 'Setting') %>%
    merge(dat_AB_setting_calibration_corrected, by = c('n','interim_fraction','MonitoringType')) %>%
    mutate(stop_futility = ifelse(proportion_sim_fut < FutilityProp_expCorr2, yes = 1, no = 0),
           stop_efficacy = ifelse(proportion_sim_eff > EfficacyProp_expCorr2, yes = 1, no = 0)) %>%
    select(Setting,Iteration,MonitoringType,stop_futility,stop_efficacy) %>%
    merge(dat_setting %>% filter(n == 600) %>% select(n,Setting,power,interim_fraction), by='Setting') %>%
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
    merge(dat_fullsample_reject_AB %>% filter(n == 600), by=c('n','power','MonitoringType','Iteration')) %>%
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
           'RejectionRate_FE_FS' = RejectionRate_FE_FS_FE),
  #n = 2000
  merge(dat_results %>% select(Setting,Iteration,P_Upstrap), 
        dat_setting %>% select(Setting,n,interim_fraction,power) , by = 'Setting') %>%
    filter(n == 2000) %>%
    merge(dat_AB_setting_calibration_corrected %>% filter(interim_fraction != 1.00 & n == 2000), by = c('n','interim_fraction')) %>%
    group_by(Setting,Iteration,MonitoringType) %>% 
    mutate(p_lessthan_fut = ifelse(is.na(FutilityP) == F & (as.numeric(P_Upstrap) < FutilityP), yes = 1, no = 0),
           p_lessthan_eff = ifelse(is.na(EfficacyP) == F & (as.numeric(P_Upstrap) < EfficacyP), yes = 1, no = 0)) %>%
    summarise(proportion_sim_fut = mean(p_lessthan_fut),
              proportion_sim_eff = mean(p_lessthan_eff)) %>%
    merge(dat_setting %>% select(Setting,n,interim_fraction), by = 'Setting') %>%
    merge(dat_AB_setting_calibration_corrected, by = c('n','interim_fraction','MonitoringType')) %>%
    mutate(stop_futility = ifelse(proportion_sim_fut < FutilityProp_expCorr2, yes = 1, no = 0),
           stop_efficacy = ifelse(proportion_sim_eff > EfficacyProp_expCorr2, yes = 1, no = 0)) %>%
    select(Setting,Iteration,MonitoringType,stop_futility,stop_efficacy) %>%
    merge(dat_setting %>% filter(n == 2000) %>% select(n,Setting,power,interim_fraction), by='Setting') %>%
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
    merge(dat_fullsample_reject_AB %>% filter(n == 2000), by=c('n','power','MonitoringType','Iteration')) %>%
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
)

###############################################################
### REPEAT WITHOUT 25% STOPPING POINT                       ###
###############################################################

#apply calibration corrections (multiplicative, exponential) to AU results
dat_arbitrary_corrected_M1_upstraped_results_no25 <- dat_results %>%
  group_by(Setting,Iteration) %>% 
  mutate(p_lessthan_fut = ifelse(as.numeric(P_Upstrap) < 0.05, yes = 1, no = 0),
         p_lessthan_eff = ifelse(as.numeric(P_Upstrap) < 0.05, yes = 1, no = 0)) %>%
  summarise(proportion_sim_fut = mean(p_lessthan_fut),
            proportion_sim_eff = mean(p_lessthan_eff)) %>%
  merge(dat_setting_calibration %>% select(Setting,interim_fraction,Prop_futility,Prop_efficacy), by = 'Setting') %>%
  mutate(stop_futility = ifelse(proportion_sim_fut < 0.05/(4+1) & interim_fraction != 0.25, yes = 1, no = 0),
         stop_efficacy = ifelse(proportion_sim_eff > 0.80+(1-0.80)/(4+1) & interim_fraction != 0.25, yes = 1, no = 0)) %>%
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
dat_arbitrary_corrected_M2_upstraped_results_no25 <- dat_results %>%
  group_by(Setting,Iteration) %>% 
  mutate(p_lessthan_fut = ifelse(as.numeric(P_Upstrap) < 0.05, yes = 1, no = 0),
         p_lessthan_eff = ifelse(as.numeric(P_Upstrap) < 0.05, yes = 1, no = 0)) %>%
  summarise(proportion_sim_fut = mean(p_lessthan_fut),
            proportion_sim_eff = mean(p_lessthan_eff)) %>%
  merge(dat_setting_calibration %>% select(Setting,interim_fraction,Prop_futility,Prop_efficacy), by = 'Setting') %>%
  mutate(stop_futility = ifelse(proportion_sim_fut < 0.05/(4) & interim_fraction != 0.25, yes = 1, no = 0),
         stop_efficacy = ifelse(proportion_sim_eff > 0.80+(1-0.80)/(4) & interim_fraction != 0.25, yes = 1, no = 0)) %>%
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
dat_arbitrary_corrected_M3_upstraped_results_no25 <- dat_results %>%
  group_by(Setting,Iteration) %>% 
  mutate(p_lessthan_fut = ifelse(as.numeric(P_Upstrap) < 0.05, yes = 1, no = 0),
         p_lessthan_eff = ifelse(as.numeric(P_Upstrap) < 0.05, yes = 1, no = 0)) %>%
  summarise(proportion_sim_fut = mean(p_lessthan_fut),
            proportion_sim_eff = mean(p_lessthan_eff)) %>%
  merge(dat_setting_calibration %>% select(Setting,interim_fraction,Prop_futility,Prop_efficacy), by = 'Setting') %>%
  mutate(stop_futility = ifelse(proportion_sim_fut < 0.05/(4-1) & interim_fraction != 0.25, yes = 1, no = 0),
         stop_efficacy = ifelse(proportion_sim_eff > 0.80+(1-0.80)/(4-1) & interim_fraction != 0.25, yes = 1, no = 0)) %>%
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
dat_arbitrary_corrected_M4_upstraped_results_no25 <- dat_results %>%
  group_by(Setting,Iteration) %>% 
  mutate(p_lessthan_fut = ifelse(as.numeric(P_Upstrap) < 0.05, yes = 1, no = 0),
         p_lessthan_eff = ifelse(as.numeric(P_Upstrap) < 0.05, yes = 1, no = 0)) %>%
  summarise(proportion_sim_fut = mean(p_lessthan_fut),
            proportion_sim_eff = mean(p_lessthan_eff)) %>%
  merge(dat_setting_calibration %>% select(Setting,interim_fraction,Prop_futility,Prop_efficacy), by = 'Setting') %>%
  mutate(stop_futility = ifelse(proportion_sim_fut < 0.05/(4-1.5) & interim_fraction != 0.25, yes = 1, no = 0),
         stop_efficacy = ifelse(proportion_sim_eff > 0.80+(1-0.80) & interim_fraction != 0.25/(4-1.5), yes = 1, no = 0)) %>%
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
dat_arbitrary_corrected_M5_upstraped_results_no25 <- dat_results %>%
  group_by(Setting,Iteration) %>% 
  mutate(p_lessthan_fut = ifelse(as.numeric(P_Upstrap) < 0.05, yes = 1, no = 0),
         p_lessthan_eff = ifelse(as.numeric(P_Upstrap) < 0.05, yes = 1, no = 0)) %>%
  summarise(proportion_sim_fut = mean(p_lessthan_fut),
            proportion_sim_eff = mean(p_lessthan_eff)) %>%
  merge(dat_setting_calibration %>% select(Setting,interim_fraction,Prop_futility,Prop_efficacy), by = 'Setting') %>%
  mutate(stop_futility = ifelse(proportion_sim_fut < 0.05/(4-2) & interim_fraction != 0.25, yes = 1, no = 0),
         stop_efficacy = ifelse(proportion_sim_eff > 0.80+(1-0.80)/(4-2) & interim_fraction != 0.25, yes = 1, no = 0)) %>%
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
dat_arbitrary_corrected_E1_upstraped_results_no25 <- dat_results %>%
  group_by(Setting,Iteration) %>% 
  mutate(p_lessthan_fut = ifelse(as.numeric(P_Upstrap) < 0.05, yes = 1, no = 0),
         p_lessthan_eff = ifelse(as.numeric(P_Upstrap) < 0.05, yes = 1, no = 0)) %>%
  summarise(proportion_sim_fut = mean(p_lessthan_fut),
            proportion_sim_eff = mean(p_lessthan_eff)) %>%
  merge(dat_setting_calibration %>% select(Setting,interim_fraction,Prop_futility,Prop_efficacy), by = 'Setting') %>%
  mutate(stop_futility = ifelse(proportion_sim_fut < 0.05^0.95 & interim_fraction != 0.25, yes = 1, no = 0),
         stop_efficacy = ifelse(proportion_sim_eff > 0.80^(1-0.95) & interim_fraction != 0.25, yes = 1, no = 0)) %>%
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
dat_arbitrary_corrected_E2_upstraped_results_no25 <- dat_results %>%
  group_by(Setting,Iteration) %>% 
  mutate(p_lessthan_fut = ifelse(as.numeric(P_Upstrap) < 0.05, yes = 1, no = 0),
         p_lessthan_eff = ifelse(as.numeric(P_Upstrap) < 0.05, yes = 1, no = 0)) %>%
  summarise(proportion_sim_fut = mean(p_lessthan_fut),
            proportion_sim_eff = mean(p_lessthan_eff)) %>%
  merge(dat_setting_calibration %>% select(Setting,interim_fraction,Prop_futility,Prop_efficacy), by = 'Setting') %>%
  mutate(stop_futility = ifelse(proportion_sim_fut < 0.05^0.70 & interim_fraction != 0.25, yes = 1, no = 0),
         stop_efficacy = ifelse(proportion_sim_eff > 0.80^(1-0.70) & interim_fraction != 0.25, yes = 1, no = 0)) %>%
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

#apply (corrected) calibrated upstrapping stopping boundaries to simulated results
dat_calibrated_corrected_M1_upstraped_results_no25 <- merge(dat_results, dat_setting_calibration_corrected %>% 
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
          select(Setting,interim_fraction,P_futility,Prop_efficacy_multCorr1,P_efficacy,Prop_futility_multCorr1) %>%
          rename(Prop_efficacy = Prop_efficacy_multCorr1,
                 Prop_futility = Prop_futility_multCorr1), 
        by = 'Setting') %>%
  mutate(stop_futility = ifelse(proportion_sim_fut < Prop_futility  & interim_fraction != 0.25, yes = 1, no = 0),
         stop_efficacy = ifelse(proportion_sim_eff > Prop_efficacy  & interim_fraction != 0.25, yes = 1, no = 0)) %>%
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
dat_calibrated_corrected_M2_upstraped_results_no25 <- merge(dat_results, dat_setting_calibration_corrected %>% 
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
          select(Setting,interim_fraction,P_futility,Prop_efficacy_multCorr2,P_efficacy,Prop_futility_multCorr2) %>%
          rename(Prop_efficacy = Prop_efficacy_multCorr2,
                 Prop_futility = Prop_futility_multCorr2), 
        by = 'Setting') %>%
  mutate(stop_futility = ifelse(proportion_sim_fut < Prop_futility & interim_fraction != 0.25, yes = 1, no = 0),
         stop_efficacy = ifelse(proportion_sim_eff > Prop_efficacy & interim_fraction != 0.25, yes = 1, no = 0)) %>%
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
dat_calibrated_corrected_M3_upstraped_results_no25 <- merge(dat_results, dat_setting_calibration_corrected %>% 
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
          select(Setting,interim_fraction,P_futility,Prop_efficacy_multCorr3,P_efficacy,Prop_futility_multCorr3) %>%
          rename(Prop_efficacy = Prop_efficacy_multCorr3,
                 Prop_futility = Prop_futility_multCorr3), 
        by = 'Setting') %>%
  mutate(stop_futility = ifelse(proportion_sim_fut < Prop_futility & interim_fraction != 0.25, yes = 1, no = 0),
         stop_efficacy = ifelse(proportion_sim_eff > Prop_efficacy & interim_fraction != 0.25, yes = 1, no = 0)) %>%
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
dat_calibrated_corrected_M4_upstraped_results_no25 <- merge(dat_results, dat_setting_calibration_corrected %>% 
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
          select(Setting,interim_fraction,P_futility,Prop_efficacy_multCorr4,P_efficacy,Prop_futility_multCorr4) %>%
          rename(Prop_efficacy = Prop_efficacy_multCorr4,
                 Prop_futility = Prop_futility_multCorr4), 
        by = 'Setting') %>%
  mutate(stop_futility = ifelse(proportion_sim_fut < Prop_futility & interim_fraction != 0.25, yes = 1, no = 0),
         stop_efficacy = ifelse(proportion_sim_eff > Prop_efficacy & interim_fraction != 0.25, yes = 1, no = 0)) %>%
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
dat_calibrated_corrected_M5_upstraped_results_no25 <- merge(dat_results, dat_setting_calibration_corrected %>% 
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
          select(Setting,interim_fraction,P_futility,Prop_efficacy_multCorr5,P_efficacy,Prop_futility_multCorr5) %>%
          rename(Prop_efficacy = Prop_efficacy_multCorr5,
                 Prop_futility = Prop_futility_multCorr5), 
        by = 'Setting') %>%
  mutate(stop_futility = ifelse(proportion_sim_fut < Prop_futility & interim_fraction != 0.25, yes = 1, no = 0),
         stop_efficacy = ifelse(proportion_sim_eff > Prop_efficacy & interim_fraction != 0.25, yes = 1, no = 0)) %>%
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
dat_calibrated_corrected_E1_upstraped_results_no25 <- merge(dat_results, dat_setting_calibration_corrected %>% 
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
          select(Setting,interim_fraction,P_futility,Prop_efficacy_expCorr1,P_efficacy,Prop_futility_expCorr1) %>%
          rename(Prop_efficacy = Prop_efficacy_expCorr1,
                 Prop_futility = Prop_futility_expCorr1), 
        by = 'Setting') %>%
  mutate(stop_futility = ifelse(proportion_sim_fut < Prop_futility & interim_fraction != 0.25, yes = 1, no = 0),
         stop_efficacy = ifelse(proportion_sim_eff > Prop_efficacy & interim_fraction != 0.25, yes = 1, no = 0)) %>%
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
dat_calibrated_corrected_E2_upstraped_results_no25 <- merge(dat_results, dat_setting_calibration_corrected %>% 
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
          select(Setting,interim_fraction,P_futility,Prop_efficacy_expCorr2,P_efficacy,Prop_futility_expCorr2) %>%
          rename(Prop_efficacy = Prop_efficacy_expCorr2,
                 Prop_futility = Prop_futility_expCorr2), 
        by = 'Setting') %>%
  mutate(stop_futility = ifelse(proportion_sim_fut < Prop_futility & interim_fraction != 0.25, yes = 1, no = 0),
         stop_efficacy = ifelse(proportion_sim_eff > Prop_efficacy & interim_fraction != 0.25, yes = 1, no = 0)) %>%
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

#apply calibration corrections (multiplicative, exponential) to GU results
dat_ABcalibrated_corrected_M1_upstraped_results_no25 <- rbind(
  #n = 40
  merge(dat_results %>% select(Setting,Iteration,P_Upstrap), 
        dat_setting %>% select(Setting,n,interim_fraction,power) , by = 'Setting') %>%
    filter(n == 40) %>%
    merge(dat_AB_setting_calibration_corrected %>% filter(interim_fraction != 1.00 & n == 40), by = c('n','interim_fraction')) %>%
    group_by(Setting,Iteration,MonitoringType) %>% 
    mutate(p_lessthan_fut = ifelse(is.na(FutilityP) == F & (as.numeric(P_Upstrap) < FutilityP), yes = 1, no = 0),
           p_lessthan_eff = ifelse(is.na(EfficacyP) == F & (as.numeric(P_Upstrap) < EfficacyP), yes = 1, no = 0)) %>%
    summarise(proportion_sim_fut = mean(p_lessthan_fut),
              proportion_sim_eff = mean(p_lessthan_eff)) %>%
    merge(dat_setting %>% select(Setting,n,interim_fraction), by = 'Setting') %>%
    merge(dat_AB_setting_calibration_corrected, by = c('n','interim_fraction','MonitoringType')) %>%
    mutate(stop_futility = ifelse(proportion_sim_fut < FutilityProp_multCorr1 & interim_fraction != 0.25, yes = 1, no = 0),
           stop_efficacy = ifelse(proportion_sim_eff > EfficacyProp_multCorr1 & interim_fraction != 0.25, yes = 1, no = 0)) %>%
    select(Setting,Iteration,MonitoringType,stop_futility,stop_efficacy) %>%
    merge(dat_setting %>% filter(n == 40) %>% select(n,Setting,power,interim_fraction), by='Setting') %>%
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
    merge(dat_fullsample_reject_AB %>% filter(n == 40), by=c('n','power','MonitoringType','Iteration')) %>%
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
           'RejectionRate_FE_FS' = RejectionRate_FE_FS_FE),
  #n = 160
  merge(dat_results %>% select(Setting,Iteration,P_Upstrap), 
        dat_setting %>% select(Setting,n,interim_fraction,power) , by = 'Setting') %>%
    filter(n == 160) %>%
    merge(dat_AB_setting_calibration_corrected %>% filter(interim_fraction != 1.00 & n == 160), by = c('n','interim_fraction')) %>%
    group_by(Setting,Iteration,MonitoringType) %>% 
    mutate(p_lessthan_fut = ifelse(is.na(FutilityP) == F & (as.numeric(P_Upstrap) < FutilityP), yes = 1, no = 0),
           p_lessthan_eff = ifelse(is.na(EfficacyP) == F & (as.numeric(P_Upstrap) < EfficacyP), yes = 1, no = 0)) %>%
    summarise(proportion_sim_fut = mean(p_lessthan_fut),
              proportion_sim_eff = mean(p_lessthan_eff)) %>%
    merge(dat_setting %>% select(Setting,n,interim_fraction), by = 'Setting') %>%
    merge(dat_AB_setting_calibration_corrected, by = c('n','interim_fraction','MonitoringType')) %>%
    mutate(stop_futility = ifelse(proportion_sim_fut < FutilityProp_multCorr1 & interim_fraction != 0.25, yes = 1, no = 0),
           stop_efficacy = ifelse(proportion_sim_eff > EfficacyProp_multCorr1 & interim_fraction != 0.25, yes = 1, no = 0)) %>%
    select(Setting,Iteration,MonitoringType,stop_futility,stop_efficacy) %>%
    merge(dat_setting %>% filter(n == 160) %>% select(n,Setting,power,interim_fraction), by='Setting') %>%
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
    merge(dat_fullsample_reject_AB %>% filter(n == 160), by=c('n','power','MonitoringType','Iteration')) %>%
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
           'RejectionRate_FE_FS' = RejectionRate_FE_FS_FE),
  #n = 600
  merge(dat_results %>% select(Setting,Iteration,P_Upstrap), 
        dat_setting %>% select(Setting,n,interim_fraction,power) , by = 'Setting') %>%
    filter(n == 600) %>%
    merge(dat_AB_setting_calibration_corrected %>% filter(interim_fraction != 1.00 & n == 600), by = c('n','interim_fraction')) %>%
    group_by(Setting,Iteration,MonitoringType) %>% 
    mutate(p_lessthan_fut = ifelse(is.na(FutilityP) == F & (as.numeric(P_Upstrap) < FutilityP), yes = 1, no = 0),
           p_lessthan_eff = ifelse(is.na(EfficacyP) == F & (as.numeric(P_Upstrap) < EfficacyP), yes = 1, no = 0)) %>%
    summarise(proportion_sim_fut = mean(p_lessthan_fut),
              proportion_sim_eff = mean(p_lessthan_eff)) %>%
    merge(dat_setting %>% select(Setting,n,interim_fraction), by = 'Setting') %>%
    merge(dat_AB_setting_calibration_corrected, by = c('n','interim_fraction','MonitoringType')) %>%
    mutate(stop_futility = ifelse(proportion_sim_fut < FutilityProp_multCorr1 & interim_fraction != 0.25, yes = 1, no = 0),
           stop_efficacy = ifelse(proportion_sim_eff > EfficacyProp_multCorr1 & interim_fraction != 0.25, yes = 1, no = 0)) %>%
    select(Setting,Iteration,MonitoringType,stop_futility,stop_efficacy) %>%
    merge(dat_setting %>% filter(n == 600) %>% select(n,Setting,power,interim_fraction), by='Setting') %>%
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
    merge(dat_fullsample_reject_AB %>% filter(n == 600), by=c('n','power','MonitoringType','Iteration')) %>%
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
           'RejectionRate_FE_FS' = RejectionRate_FE_FS_FE),
  #n = 2000
  merge(dat_results %>% select(Setting,Iteration,P_Upstrap), 
        dat_setting %>% select(Setting,n,interim_fraction,power) , by = 'Setting') %>%
    filter(n == 2000) %>%
    merge(dat_AB_setting_calibration_corrected %>% filter(interim_fraction != 1.00 & n == 2000), by = c('n','interim_fraction')) %>%
    group_by(Setting,Iteration,MonitoringType) %>% 
    mutate(p_lessthan_fut = ifelse(is.na(FutilityP) == F & (as.numeric(P_Upstrap) < FutilityP), yes = 1, no = 0),
           p_lessthan_eff = ifelse(is.na(EfficacyP) == F & (as.numeric(P_Upstrap) < EfficacyP), yes = 1, no = 0)) %>%
    summarise(proportion_sim_fut = mean(p_lessthan_fut),
              proportion_sim_eff = mean(p_lessthan_eff)) %>%
    merge(dat_setting %>% select(Setting,n,interim_fraction), by = 'Setting') %>%
    merge(dat_AB_setting_calibration_corrected, by = c('n','interim_fraction','MonitoringType')) %>%
    mutate(stop_futility = ifelse(proportion_sim_fut < FutilityProp_multCorr1 & interim_fraction != 0.25, yes = 1, no = 0),
           stop_efficacy = ifelse(proportion_sim_eff > EfficacyProp_multCorr1 & interim_fraction != 0.25, yes = 1, no = 0)) %>%
    select(Setting,Iteration,MonitoringType,stop_futility,stop_efficacy) %>%
    merge(dat_setting %>% filter(n == 2000) %>% select(n,Setting,power,interim_fraction), by='Setting') %>%
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
    merge(dat_fullsample_reject_AB %>% filter(n == 2000), by=c('n','power','MonitoringType','Iteration')) %>%
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
)
dat_ABcalibrated_corrected_M2_upstraped_results_no25 <- rbind(
  #n = 40
  merge(dat_results %>% select(Setting,Iteration,P_Upstrap), 
        dat_setting %>% select(Setting,n,interim_fraction,power) , by = 'Setting') %>%
    filter(n == 40) %>%
    merge(dat_AB_setting_calibration_corrected %>% filter(interim_fraction != 1.00 & n == 40), by = c('n','interim_fraction')) %>%
    group_by(Setting,Iteration,MonitoringType) %>% 
    mutate(p_lessthan_fut = ifelse(is.na(FutilityP) == F & (as.numeric(P_Upstrap) < FutilityP), yes = 1, no = 0),
           p_lessthan_eff = ifelse(is.na(EfficacyP) == F & (as.numeric(P_Upstrap) < EfficacyP), yes = 1, no = 0)) %>%
    summarise(proportion_sim_fut = mean(p_lessthan_fut),
              proportion_sim_eff = mean(p_lessthan_eff)) %>%
    merge(dat_setting %>% select(Setting,n,interim_fraction), by = 'Setting') %>%
    merge(dat_AB_setting_calibration_corrected, by = c('n','interim_fraction','MonitoringType')) %>%
    mutate(stop_futility = ifelse(proportion_sim_fut < FutilityProp_multCorr2 & interim_fraction != 0.25, yes = 1, no = 0),
           stop_efficacy = ifelse(proportion_sim_eff > EfficacyProp_multCorr2 & interim_fraction != 0.25, yes = 1, no = 0)) %>%
    select(Setting,Iteration,MonitoringType,stop_futility,stop_efficacy) %>%
    merge(dat_setting %>% filter(n == 40) %>% select(n,Setting,power,interim_fraction), by='Setting') %>%
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
    merge(dat_fullsample_reject_AB %>% filter(n == 40), by=c('n','power','MonitoringType','Iteration')) %>%
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
           'RejectionRate_FE_FS' = RejectionRate_FE_FS_FE),
  #n = 160
  merge(dat_results %>% select(Setting,Iteration,P_Upstrap), 
        dat_setting %>% select(Setting,n,interim_fraction,power) , by = 'Setting') %>%
    filter(n == 160) %>%
    merge(dat_AB_setting_calibration_corrected %>% filter(interim_fraction != 1.00 & n == 160), by = c('n','interim_fraction')) %>%
    group_by(Setting,Iteration,MonitoringType) %>% 
    mutate(p_lessthan_fut = ifelse(is.na(FutilityP) == F & (as.numeric(P_Upstrap) < FutilityP), yes = 1, no = 0),
           p_lessthan_eff = ifelse(is.na(EfficacyP) == F & (as.numeric(P_Upstrap) < EfficacyP), yes = 1, no = 0)) %>%
    summarise(proportion_sim_fut = mean(p_lessthan_fut),
              proportion_sim_eff = mean(p_lessthan_eff)) %>%
    merge(dat_setting %>% select(Setting,n,interim_fraction), by = 'Setting') %>%
    merge(dat_AB_setting_calibration_corrected, by = c('n','interim_fraction','MonitoringType')) %>%
    mutate(stop_futility = ifelse(proportion_sim_fut < FutilityProp_multCorr2 & interim_fraction != 0.25, yes = 1, no = 0),
           stop_efficacy = ifelse(proportion_sim_eff > EfficacyProp_multCorr2 & interim_fraction != 0.25, yes = 1, no = 0)) %>%
    select(Setting,Iteration,MonitoringType,stop_futility,stop_efficacy) %>%
    merge(dat_setting %>% filter(n == 160) %>% select(n,Setting,power,interim_fraction), by='Setting') %>%
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
    merge(dat_fullsample_reject_AB %>% filter(n == 160), by=c('n','power','MonitoringType','Iteration')) %>%
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
           'RejectionRate_FE_FS' = RejectionRate_FE_FS_FE),
  #n = 600
  merge(dat_results %>% select(Setting,Iteration,P_Upstrap), 
        dat_setting %>% select(Setting,n,interim_fraction,power) , by = 'Setting') %>%
    filter(n == 600) %>%
    merge(dat_AB_setting_calibration_corrected %>% filter(interim_fraction != 1.00 & n == 600), by = c('n','interim_fraction')) %>%
    group_by(Setting,Iteration,MonitoringType) %>% 
    mutate(p_lessthan_fut = ifelse(is.na(FutilityP) == F & (as.numeric(P_Upstrap) < FutilityP), yes = 1, no = 0),
           p_lessthan_eff = ifelse(is.na(EfficacyP) == F & (as.numeric(P_Upstrap) < EfficacyP), yes = 1, no = 0)) %>%
    summarise(proportion_sim_fut = mean(p_lessthan_fut),
              proportion_sim_eff = mean(p_lessthan_eff)) %>%
    merge(dat_setting %>% select(Setting,n,interim_fraction), by = 'Setting') %>%
    merge(dat_AB_setting_calibration_corrected, by = c('n','interim_fraction','MonitoringType')) %>%
    mutate(stop_futility = ifelse(proportion_sim_fut < FutilityProp_multCorr2 & interim_fraction != 0.25, yes = 1, no = 0),
           stop_efficacy = ifelse(proportion_sim_eff > EfficacyProp_multCorr2 & interim_fraction != 0.25, yes = 1, no = 0)) %>%
    select(Setting,Iteration,MonitoringType,stop_futility,stop_efficacy) %>%
    merge(dat_setting %>% filter(n == 600) %>% select(n,Setting,power,interim_fraction), by='Setting') %>%
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
    merge(dat_fullsample_reject_AB %>% filter(n == 600), by=c('n','power','MonitoringType','Iteration')) %>%
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
           'RejectionRate_FE_FS' = RejectionRate_FE_FS_FE),
  #n = 2000
  merge(dat_results %>% select(Setting,Iteration,P_Upstrap), 
        dat_setting %>% select(Setting,n,interim_fraction,power) , by = 'Setting') %>%
    filter(n == 2000) %>%
    merge(dat_AB_setting_calibration_corrected %>% filter(interim_fraction != 1.00 & n == 2000), by = c('n','interim_fraction')) %>%
    group_by(Setting,Iteration,MonitoringType) %>% 
    mutate(p_lessthan_fut = ifelse(is.na(FutilityP) == F & (as.numeric(P_Upstrap) < FutilityP), yes = 1, no = 0),
           p_lessthan_eff = ifelse(is.na(EfficacyP) == F & (as.numeric(P_Upstrap) < EfficacyP), yes = 1, no = 0)) %>%
    summarise(proportion_sim_fut = mean(p_lessthan_fut),
              proportion_sim_eff = mean(p_lessthan_eff)) %>%
    merge(dat_setting %>% select(Setting,n,interim_fraction), by = 'Setting') %>%
    merge(dat_AB_setting_calibration_corrected, by = c('n','interim_fraction','MonitoringType')) %>%
    mutate(stop_futility = ifelse(proportion_sim_fut < FutilityProp_multCorr2 & interim_fraction != 0.25, yes = 1, no = 0),
           stop_efficacy = ifelse(proportion_sim_eff > EfficacyProp_multCorr2 & interim_fraction != 0.25, yes = 1, no = 0)) %>%
    select(Setting,Iteration,MonitoringType,stop_futility,stop_efficacy) %>%
    merge(dat_setting %>% filter(n == 2000) %>% select(n,Setting,power,interim_fraction), by='Setting') %>%
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
    merge(dat_fullsample_reject_AB %>% filter(n == 2000), by=c('n','power','MonitoringType','Iteration')) %>%
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
)
dat_ABcalibrated_corrected_M3_upstraped_results_no25 <- rbind(
  #n = 40
  merge(dat_results %>% select(Setting,Iteration,P_Upstrap), 
        dat_setting %>% select(Setting,n,interim_fraction,power) , by = 'Setting') %>%
    filter(n == 40) %>%
    merge(dat_AB_setting_calibration_corrected %>% filter(interim_fraction != 1.00 & n == 40), by = c('n','interim_fraction')) %>%
    group_by(Setting,Iteration,MonitoringType) %>% 
    mutate(p_lessthan_fut = ifelse(is.na(FutilityP) == F & (as.numeric(P_Upstrap) < FutilityP), yes = 1, no = 0),
           p_lessthan_eff = ifelse(is.na(EfficacyP) == F & (as.numeric(P_Upstrap) < EfficacyP), yes = 1, no = 0)) %>%
    summarise(proportion_sim_fut = mean(p_lessthan_fut),
              proportion_sim_eff = mean(p_lessthan_eff)) %>%
    merge(dat_setting %>% select(Setting,n,interim_fraction), by = 'Setting') %>%
    merge(dat_AB_setting_calibration_corrected, by = c('n','interim_fraction','MonitoringType')) %>%
    mutate(stop_futility = ifelse(proportion_sim_fut < FutilityProp_multCorr3 & interim_fraction != 0.25, yes = 1, no = 0),
           stop_efficacy = ifelse(proportion_sim_eff > EfficacyProp_multCorr3 & interim_fraction != 0.25, yes = 1, no = 0)) %>%
    select(Setting,Iteration,MonitoringType,stop_futility,stop_efficacy) %>%
    merge(dat_setting %>% filter(n == 40) %>% select(n,Setting,power,interim_fraction), by='Setting') %>%
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
    merge(dat_fullsample_reject_AB %>% filter(n == 40), by=c('n','power','MonitoringType','Iteration')) %>%
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
           'RejectionRate_FE_FS' = RejectionRate_FE_FS_FE),
  #n = 160
  merge(dat_results %>% select(Setting,Iteration,P_Upstrap), 
        dat_setting %>% select(Setting,n,interim_fraction,power) , by = 'Setting') %>%
    filter(n == 160) %>%
    merge(dat_AB_setting_calibration_corrected %>% filter(interim_fraction != 1.00 & n == 160), by = c('n','interim_fraction')) %>%
    group_by(Setting,Iteration,MonitoringType) %>% 
    mutate(p_lessthan_fut = ifelse(is.na(FutilityP) == F & (as.numeric(P_Upstrap) < FutilityP), yes = 1, no = 0),
           p_lessthan_eff = ifelse(is.na(EfficacyP) == F & (as.numeric(P_Upstrap) < EfficacyP), yes = 1, no = 0)) %>%
    summarise(proportion_sim_fut = mean(p_lessthan_fut),
              proportion_sim_eff = mean(p_lessthan_eff)) %>%
    merge(dat_setting %>% select(Setting,n,interim_fraction), by = 'Setting') %>%
    merge(dat_AB_setting_calibration_corrected, by = c('n','interim_fraction','MonitoringType')) %>%
    mutate(stop_futility = ifelse(proportion_sim_fut < FutilityProp_multCorr3 & interim_fraction != 0.25, yes = 1, no = 0),
           stop_efficacy = ifelse(proportion_sim_eff > EfficacyProp_multCorr3 & interim_fraction != 0.25, yes = 1, no = 0)) %>%
    select(Setting,Iteration,MonitoringType,stop_futility,stop_efficacy) %>%
    merge(dat_setting %>% filter(n == 160) %>% select(n,Setting,power,interim_fraction), by='Setting') %>%
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
    merge(dat_fullsample_reject_AB %>% filter(n == 160), by=c('n','power','MonitoringType','Iteration')) %>%
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
           'RejectionRate_FE_FS' = RejectionRate_FE_FS_FE),
  #n = 600
  merge(dat_results %>% select(Setting,Iteration,P_Upstrap), 
        dat_setting %>% select(Setting,n,interim_fraction,power) , by = 'Setting') %>%
    filter(n == 600) %>%
    merge(dat_AB_setting_calibration_corrected %>% filter(interim_fraction != 1.00 & n == 600), by = c('n','interim_fraction')) %>%
    group_by(Setting,Iteration,MonitoringType) %>% 
    mutate(p_lessthan_fut = ifelse(is.na(FutilityP) == F & (as.numeric(P_Upstrap) < FutilityP), yes = 1, no = 0),
           p_lessthan_eff = ifelse(is.na(EfficacyP) == F & (as.numeric(P_Upstrap) < EfficacyP), yes = 1, no = 0)) %>%
    summarise(proportion_sim_fut = mean(p_lessthan_fut),
              proportion_sim_eff = mean(p_lessthan_eff)) %>%
    merge(dat_setting %>% select(Setting,n,interim_fraction), by = 'Setting') %>%
    merge(dat_AB_setting_calibration_corrected, by = c('n','interim_fraction','MonitoringType')) %>%
    mutate(stop_futility = ifelse(proportion_sim_fut < FutilityProp_multCorr3 & interim_fraction != 0.25, yes = 1, no = 0),
           stop_efficacy = ifelse(proportion_sim_eff > EfficacyProp_multCorr3 & interim_fraction != 0.25, yes = 1, no = 0)) %>%
    select(Setting,Iteration,MonitoringType,stop_futility,stop_efficacy) %>%
    merge(dat_setting %>% filter(n == 600) %>% select(n,Setting,power,interim_fraction), by='Setting') %>%
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
    merge(dat_fullsample_reject_AB %>% filter(n == 600), by=c('n','power','MonitoringType','Iteration')) %>%
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
           'RejectionRate_FE_FS' = RejectionRate_FE_FS_FE),
  #n = 2000
  merge(dat_results %>% select(Setting,Iteration,P_Upstrap), 
        dat_setting %>% select(Setting,n,interim_fraction,power) , by = 'Setting') %>%
    filter(n == 2000) %>%
    merge(dat_AB_setting_calibration_corrected %>% filter(interim_fraction != 1.00 & n == 2000), by = c('n','interim_fraction')) %>%
    group_by(Setting,Iteration,MonitoringType) %>% 
    mutate(p_lessthan_fut = ifelse(is.na(FutilityP) == F & (as.numeric(P_Upstrap) < FutilityP), yes = 1, no = 0),
           p_lessthan_eff = ifelse(is.na(EfficacyP) == F & (as.numeric(P_Upstrap) < EfficacyP), yes = 1, no = 0)) %>%
    summarise(proportion_sim_fut = mean(p_lessthan_fut),
              proportion_sim_eff = mean(p_lessthan_eff)) %>%
    merge(dat_setting %>% select(Setting,n,interim_fraction), by = 'Setting') %>%
    merge(dat_AB_setting_calibration_corrected, by = c('n','interim_fraction','MonitoringType')) %>%
    mutate(stop_futility = ifelse(proportion_sim_fut < FutilityProp_multCorr3 & interim_fraction != 0.25, yes = 1, no = 0),
           stop_efficacy = ifelse(proportion_sim_eff > EfficacyProp_multCorr3 & interim_fraction != 0.25, yes = 1, no = 0)) %>%
    select(Setting,Iteration,MonitoringType,stop_futility,stop_efficacy) %>%
    merge(dat_setting %>% filter(n == 2000) %>% select(n,Setting,power,interim_fraction), by='Setting') %>%
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
    merge(dat_fullsample_reject_AB %>% filter(n == 2000), by=c('n','power','MonitoringType','Iteration')) %>%
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
)
dat_ABcalibrated_corrected_M4_upstraped_results_no25 <- rbind(
  #n = 40
  merge(dat_results %>% select(Setting,Iteration,P_Upstrap), 
        dat_setting %>% select(Setting,n,interim_fraction,power) , by = 'Setting') %>%
    filter(n == 40) %>%
    merge(dat_AB_setting_calibration_corrected %>% filter(interim_fraction != 1.00 & n == 40), by = c('n','interim_fraction')) %>%
    group_by(Setting,Iteration,MonitoringType) %>% 
    mutate(p_lessthan_fut = ifelse(is.na(FutilityP) == F & (as.numeric(P_Upstrap) < FutilityP), yes = 1, no = 0),
           p_lessthan_eff = ifelse(is.na(EfficacyP) == F & (as.numeric(P_Upstrap) < EfficacyP), yes = 1, no = 0)) %>%
    summarise(proportion_sim_fut = mean(p_lessthan_fut),
              proportion_sim_eff = mean(p_lessthan_eff)) %>%
    merge(dat_setting %>% select(Setting,n,interim_fraction), by = 'Setting') %>%
    merge(dat_AB_setting_calibration_corrected, by = c('n','interim_fraction','MonitoringType')) %>%
    mutate(stop_futility = ifelse(proportion_sim_fut < FutilityProp_multCorr4 & interim_fraction != 0.25, yes = 1, no = 0),
           stop_efficacy = ifelse(proportion_sim_eff > EfficacyProp_multCorr4 & interim_fraction != 0.25, yes = 1, no = 0)) %>%
    select(Setting,Iteration,MonitoringType,stop_futility,stop_efficacy) %>%
    merge(dat_setting %>% filter(n == 40) %>% select(n,Setting,power,interim_fraction), by='Setting') %>%
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
    merge(dat_fullsample_reject_AB %>% filter(n == 40), by=c('n','power','MonitoringType','Iteration')) %>%
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
           'RejectionRate_FE_FS' = RejectionRate_FE_FS_FE),
  #n = 160
  merge(dat_results %>% select(Setting,Iteration,P_Upstrap), 
        dat_setting %>% select(Setting,n,interim_fraction,power) , by = 'Setting') %>%
    filter(n == 160) %>%
    merge(dat_AB_setting_calibration_corrected %>% filter(interim_fraction != 1.00 & n == 160), by = c('n','interim_fraction')) %>%
    group_by(Setting,Iteration,MonitoringType) %>% 
    mutate(p_lessthan_fut = ifelse(is.na(FutilityP) == F & (as.numeric(P_Upstrap) < FutilityP), yes = 1, no = 0),
           p_lessthan_eff = ifelse(is.na(EfficacyP) == F & (as.numeric(P_Upstrap) < EfficacyP), yes = 1, no = 0)) %>%
    summarise(proportion_sim_fut = mean(p_lessthan_fut),
              proportion_sim_eff = mean(p_lessthan_eff)) %>%
    merge(dat_setting %>% select(Setting,n,interim_fraction), by = 'Setting') %>%
    merge(dat_AB_setting_calibration_corrected, by = c('n','interim_fraction','MonitoringType')) %>%
    mutate(stop_futility = ifelse(proportion_sim_fut < FutilityProp_multCorr4 & interim_fraction != 0.25, yes = 1, no = 0),
           stop_efficacy = ifelse(proportion_sim_eff > EfficacyProp_multCorr4 & interim_fraction != 0.25, yes = 1, no = 0)) %>%
    select(Setting,Iteration,MonitoringType,stop_futility,stop_efficacy) %>%
    merge(dat_setting %>% filter(n == 160) %>% select(n,Setting,power,interim_fraction), by='Setting') %>%
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
    merge(dat_fullsample_reject_AB %>% filter(n == 160), by=c('n','power','MonitoringType','Iteration')) %>%
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
           'RejectionRate_FE_FS' = RejectionRate_FE_FS_FE),
  #n = 600
  merge(dat_results %>% select(Setting,Iteration,P_Upstrap), 
        dat_setting %>% select(Setting,n,interim_fraction,power) , by = 'Setting') %>%
    filter(n == 600) %>%
    merge(dat_AB_setting_calibration_corrected %>% filter(interim_fraction != 1.00 & n == 600), by = c('n','interim_fraction')) %>%
    group_by(Setting,Iteration,MonitoringType) %>% 
    mutate(p_lessthan_fut = ifelse(is.na(FutilityP) == F & (as.numeric(P_Upstrap) < FutilityP), yes = 1, no = 0),
           p_lessthan_eff = ifelse(is.na(EfficacyP) == F & (as.numeric(P_Upstrap) < EfficacyP), yes = 1, no = 0)) %>%
    summarise(proportion_sim_fut = mean(p_lessthan_fut),
              proportion_sim_eff = mean(p_lessthan_eff)) %>%
    merge(dat_setting %>% select(Setting,n,interim_fraction), by = 'Setting') %>%
    merge(dat_AB_setting_calibration_corrected, by = c('n','interim_fraction','MonitoringType')) %>%
    mutate(stop_futility = ifelse(proportion_sim_fut < FutilityProp_multCorr4 & interim_fraction != 0.25, yes = 1, no = 0),
           stop_efficacy = ifelse(proportion_sim_eff > EfficacyProp_multCorr4 & interim_fraction != 0.25, yes = 1, no = 0)) %>%
    select(Setting,Iteration,MonitoringType,stop_futility,stop_efficacy) %>%
    merge(dat_setting %>% filter(n == 600) %>% select(n,Setting,power,interim_fraction), by='Setting') %>%
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
    merge(dat_fullsample_reject_AB %>% filter(n == 600), by=c('n','power','MonitoringType','Iteration')) %>%
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
           'RejectionRate_FE_FS' = RejectionRate_FE_FS_FE),
  #n = 2000
  merge(dat_results %>% select(Setting,Iteration,P_Upstrap), 
        dat_setting %>% select(Setting,n,interim_fraction,power) , by = 'Setting') %>%
    filter(n == 2000) %>%
    merge(dat_AB_setting_calibration_corrected %>% filter(interim_fraction != 1.00 & n == 2000), by = c('n','interim_fraction')) %>%
    group_by(Setting,Iteration,MonitoringType) %>% 
    mutate(p_lessthan_fut = ifelse(is.na(FutilityP) == F & (as.numeric(P_Upstrap) < FutilityP), yes = 1, no = 0),
           p_lessthan_eff = ifelse(is.na(EfficacyP) == F & (as.numeric(P_Upstrap) < EfficacyP), yes = 1, no = 0)) %>%
    summarise(proportion_sim_fut = mean(p_lessthan_fut),
              proportion_sim_eff = mean(p_lessthan_eff)) %>%
    merge(dat_setting %>% select(Setting,n,interim_fraction), by = 'Setting') %>%
    merge(dat_AB_setting_calibration_corrected, by = c('n','interim_fraction','MonitoringType')) %>%
    mutate(stop_futility = ifelse(proportion_sim_fut < FutilityProp_multCorr4 & interim_fraction != 0.25, yes = 1, no = 0),
           stop_efficacy = ifelse(proportion_sim_eff > EfficacyProp_multCorr4 & interim_fraction != 0.25, yes = 1, no = 0)) %>%
    select(Setting,Iteration,MonitoringType,stop_futility,stop_efficacy) %>%
    merge(dat_setting %>% filter(n == 2000) %>% select(n,Setting,power,interim_fraction), by='Setting') %>%
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
    merge(dat_fullsample_reject_AB %>% filter(n == 2000), by=c('n','power','MonitoringType','Iteration')) %>%
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
)
dat_ABcalibrated_corrected_M5_upstraped_results_no25 <- rbind(
  #n = 40
  merge(dat_results %>% select(Setting,Iteration,P_Upstrap), 
        dat_setting %>% select(Setting,n,interim_fraction,power) , by = 'Setting') %>%
    filter(n == 40) %>%
    merge(dat_AB_setting_calibration_corrected %>% filter(interim_fraction != 1.00 & n == 40), by = c('n','interim_fraction')) %>%
    group_by(Setting,Iteration,MonitoringType) %>% 
    mutate(p_lessthan_fut = ifelse(is.na(FutilityP) == F & (as.numeric(P_Upstrap) < FutilityP), yes = 1, no = 0),
           p_lessthan_eff = ifelse(is.na(EfficacyP) == F & (as.numeric(P_Upstrap) < EfficacyP), yes = 1, no = 0)) %>%
    summarise(proportion_sim_fut = mean(p_lessthan_fut),
              proportion_sim_eff = mean(p_lessthan_eff)) %>%
    merge(dat_setting %>% select(Setting,n,interim_fraction), by = 'Setting') %>%
    merge(dat_AB_setting_calibration_corrected, by = c('n','interim_fraction','MonitoringType')) %>%
    mutate(stop_futility = ifelse(proportion_sim_fut < FutilityProp_multCorr5 & interim_fraction != 0.25, yes = 1, no = 0),
           stop_efficacy = ifelse(proportion_sim_eff > EfficacyProp_multCorr5 & interim_fraction != 0.25, yes = 1, no = 0)) %>%
    select(Setting,Iteration,MonitoringType,stop_futility,stop_efficacy) %>%
    merge(dat_setting %>% filter(n == 40) %>% select(n,Setting,power,interim_fraction), by='Setting') %>%
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
    merge(dat_fullsample_reject_AB %>% filter(n == 40), by=c('n','power','MonitoringType','Iteration')) %>%
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
           'RejectionRate_FE_FS' = RejectionRate_FE_FS_FE),
  #n = 160
  merge(dat_results %>% select(Setting,Iteration,P_Upstrap), 
        dat_setting %>% select(Setting,n,interim_fraction,power) , by = 'Setting') %>%
    filter(n == 160) %>%
    merge(dat_AB_setting_calibration_corrected %>% filter(interim_fraction != 1.00 & n == 160), by = c('n','interim_fraction')) %>%
    group_by(Setting,Iteration,MonitoringType) %>% 
    mutate(p_lessthan_fut = ifelse(is.na(FutilityP) == F & (as.numeric(P_Upstrap) < FutilityP), yes = 1, no = 0),
           p_lessthan_eff = ifelse(is.na(EfficacyP) == F & (as.numeric(P_Upstrap) < EfficacyP), yes = 1, no = 0)) %>%
    summarise(proportion_sim_fut = mean(p_lessthan_fut),
              proportion_sim_eff = mean(p_lessthan_eff)) %>%
    merge(dat_setting %>% select(Setting,n,interim_fraction), by = 'Setting') %>%
    merge(dat_AB_setting_calibration_corrected, by = c('n','interim_fraction','MonitoringType')) %>%
    mutate(stop_futility = ifelse(proportion_sim_fut < FutilityProp_multCorr5 & interim_fraction != 0.25, yes = 1, no = 0),
           stop_efficacy = ifelse(proportion_sim_eff > EfficacyProp_multCorr5 & interim_fraction != 0.25, yes = 1, no = 0)) %>%
    select(Setting,Iteration,MonitoringType,stop_futility,stop_efficacy) %>%
    merge(dat_setting %>% filter(n == 160) %>% select(n,Setting,power,interim_fraction), by='Setting') %>%
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
    merge(dat_fullsample_reject_AB %>% filter(n == 160), by=c('n','power','MonitoringType','Iteration')) %>%
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
           'RejectionRate_FE_FS' = RejectionRate_FE_FS_FE),
  #n = 600
  merge(dat_results %>% select(Setting,Iteration,P_Upstrap), 
        dat_setting %>% select(Setting,n,interim_fraction,power) , by = 'Setting') %>%
    filter(n == 600) %>%
    merge(dat_AB_setting_calibration_corrected %>% filter(interim_fraction != 1.00 & n == 600), by = c('n','interim_fraction')) %>%
    group_by(Setting,Iteration,MonitoringType) %>% 
    mutate(p_lessthan_fut = ifelse(is.na(FutilityP) == F & (as.numeric(P_Upstrap) < FutilityP), yes = 1, no = 0),
           p_lessthan_eff = ifelse(is.na(EfficacyP) == F & (as.numeric(P_Upstrap) < EfficacyP), yes = 1, no = 0)) %>%
    summarise(proportion_sim_fut = mean(p_lessthan_fut),
              proportion_sim_eff = mean(p_lessthan_eff)) %>%
    merge(dat_setting %>% select(Setting,n,interim_fraction), by = 'Setting') %>%
    merge(dat_AB_setting_calibration_corrected, by = c('n','interim_fraction','MonitoringType')) %>%
    mutate(stop_futility = ifelse(proportion_sim_fut < FutilityProp_multCorr5 & interim_fraction != 0.25, yes = 1, no = 0),
           stop_efficacy = ifelse(proportion_sim_eff > EfficacyProp_multCorr5 & interim_fraction != 0.25, yes = 1, no = 0)) %>%
    select(Setting,Iteration,MonitoringType,stop_futility,stop_efficacy) %>%
    merge(dat_setting %>% filter(n == 600) %>% select(n,Setting,power,interim_fraction), by='Setting') %>%
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
    merge(dat_fullsample_reject_AB %>% filter(n == 600), by=c('n','power','MonitoringType','Iteration')) %>%
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
           'RejectionRate_FE_FS' = RejectionRate_FE_FS_FE),
  #n = 2000
  merge(dat_results %>% select(Setting,Iteration,P_Upstrap), 
        dat_setting %>% select(Setting,n,interim_fraction,power) , by = 'Setting') %>%
    filter(n == 2000) %>%
    merge(dat_AB_setting_calibration_corrected %>% filter(interim_fraction != 1.00 & n == 2000), by = c('n','interim_fraction')) %>%
    group_by(Setting,Iteration,MonitoringType) %>% 
    mutate(p_lessthan_fut = ifelse(is.na(FutilityP) == F & (as.numeric(P_Upstrap) < FutilityP), yes = 1, no = 0),
           p_lessthan_eff = ifelse(is.na(EfficacyP) == F & (as.numeric(P_Upstrap) < EfficacyP), yes = 1, no = 0)) %>%
    summarise(proportion_sim_fut = mean(p_lessthan_fut),
              proportion_sim_eff = mean(p_lessthan_eff)) %>%
    merge(dat_setting %>% select(Setting,n,interim_fraction), by = 'Setting') %>%
    merge(dat_AB_setting_calibration_corrected, by = c('n','interim_fraction','MonitoringType')) %>%
    mutate(stop_futility = ifelse(proportion_sim_fut < FutilityProp_multCorr5 & interim_fraction != 0.25, yes = 1, no = 0),
           stop_efficacy = ifelse(proportion_sim_eff > EfficacyProp_multCorr5 & interim_fraction != 0.25, yes = 1, no = 0)) %>%
    select(Setting,Iteration,MonitoringType,stop_futility,stop_efficacy) %>%
    merge(dat_setting %>% filter(n == 2000) %>% select(n,Setting,power,interim_fraction), by='Setting') %>%
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
    merge(dat_fullsample_reject_AB %>% filter(n == 2000), by=c('n','power','MonitoringType','Iteration')) %>%
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
)
dat_ABcalibrated_corrected_E1_upstraped_results_no25 <- rbind(
  #n = 40
  merge(dat_results %>% select(Setting,Iteration,P_Upstrap), 
        dat_setting %>% select(Setting,n,interim_fraction,power) , by = 'Setting') %>%
    filter(n == 40) %>%
    merge(dat_AB_setting_calibration_corrected %>% filter(interim_fraction != 1.00 & n == 40), by = c('n','interim_fraction')) %>%
    group_by(Setting,Iteration,MonitoringType) %>% 
    mutate(p_lessthan_fut = ifelse(is.na(FutilityP) == F & (as.numeric(P_Upstrap) < FutilityP), yes = 1, no = 0),
           p_lessthan_eff = ifelse(is.na(EfficacyP) == F & (as.numeric(P_Upstrap) < EfficacyP), yes = 1, no = 0)) %>%
    summarise(proportion_sim_fut = mean(p_lessthan_fut),
              proportion_sim_eff = mean(p_lessthan_eff)) %>%
    merge(dat_setting %>% select(Setting,n,interim_fraction), by = 'Setting') %>%
    merge(dat_AB_setting_calibration_corrected, by = c('n','interim_fraction','MonitoringType')) %>%
    mutate(stop_futility = ifelse(proportion_sim_fut < FutilityProp_expCorr1 & interim_fraction != 0.25, yes = 1, no = 0),
           stop_efficacy = ifelse(proportion_sim_eff > EfficacyProp_expCorr1 & interim_fraction != 0.25, yes = 1, no = 0)) %>%
    select(Setting,Iteration,MonitoringType,stop_futility,stop_efficacy) %>%
    merge(dat_setting %>% filter(n == 40) %>% select(n,Setting,power,interim_fraction), by='Setting') %>%
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
    merge(dat_fullsample_reject_AB %>% filter(n == 40), by=c('n','power','MonitoringType','Iteration')) %>%
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
           'RejectionRate_FE_FS' = RejectionRate_FE_FS_FE),
  #n = 160
  merge(dat_results %>% select(Setting,Iteration,P_Upstrap), 
        dat_setting %>% select(Setting,n,interim_fraction,power) , by = 'Setting') %>%
    filter(n == 160) %>%
    merge(dat_AB_setting_calibration_corrected %>% filter(interim_fraction != 1.00 & n == 160), by = c('n','interim_fraction')) %>%
    group_by(Setting,Iteration,MonitoringType) %>% 
    mutate(p_lessthan_fut = ifelse(is.na(FutilityP) == F & (as.numeric(P_Upstrap) < FutilityP), yes = 1, no = 0),
           p_lessthan_eff = ifelse(is.na(EfficacyP) == F & (as.numeric(P_Upstrap) < EfficacyP), yes = 1, no = 0)) %>%
    summarise(proportion_sim_fut = mean(p_lessthan_fut),
              proportion_sim_eff = mean(p_lessthan_eff)) %>%
    merge(dat_setting %>% select(Setting,n,interim_fraction), by = 'Setting') %>%
    merge(dat_AB_setting_calibration_corrected, by = c('n','interim_fraction','MonitoringType')) %>%
    mutate(stop_futility = ifelse(proportion_sim_fut < FutilityProp_expCorr1 & interim_fraction != 0.25, yes = 1, no = 0),
           stop_efficacy = ifelse(proportion_sim_eff > EfficacyProp_expCorr1 & interim_fraction != 0.25, yes = 1, no = 0)) %>%
    select(Setting,Iteration,MonitoringType,stop_futility,stop_efficacy) %>%
    merge(dat_setting %>% filter(n == 160) %>% select(n,Setting,power,interim_fraction), by='Setting') %>%
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
    merge(dat_fullsample_reject_AB %>% filter(n == 160), by=c('n','power','MonitoringType','Iteration')) %>%
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
           'RejectionRate_FE_FS' = RejectionRate_FE_FS_FE),
  #n = 600
  merge(dat_results %>% select(Setting,Iteration,P_Upstrap), 
        dat_setting %>% select(Setting,n,interim_fraction,power) , by = 'Setting') %>%
    filter(n == 600) %>%
    merge(dat_AB_setting_calibration_corrected %>% filter(interim_fraction != 1.00 & n == 600), by = c('n','interim_fraction')) %>%
    group_by(Setting,Iteration,MonitoringType) %>% 
    mutate(p_lessthan_fut = ifelse(is.na(FutilityP) == F & (as.numeric(P_Upstrap) < FutilityP), yes = 1, no = 0),
           p_lessthan_eff = ifelse(is.na(EfficacyP) == F & (as.numeric(P_Upstrap) < EfficacyP), yes = 1, no = 0)) %>%
    summarise(proportion_sim_fut = mean(p_lessthan_fut),
              proportion_sim_eff = mean(p_lessthan_eff)) %>%
    merge(dat_setting %>% select(Setting,n,interim_fraction), by = 'Setting') %>%
    merge(dat_AB_setting_calibration_corrected, by = c('n','interim_fraction','MonitoringType')) %>%
    mutate(stop_futility = ifelse(proportion_sim_fut < FutilityProp_expCorr1 & interim_fraction != 0.25, yes = 1, no = 0),
           stop_efficacy = ifelse(proportion_sim_eff > EfficacyProp_expCorr1 & interim_fraction != 0.25, yes = 1, no = 0)) %>%
    select(Setting,Iteration,MonitoringType,stop_futility,stop_efficacy) %>%
    merge(dat_setting %>% filter(n == 600) %>% select(n,Setting,power,interim_fraction), by='Setting') %>%
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
    merge(dat_fullsample_reject_AB %>% filter(n == 600), by=c('n','power','MonitoringType','Iteration')) %>%
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
           'RejectionRate_FE_FS' = RejectionRate_FE_FS_FE),
  #n = 2000
  merge(dat_results %>% select(Setting,Iteration,P_Upstrap), 
        dat_setting %>% select(Setting,n,interim_fraction,power) , by = 'Setting') %>%
    filter(n == 2000) %>%
    merge(dat_AB_setting_calibration_corrected %>% filter(interim_fraction != 1.00 & n == 2000), by = c('n','interim_fraction')) %>%
    group_by(Setting,Iteration,MonitoringType) %>% 
    mutate(p_lessthan_fut = ifelse(is.na(FutilityP) == F & (as.numeric(P_Upstrap) < FutilityP), yes = 1, no = 0),
           p_lessthan_eff = ifelse(is.na(EfficacyP) == F & (as.numeric(P_Upstrap) < EfficacyP), yes = 1, no = 0)) %>%
    summarise(proportion_sim_fut = mean(p_lessthan_fut),
              proportion_sim_eff = mean(p_lessthan_eff)) %>%
    merge(dat_setting %>% select(Setting,n,interim_fraction), by = 'Setting') %>%
    merge(dat_AB_setting_calibration_corrected, by = c('n','interim_fraction','MonitoringType')) %>%
    mutate(stop_futility = ifelse(proportion_sim_fut < FutilityProp_expCorr1 & interim_fraction != 0.25, yes = 1, no = 0),
           stop_efficacy = ifelse(proportion_sim_eff > EfficacyProp_expCorr1 & interim_fraction != 0.25, yes = 1, no = 0)) %>%
    select(Setting,Iteration,MonitoringType,stop_futility,stop_efficacy) %>%
    merge(dat_setting %>% filter(n == 2000) %>% select(n,Setting,power,interim_fraction), by='Setting') %>%
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
    merge(dat_fullsample_reject_AB %>% filter(n == 2000), by=c('n','power','MonitoringType','Iteration')) %>%
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
)
dat_ABcalibrated_corrected_E2_upstraped_results_no25 <- rbind(
  #n = 40
  merge(dat_results %>% select(Setting,Iteration,P_Upstrap), 
        dat_setting %>% select(Setting,n,interim_fraction,power) , by = 'Setting') %>%
    filter(n == 40) %>%
    merge(dat_AB_setting_calibration_corrected %>% filter(interim_fraction != 1.00 & n == 40), by = c('n','interim_fraction')) %>%
    group_by(Setting,Iteration,MonitoringType) %>% 
    mutate(p_lessthan_fut = ifelse(is.na(FutilityP) == F & (as.numeric(P_Upstrap) < FutilityP), yes = 1, no = 0),
           p_lessthan_eff = ifelse(is.na(EfficacyP) == F & (as.numeric(P_Upstrap) < EfficacyP), yes = 1, no = 0)) %>%
    summarise(proportion_sim_fut = mean(p_lessthan_fut),
              proportion_sim_eff = mean(p_lessthan_eff)) %>%
    merge(dat_setting %>% select(Setting,n,interim_fraction), by = 'Setting') %>%
    merge(dat_AB_setting_calibration_corrected, by = c('n','interim_fraction','MonitoringType')) %>%
    mutate(stop_futility = ifelse(proportion_sim_fut < FutilityProp_expCorr2 & interim_fraction != 0.25, yes = 1, no = 0),
           stop_efficacy = ifelse(proportion_sim_eff > EfficacyProp_expCorr2 & interim_fraction != 0.25, yes = 1, no = 0)) %>%
    select(Setting,Iteration,MonitoringType,stop_futility,stop_efficacy) %>%
    merge(dat_setting %>% filter(n == 40) %>% select(n,Setting,power,interim_fraction), by='Setting') %>%
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
    merge(dat_fullsample_reject_AB %>% filter(n == 40), by=c('n','power','MonitoringType','Iteration')) %>%
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
           'RejectionRate_FE_FS' = RejectionRate_FE_FS_FE),
  #n = 160
  merge(dat_results %>% select(Setting,Iteration,P_Upstrap), 
        dat_setting %>% select(Setting,n,interim_fraction,power) , by = 'Setting') %>%
    filter(n == 160) %>%
    merge(dat_AB_setting_calibration_corrected %>% filter(interim_fraction != 1.00 & n == 160), by = c('n','interim_fraction')) %>%
    group_by(Setting,Iteration,MonitoringType) %>% 
    mutate(p_lessthan_fut = ifelse(is.na(FutilityP) == F & (as.numeric(P_Upstrap) < FutilityP), yes = 1, no = 0),
           p_lessthan_eff = ifelse(is.na(EfficacyP) == F & (as.numeric(P_Upstrap) < EfficacyP), yes = 1, no = 0)) %>%
    summarise(proportion_sim_fut = mean(p_lessthan_fut),
              proportion_sim_eff = mean(p_lessthan_eff)) %>%
    merge(dat_setting %>% select(Setting,n,interim_fraction), by = 'Setting') %>%
    merge(dat_AB_setting_calibration_corrected, by = c('n','interim_fraction','MonitoringType')) %>%
    mutate(stop_futility = ifelse(proportion_sim_fut < FutilityProp_expCorr2 & interim_fraction != 0.25, yes = 1, no = 0),
           stop_efficacy = ifelse(proportion_sim_eff > EfficacyProp_expCorr2 & interim_fraction != 0.25, yes = 1, no = 0)) %>%
    select(Setting,Iteration,MonitoringType,stop_futility,stop_efficacy) %>%
    merge(dat_setting %>% filter(n == 160) %>% select(n,Setting,power,interim_fraction), by='Setting') %>%
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
    merge(dat_fullsample_reject_AB %>% filter(n == 160), by=c('n','power','MonitoringType','Iteration')) %>%
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
           'RejectionRate_FE_FS' = RejectionRate_FE_FS_FE),
  #n = 600
  merge(dat_results %>% select(Setting,Iteration,P_Upstrap), 
        dat_setting %>% select(Setting,n,interim_fraction,power) , by = 'Setting') %>%
    filter(n == 600) %>%
    merge(dat_AB_setting_calibration_corrected %>% filter(interim_fraction != 1.00 & n == 600), by = c('n','interim_fraction')) %>%
    group_by(Setting,Iteration,MonitoringType) %>% 
    mutate(p_lessthan_fut = ifelse(is.na(FutilityP) == F & (as.numeric(P_Upstrap) < FutilityP), yes = 1, no = 0),
           p_lessthan_eff = ifelse(is.na(EfficacyP) == F & (as.numeric(P_Upstrap) < EfficacyP), yes = 1, no = 0)) %>%
    summarise(proportion_sim_fut = mean(p_lessthan_fut),
              proportion_sim_eff = mean(p_lessthan_eff)) %>%
    merge(dat_setting %>% select(Setting,n,interim_fraction), by = 'Setting') %>%
    merge(dat_AB_setting_calibration_corrected, by = c('n','interim_fraction','MonitoringType')) %>%
    mutate(stop_futility = ifelse(proportion_sim_fut < FutilityProp_expCorr2 & interim_fraction != 0.25, yes = 1, no = 0),
           stop_efficacy = ifelse(proportion_sim_eff > EfficacyProp_expCorr2 & interim_fraction != 0.25, yes = 1, no = 0)) %>%
    select(Setting,Iteration,MonitoringType,stop_futility,stop_efficacy) %>%
    merge(dat_setting %>% filter(n == 600) %>% select(n,Setting,power,interim_fraction), by='Setting') %>%
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
    merge(dat_fullsample_reject_AB %>% filter(n == 600), by=c('n','power','MonitoringType','Iteration')) %>%
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
           'RejectionRate_FE_FS' = RejectionRate_FE_FS_FE),
  #n = 2000
  merge(dat_results %>% select(Setting,Iteration,P_Upstrap), 
        dat_setting %>% select(Setting,n,interim_fraction,power) , by = 'Setting') %>%
    filter(n == 2000) %>%
    merge(dat_AB_setting_calibration_corrected %>% filter(interim_fraction != 1.00 & n == 2000), by = c('n','interim_fraction')) %>%
    group_by(Setting,Iteration,MonitoringType) %>% 
    mutate(p_lessthan_fut = ifelse(is.na(FutilityP) == F & (as.numeric(P_Upstrap) < FutilityP), yes = 1, no = 0),
           p_lessthan_eff = ifelse(is.na(EfficacyP) == F & (as.numeric(P_Upstrap) < EfficacyP), yes = 1, no = 0)) %>%
    summarise(proportion_sim_fut = mean(p_lessthan_fut),
              proportion_sim_eff = mean(p_lessthan_eff)) %>%
    merge(dat_setting %>% select(Setting,n,interim_fraction), by = 'Setting') %>%
    merge(dat_AB_setting_calibration_corrected, by = c('n','interim_fraction','MonitoringType')) %>%
    mutate(stop_futility = ifelse(proportion_sim_fut < FutilityProp_expCorr2 & interim_fraction != 0.25, yes = 1, no = 0),
           stop_efficacy = ifelse(proportion_sim_eff > EfficacyProp_expCorr2 & interim_fraction != 0.25, yes = 1, no = 0)) %>%
    select(Setting,Iteration,MonitoringType,stop_futility,stop_efficacy) %>%
    merge(dat_setting %>% filter(n == 2000) %>% select(n,Setting,power,interim_fraction), by='Setting') %>%
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
    merge(dat_fullsample_reject_AB %>% filter(n == 2000), by=c('n','power','MonitoringType','Iteration')) %>%
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
)

###############################################################
### CONSOLIDATE AND STANDARDIZE DATA FORMATTING             ###
###############################################################

#combine data from all AU calibration corrections
dat_final_calibration_correction_AU <- rbind(
  dat_arbitrary_corrected_M1_upstraped_results %>% mutate(Method = 'Arbitrary Upstrap M1'),
  dat_arbitrary_corrected_M2_upstraped_results %>% mutate(Method = 'Arbitrary Upstrap M2'),
  dat_arbitrary_corrected_M3_upstraped_results %>% mutate(Method = 'Arbitrary Upstrap M3'),
  dat_arbitrary_corrected_M4_upstraped_results %>% mutate(Method = 'Arbitrary Upstrap M4'),
  dat_arbitrary_corrected_M5_upstraped_results %>% mutate(Method = 'Arbitrary Upstrap M5'),
  dat_arbitrary_corrected_E1_upstraped_results %>% mutate(Method = 'Arbitrary Upstrap E1'),
  dat_arbitrary_corrected_E2_upstraped_results %>% mutate(Method = 'Arbitrary Upstrap E2'),
  dat_final_output %>% filter(Method %in% c('OBrien Fleming GSD','Arbitrary Upstrap'))
)

dat_final_calibration_correction_AU_no25 <- rbind(
  dat_arbitrary_corrected_M1_upstraped_results_no25 %>% mutate(Method = '0.50-0.75 Arbitrary Upstrap M1'),
  dat_arbitrary_corrected_M2_upstraped_results_no25 %>% mutate(Method = '0.50-0.75 Arbitrary Upstrap M2'),
  dat_arbitrary_corrected_M3_upstraped_results_no25 %>% mutate(Method = '0.50-0.75 Arbitrary Upstrap M3'),
  dat_arbitrary_corrected_M4_upstraped_results_no25 %>% mutate(Method = '0.50-0.75 Arbitrary Upstrap M4'),
  dat_arbitrary_corrected_M5_upstraped_results_no25 %>% mutate(Method = '0.50-0.75 Arbitrary Upstrap M5'),
  dat_arbitrary_corrected_E1_upstraped_results_no25 %>% mutate(Method = '0.50-0.75 Arbitrary Upstrap E1'),
  dat_arbitrary_corrected_E2_upstraped_results_no25 %>% mutate(Method = '0.50-0.75 Arbitrary Upstrap E2'),
  dat_final_output %>% filter(Method == 'OBrien Fleming GSD'),
  dat_final_calibration_no25 %>% filter(Method == '0.50-0.75 Arbitrary Upstrap')
)

#combine data from all CU calibration corrections
dat_final_calibration_correction_CU <- rbind(
  dat_calibrated_corrected_M1_upstraped_results %>% mutate(Method = 'Calibrated Upstrap M1'),
  dat_calibrated_corrected_M2_upstraped_results %>% mutate(Method = 'Calibrated Upstrap M2'),
  dat_calibrated_corrected_M3_upstraped_results %>% mutate(Method = 'Calibrated Upstrap M3'),
  dat_calibrated_corrected_M4_upstraped_results %>% mutate(Method = 'Calibrated Upstrap M4'),
  dat_calibrated_corrected_M5_upstraped_results %>% mutate(Method = 'Calibrated Upstrap M5'),
  dat_calibrated_corrected_E1_upstraped_results %>% mutate(Method = 'Calibrated Upstrap E1'),
  dat_calibrated_corrected_E2_upstraped_results %>% mutate(Method = 'Calibrated Upstrap E2'),
  dat_final_output %>% filter(Method %in% c('OBrien Fleming GSD','Calibrated Upstrap'))
)

dat_final_calibration_correction_CU_no25 <- rbind(
  dat_calibrated_corrected_M1_upstraped_results_no25 %>% mutate(Method = '0.50-0.75 Calibrated Upstrap M1'),
  dat_calibrated_corrected_M2_upstraped_results_no25 %>% mutate(Method = '0.50-0.75 Calibrated Upstrap M2'),
  dat_calibrated_corrected_M3_upstraped_results_no25 %>% mutate(Method = '0.50-0.75 Calibrated Upstrap M3'),
  dat_calibrated_corrected_M4_upstraped_results_no25 %>% mutate(Method = '0.50-0.75 Calibrated Upstrap M4'),
  dat_calibrated_corrected_M5_upstraped_results_no25 %>% mutate(Method = '0.50-0.75 Calibrated Upstrap M5'),
  dat_calibrated_corrected_E1_upstraped_results_no25 %>% mutate(Method = '0.50-0.75 Calibrated Upstrap E1'),
  dat_calibrated_corrected_E2_upstraped_results_no25 %>% mutate(Method = '0.50-0.75 Calibrated Upstrap E2'),
  dat_final_output %>% filter(Method == 'OBrien Fleming GSD'),
  dat_final_calibration_no25 %>% filter(Method == '0.50-0.75 Calibrated Upstrap')
)

#combine data from all GU calibration corrections
dat_final_calibration_correction_GU <- rbind(
  dat_ABcalibrated_corrected_M1_upstraped_results %>% mutate(Method = 'AB Upstrap M1'),
  dat_ABcalibrated_corrected_M2_upstraped_results %>% mutate(Method = 'AB Upstrap M2'),
  dat_ABcalibrated_corrected_M3_upstraped_results %>% mutate(Method = 'AB Upstrap M3'),
  dat_ABcalibrated_corrected_M4_upstraped_results %>% mutate(Method = 'AB Upstrap M4'),
  dat_ABcalibrated_corrected_M5_upstraped_results %>% mutate(Method = 'AB Upstrap M5'),
  dat_ABcalibrated_corrected_E1_upstraped_results %>% mutate(Method = 'AB Upstrap E1'),
  dat_ABcalibrated_corrected_E2_upstraped_results %>% mutate(Method = 'AB Upstrap E2'),
  dat_final_output %>% filter(Method %in% c('OBrien Fleming GSD','Alpha Beta Spending'))
)

dat_final_calibration_correction_GU_no25 <- rbind(
  dat_ABcalibrated_corrected_M1_upstraped_results_no25 %>% mutate(Method = '0.50-0.75 AB Upstrap M1'),
  dat_ABcalibrated_corrected_M2_upstraped_results_no25 %>% mutate(Method = '0.50-0.75 AB Upstrap M2'),
  dat_ABcalibrated_corrected_M3_upstraped_results_no25 %>% mutate(Method = '0.50-0.75 AB Upstrap M3'),
  dat_ABcalibrated_corrected_M4_upstraped_results_no25 %>% mutate(Method = '0.50-0.75 AB Upstrap M4'),
  dat_ABcalibrated_corrected_M5_upstraped_results_no25 %>% mutate(Method = '0.50-0.75 AB Upstrap M5'),
  dat_ABcalibrated_corrected_E1_upstraped_results_no25 %>% mutate(Method = '0.50-0.75 AB Upstrap E1'),
  dat_ABcalibrated_corrected_E2_upstraped_results_no25 %>% mutate(Method = '0.50-0.75 AB Upstrap E2'),
  dat_final_output %>% filter(Method == 'OBrien Fleming GSD'),
  dat_final_calibration_no25 %>% filter(Method == '0.50-0.75 AB Calbirated Upstrap')
)

###############################################################
### SAVE OUTPUT DATASETS                                    ###
###############################################################

save(dat_setting_calibration_corrected, file = './01 Data/Data Processed/dat_setting_calibration_corrected.RData')
#load('./01 Data/Data Processed/dat_setting_calibration_corrected.RData')
save(dat_AB_setting_calibration_corrected, file = './01 Data/Data Processed/dat_AB_setting_calibration_corrected.RData')
#load('./01 Data/Data Processed/dat_AB_setting_calibration_corrected.RData')
save(dat_final_calibration_correction_AU, file = './01 Data/Data Processed/dat_final_calibration_correction_AU.RData')
#load('./01 Data/Data Processed/dat_final_calibration_correction_AU.RData')
save(dat_final_calibration_correction_AU_no25, file = './01 Data/Data Processed/dat_final_calibration_correction_AU_no25.RData')
#load('./01 Data/Data Processed/dat_final_calibration_correction_AU_no25.RData')
save(dat_final_calibration_correction_CU, file = './01 Data/Data Processed/dat_final_calibration_correction_CU.RData')
#load('./01 Data/Data Processed/dat_final_calibration_correction_CU.RData')
save(dat_final_calibration_correction_CU_no25, file = './01 Data/Data Processed/dat_final_calibration_correction_CU_no25.RData')
#load('./01 Data/Data Processed/dat_final_calibration_correction_CU_no25.RData')
save(dat_final_calibration_correction_GU, file = './01 Data/Data Processed/dat_final_calibration_correction_GU.RData')
#load('./01 Data/Data Processed/dat_final_calibration_correction_GU.RData')
save(dat_final_calibration_correction_GU_no25, file = './01 Data/Data Processed/dat_final_calibration_correction_GU_no25.RData')
#load('./01 Data/Data Processed/dat_final_calibration_correction_GU_no25.RData')

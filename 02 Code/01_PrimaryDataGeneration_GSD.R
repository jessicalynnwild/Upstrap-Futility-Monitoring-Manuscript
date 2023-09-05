###############################################################
### file: 01_PrimaryDataGeneration_GSD                      ###
### authors: Jess Wild                                      ###
### creation date:    07/20/2023                            ###
### latest edit date: 07/20/2023                            ###
### description: code to perform GSD interim monitoring     ###
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
load('./01 Data/Data Processed/dat_GSDbounds.RData')

###############################################################
### APPLY GSD CRITERIA TO SIMULATION RESULTS                ###
###############################################################

#apply group sequential stopping boundaries to simulated results
dat_GSD_results <- merge(rbind(dat_results %>% 
                                 select(Setting,Iteration,P_Observed_Full) %>% 
                                 unique() %>% 
                                 filter(Setting %in% 1:16) %>%
                                 mutate(interim_fraction = 1.00) %>%
                                 rename('P_Observed' = P_Observed_Full) %>%
                                 merge(dat_setting %>% select(Setting,n,power) %>% filter(Setting %in% 1:16), by = 'Setting') %>%
                                 select(Iteration,n,interim_fraction,power,P_Observed),
                               dat_results %>% 
                                 select(Setting,Iteration,P_Observed_Int) %>% 
                                 unique() %>%
                                 merge(dat_setting %>% select(Setting,n,power,interim_fraction), by='Setting') %>%
                                 rename('P_Observed' = P_Observed_Int) %>%
                                 select(Iteration,n,interim_fraction,power,P_Observed)),dat_GSDbounds, by = 'interim_fraction') %>%
  mutate(ObservedCriticalValue = sqrt(qchisq(p = as.numeric(as.character(P_Observed)), df=1, lower.tail=F)),
         stop_futility = ifelse(is.na(FutilityCriticalValue)==F & MonitoringType != 'EO' & interim_fraction != 1 & ObservedCriticalValue <= FutilityCriticalValue, yes = 1, no = 0),
         stop_efficacy = ifelse(is.na(FutilityCriticalValue)==F & MonitoringType != 'FO' & ObservedCriticalValue >= EfficacyCriticalValue, yes = 1, no = 
                                  ifelse(MonitoringType == 'FO' & interim_fraction == 1 & ObservedCriticalValue >= EfficacyCriticalValue, yes = 1, no = 0)))

###############################################################
### CALCULATE EVALUATION METRICS                            ###
###############################################################

#use stopping boundary results to get expected sample size and full trial result statistics for calibrated upstrapping
dat_final_results_GSD <- dat_GSD_results %>%
  select(-P_Observed,-FutilityCriticalValue,-EfficacyCriticalValue,-ObservedCriticalValue) %>%
  mutate(interim_fraction = factor(interim_fraction)) %>%
  reshape(idvar = c("n","Iteration","power","MonitoringType","BoundType"), 
          timevar = "interim_fraction", 
          direction = "wide") %>%
  select(-stop_futility.1) %>%
  rename(Reject_FS = stop_efficacy.1) %>%
  mutate(ExpectedN = ifelse(MonitoringType == 'FO' & stop_futility.0.25==1, yes = 0.25*n, no = 
                              ifelse(MonitoringType == 'FO' & stop_futility.0.5==1, yes = 0.5*n, no = 
                                       ifelse(MonitoringType == 'FO' & stop_futility.0.75==1, yes = 0.75*n, no = 
                                                ifelse(MonitoringType == 'EO' & stop_efficacy.0.25==1, yes = 0.25*n, no = 
                                                         ifelse(MonitoringType == 'EO' & stop_efficacy.0.5==1, yes = 0.5*n, no = 
                                                                  ifelse(MonitoringType == 'EO' & stop_efficacy.0.75==1, yes = 0.75*n, no = 
                                                                           ifelse(MonitoringType == 'FE' & stop_efficacy.0.25==1|stop_futility.0.25==1, yes = 0.25*n, no = 
                                                                                    ifelse(MonitoringType == 'FE' & stop_efficacy.0.5==1|stop_futility.0.5==1, yes = 0.5*n, no = 
                                                                                             ifelse(MonitoringType == 'FE' & stop_efficacy.0.75==1|stop_futility.0.75==1, yes = 0.75*n, no = n))))))))),
         Decision = ifelse(ExpectedN<n, yes = 1, no = 0),
         Decision_25 = ifelse(ExpectedN==0.25*n,yes = 1, no = 0),
         Decision_50 = ifelse(ExpectedN==0.5*n,yes = 1, no = 0),
         Decision_75 = ifelse(ExpectedN==0.75*n,yes = 1, no = 0),
         DecisionType = ifelse(MonitoringType == 'FO' & Decision == 1, yes = 'F', no = 
                                 ifelse(MonitoringType == 'EO' & Decision == 1, yes = 'E', no =
                                          ifelse(MonitoringType == 'FE' & Decision == 1 & ExpectedN/n == 0.25 & stop_futility.0.25 == 1, yes = 'F', no = 
                                                   ifelse(MonitoringType == 'FE' & Decision == 1 & ExpectedN/n == 0.50 & stop_futility.0.5 == 1, yes = 'F', no = 
                                                            ifelse(MonitoringType == 'FE' & Decision == 1 & ExpectedN/n == 0.75 & stop_futility.0.75 == 1, yes = 'F', no = 
                                                                     ifelse(MonitoringType == 'FE' & Decision == 1, yes = 'E', no = NA)))))),
         Decision_F = ifelse(is.na(DecisionType) == F & DecisionType == 'F', yes = 1, no = 0),
         Decision_E = ifelse(is.na(DecisionType) == F & DecisionType == 'E', yes = 1, no = 0)) %>%
  mutate(RejectionRate_FS = ifelse((Decision_E == 1)|(Reject_FS == 1 & Decision_F == 0),yes=1,no=0)) %>%
  select(n,power,MonitoringType,BoundType,Iteration,ExpectedN,Decision,DecisionType,Decision_F,
         Decision_E,Decision_25,Decision_50,Decision_75,RejectionRate_FS) %>%
  group_by(n,power,MonitoringType,BoundType) %>%
  summarise(ExpectedN_mean = mean(ExpectedN),
            ExpectedN_sd = sd(ExpectedN),
            Decision = mean(Decision),
            Decision_25 = mean(Decision_25),
            Decision_50 = mean(Decision_50),
            Decision_75 = mean(Decision_75),
            Decision_F = mean(Decision_F),
            Decision_E = mean(Decision_E),
            RejectionRate_FS = mean(RejectionRate_FS)) %>%
  mutate(MonitoringType = factor(MonitoringType)) %>%
  pivot_wider(
    id_cols = c('n','power'),
    names_from = c('MonitoringType','BoundType'),
    values_from = c('ExpectedN_mean','ExpectedN_sd','Decision','Decision_25','Decision_50','Decision_75','Decision_F','Decision_E','RejectionRate_FS')
  ) %>%
  select(-contains(c('Decision_F_FO','Decision_F_EO','Decision_E_FO','Decision_E_EO'))) %>%
  select(n,power,
         ExpectedN_mean_FO_OBF,ExpectedN_sd_FO_OBF,Decision_FO_OBF,Decision_25_FO_OBF,Decision_50_FO_OBF,Decision_75_FO_OBF,RejectionRate_FS_FO_OBF,
         ExpectedN_mean_EO_OBF,ExpectedN_sd_EO_OBF,Decision_EO_OBF,Decision_25_EO_OBF,Decision_50_EO_OBF,Decision_75_EO_OBF,RejectionRate_FS_EO_OBF,
         ExpectedN_mean_FE_OBF,ExpectedN_sd_FE_OBF,Decision_FE_OBF,Decision_25_FE_OBF,Decision_50_FE_OBF,Decision_75_FE_OBF,Decision_F_FE_OBF,Decision_E_FE_OBF,RejectionRate_FS_FE_OBF,
         ExpectedN_mean_FO_PO,ExpectedN_sd_FO_PO,Decision_FO_PO,Decision_25_FO_PO,Decision_50_FO_PO,Decision_75_FO_PO,RejectionRate_FS_FO_PO,
         ExpectedN_mean_EO_PO,ExpectedN_sd_EO_PO,Decision_EO_PO,Decision_25_EO_PO,Decision_50_EO_PO,Decision_75_EO_PO,RejectionRate_FS_EO_PO,
         ExpectedN_mean_FE_PO,ExpectedN_sd_FE_PO,Decision_FE_PO,Decision_F_FE_PO,Decision_E_FE_PO,Decision_25_FE_PO,Decision_50_FE_PO,Decision_75_FE_PO,RejectionRate_FS_FE_PO)

###############################################################
### SAVE OUTPUT DATASETS                                    ###
###############################################################

save(dat_GSD_results, file = './01 Data/Data Processed/dat_GSD_results.RData')
#load('./01 Data/Data Processed/dat_GSD_results.RData')
save(dat_final_results_GSD, file = './01 Data/Data Processed/dat_final_results_GSD.RData')
#load('./01 Data/Data Processed/dat_final_results_GSD.RData')

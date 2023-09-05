###############################################################
### file: 03_DataHarmonization                              ###
### authors: Jess Wild                                      ###
### creation date:    07/20/2022                            ###
### latest edit date: 08/01/2023                            ###
### description: code to consolidate all generated final    ###
### results datasets                                        ###
###############################################################

###############################################################
### SETUP WORKSPACE                                         ###
###############################################################

#clear environment
rm(list = ls())

#set working directory
setwd('C:/Users/wildje/Desktop/DIssertation/Paper 1')

#import packages
library(tidyverse)

#load input datasets
load('./01 Data/Data Processed/dat_pow_err.RData')
load('./01 Data/Data Processed/dat_setting.RData')
load('./01 Data/Data Processed/dat_final_results_arbitrary.RData')
load('./01 Data/Data Processed/dat_final_results_arbitrary_0.25.RData')
load('./01 Data/Data Processed/dat_final_results_arbitrary_0.50.RData')
load('./01 Data/Data Processed/dat_final_results_arbitrary_0.75.RData')
load('./01 Data/Data Processed/dat_final_results_calibrated.RData')
load('./01 Data/Data Processed/dat_final_results_calibrated_0.25.RData')
load('./01 Data/Data Processed/dat_final_results_calibrated_0.50.RData')
load('./01 Data/Data Processed/dat_final_results_calibrated_0.75.RData')
load('./01 Data/Data Processed/dat_final_results_GSD.RData')
load('./01 Data/Data Processed/dat_final_results_AB.RData')
load('./01 Data/Data Processed/dat_final_results_AB_0.25.RData')
load('./01 Data/Data Processed/dat_final_results_AB_0.50.RData')
load('./01 Data/Data Processed/dat_final_results_AB_0.75.RData')
load('./01 Data/Data Processed/dat_final_calibration_no25.RData')
load('./01 Data/Data Processed/dat_final_calibration_correction_AU.RData')
load('./01 Data/Data Processed/dat_final_calibration_correction_CU.RData')
load('./01 Data/Data Processed/dat_final_calibration_correction_GU.RData')
load('./01 Data/Data Processed/dat_final_calibration_correction_AU_no25.RData')
load('./01 Data/Data Processed/dat_final_calibration_correction_CU_no25.RData')
load('./01 Data/Data Processed/dat_final_calibration_correction_GU_no25.RData')
load('./01 Data/Data Processed/dat_condpow_processed.RData')
load('./01 Data/Data Processed/dat_condpow_processed_no25.RData')

###############################################################
### COMBINE ALL SIMULATION RESULTS                          ###
###############################################################

#combine all broad level simulation results
dat_final_output <- rbind(
  dat_final_results_arbitrary %>% mutate(Method = 'Arbitrary Upstrap'),
  dat_final_results_arbitrary_0.25 %>% mutate(Method = 'Arbitrary Upstrap (0.25 Stopping Point)'),
  dat_final_results_arbitrary_0.50 %>% mutate(Method = 'Arbitrary Upstrap (0.50 Stopping Point)'),
  dat_final_results_arbitrary_0.75 %>% mutate(Method = 'Arbitrary Upstrap (0.75 Stopping Point)'),
  dat_final_results_calibrated %>% mutate(Method = 'Calibrated Upstrap'),
  dat_final_results_calibrated_0.25 %>% mutate(Method = 'Calibrated Upstrap (0.25 Stopping Point)'),
  dat_final_results_calibrated_0.50 %>% mutate(Method = 'Calibrated Upstrap (0.50 Stopping Point)'),
  dat_final_results_calibrated_0.75 %>% mutate(Method = 'Calibrated Upstrap (0.75 Stopping Point)'),
  dat_final_results_GSD %>% mutate(Method = 'OBrien Fleming GSD') %>%
    select(n,power,
           ExpectedN_mean_FO_OBF,ExpectedN_sd_FO_OBF,ExpectedN_mean_EO_OBF,ExpectedN_sd_EO_OBF,
           ExpectedN_mean_FE_OBF,ExpectedN_sd_FE_OBF,Decision_FO_OBF,Decision_EO_OBF,Decision_FE_OBF,
           Decision_F_FE_OBF,Decision_E_FE_OBF,Decision_25_FO_OBF,Decision_25_EO_OBF,Decision_25_FE_OBF,
           Decision_50_FO_OBF,Decision_50_EO_OBF,Decision_50_FE_OBF,Decision_75_FO_OBF,Decision_75_EO_OBF,
           Decision_75_FE_OBF,RejectionRate_FS_FO_OBF,RejectionRate_FS_EO_OBF,RejectionRate_FS_FE_OBF,Method) %>%
    rename('ExpectedN_FO_mean' = ExpectedN_mean_FO_OBF,
           'ExpectedN_FO_sd' = ExpectedN_sd_FO_OBF,
           'ExpectedN_EO_mean' = ExpectedN_mean_EO_OBF,
           'ExpectedN_EO_sd' = ExpectedN_sd_EO_OBF,
           'ExpectedN_FE_mean' = ExpectedN_mean_FE_OBF,
           'ExpectedN_FE_sd' = ExpectedN_sd_FE_OBF,
           'DecisionFO' = Decision_FO_OBF,
           'DecisionEO' = Decision_EO_OBF,
           'DecisionFE' = Decision_FE_OBF,
           'DecisionFO_25' = Decision_25_FO_OBF,
           'DecisionEO_25' = Decision_25_EO_OBF,
           'DecisionFE_25' = Decision_25_FE_OBF,
           'DecisionFO_50' = Decision_50_FO_OBF,
           'DecisionEO_50' = Decision_50_EO_OBF,
           'DecisionFE_50' = Decision_50_FE_OBF,
           'DecisionFO_75' = Decision_75_FO_OBF,
           'DecisionEO_75' = Decision_75_EO_OBF,
           'DecisionFE_75' = Decision_75_FE_OBF,
           'DecisionFE_F' = Decision_F_FE_OBF,
           'DecisionFE_E' = Decision_E_FE_OBF,
           'RejectionRate_FO_FS' = RejectionRate_FS_FO_OBF,
           'RejectionRate_EO_FS' = RejectionRate_FS_EO_OBF,
           'RejectionRate_FE_FS' = RejectionRate_FS_FE_OBF),
  dat_final_results_GSD %>% mutate(Method = 'Pocock GSD') %>%
    select(n,power,
           ExpectedN_mean_FO_PO,ExpectedN_sd_FO_PO,ExpectedN_mean_EO_PO,ExpectedN_sd_EO_PO,
           ExpectedN_mean_FE_PO,ExpectedN_sd_FE_PO,Decision_FO_PO,Decision_EO_PO,Decision_FE_PO,
           Decision_F_FE_PO,Decision_E_FE_PO,Decision_25_FO_PO,Decision_25_EO_PO,Decision_25_FE_PO,
           Decision_50_FO_PO,Decision_50_EO_PO,Decision_50_FE_PO,Decision_75_FO_PO,Decision_75_EO_PO,
           Decision_75_FE_PO,RejectionRate_FS_FO_PO,RejectionRate_FS_EO_PO,RejectionRate_FS_FE_PO,Method) %>%
    rename('ExpectedN_FO_mean' = ExpectedN_mean_FO_PO,
           'ExpectedN_FO_sd' = ExpectedN_sd_FO_PO,
           'ExpectedN_EO_mean' = ExpectedN_mean_EO_PO,
           'ExpectedN_EO_sd' = ExpectedN_sd_EO_PO,
           'ExpectedN_FE_mean' = ExpectedN_mean_FE_PO,
           'ExpectedN_FE_sd' = ExpectedN_sd_FE_PO,
           'DecisionFO' = Decision_FO_PO,
           'DecisionEO' = Decision_EO_PO,
           'DecisionFE' = Decision_FE_PO,
           'DecisionFE_F' = Decision_F_FE_PO,
           'DecisionFE_E' = Decision_E_FE_PO,
           'DecisionFO_25' = Decision_25_FO_PO,
           'DecisionEO_25' = Decision_25_EO_PO,
           'DecisionFE_25' = Decision_25_FE_PO,
           'DecisionFO_50' = Decision_50_FO_PO,
           'DecisionEO_50' = Decision_50_EO_PO,
           'DecisionFE_50' = Decision_50_FE_PO,
           'DecisionFO_75' = Decision_75_FO_PO,
           'DecisionEO_75' = Decision_75_EO_PO,
           'DecisionFE_75' = Decision_75_FE_PO,
           'RejectionRate_FO_FS' = RejectionRate_FS_FO_PO,
           'RejectionRate_EO_FS' = RejectionRate_FS_EO_PO,
           'RejectionRate_FE_FS' = RejectionRate_FS_FE_PO),
  dat_final_results_AB  %>% mutate(Method = 'Alpha Beta Spending'),
  dat_final_results_AB_0.25 %>% mutate(Method = 'Alpha Beta Spending (0.25 Stopping Point)'),
  dat_final_results_AB_0.50 %>% mutate(Method = 'Alpha Beta Spending (0.50 Stopping Point)'),
  dat_final_results_AB_0.75 %>% mutate(Method = 'Alpha Beta Spending (0.75 Stopping Point)'),
  dat_final_calibration_no25 %>% filter(Method != 'OBrien Fleming GSD'),
  dat_final_calibration_correction_AU %>% filter(Method != 'OBrien Fleming GSD' & Method != 'Arbitrary Upstrap'),
  dat_final_calibration_correction_CU %>% filter(Method != 'OBrien Fleming GSD' & Method != 'Calibrated Upstrap'),
  dat_final_calibration_correction_GU %>% filter(Method != 'OBrien Fleming GSD' & Method != 'Alpha Beta Spending'),
  dat_final_calibration_correction_AU_no25 %>% filter(Method != 'OBrien Fleming GSD' & Method != '0.50-0.75 Arbitrary Upstrap'),
  dat_final_calibration_correction_CU_no25 %>% filter(Method != 'OBrien Fleming GSD' & Method != '0.50-0.75 Calibrated Upstrap'),
  dat_final_calibration_correction_GU_no25 %>% filter(Method != 'OBrien Fleming GSD' & Method != '0.50-0.75 Alpha Beta Spending')
    ) %>% 
  select(-contains(c('EO','FE'))) %>%
  rbind(dat_condpow_processed,
        dat_condpow_processed_no25)

###############################################################
### COMBINE SIMULATION RESULTS: FIGURE VERSION              ###
###############################################################

#reformat for plotting 
dat_final_output_fig <- dat_final_output %>% 
  filter(power %in% c(0.05,0.80)) %>%
  mutate(Method = factor(Method, 
                         labels = c('AU',
                                    '0.25 AU',
                                    '0.50 AU',
                                    '0.75 AU',
                                    '0.50-0.75 AU',
                                    'M1 AU',
                                    'M2 AU',
                                    'M3 AU',
                                    'M4 AU',
                                    'M5 AU',
                                    'E1 AU',
                                    'E2 AU',
                                    '0.50-0.75 M1 AU',
                                    '0.50-0.75 M2 AU',
                                    '0.50-0.75 M3 AU',
                                    '0.50-0.75 M4 AU',
                                    '0.50-0.75 M5 AU',
                                    '0.50-0.75 E1 AU',
                                    '0.50-0.75 E2 AU',
                                    'CU',
                                    '0.25 CU',
                                    '0.50 CU',
                                    '0.75 CU',
                                    '0.50-0.75 CU',
                                    'M1 CU',
                                    'M2 CU',
                                    'M3 CU',
                                    'M4 CU',
                                    'M5 CU',
                                    'E1 CU',
                                    'E2 CU',
                                    '0.50-0.75 M1 CU',
                                    '0.50-0.75 M2 CU',
                                    '0.50-0.75 M3 CU',
                                    '0.50-0.75 M4 CU',
                                    '0.50-0.75 M5 CU',
                                    '0.50-0.75 E1 CU',
                                    '0.50-0.75 E2 CU',
                                    'GU',
                                    '0.25 GU',
                                    '0.50 GU',
                                    '0.75 GU',
                                    '0.50-0.75 GU',
                                    'M1 GU',
                                    'M2 GU',
                                    'M3 GU',
                                    'M4 GU',
                                    'M5 GU',
                                    'E1 GU',
                                    'E2 GU',
                                    '0.50-0.75 M1 GU',
                                    '0.50-0.75 M2 GU',
                                    '0.50-0.75 M3 GU',
                                    '0.50-0.75 M4 GU',
                                    '0.50-0.75 M5 GU',
                                    '0.50-0.75 E1 GU',
                                    '0.50-0.75 E2 GU',
                                    'OBF',
                                    'PO',
                                    'CP 1%',
                                    'CP 5%',
                                    'CP 10%',
                                    'CP 20%',
                                    '0.50-0.75 CP 1%',
                                    '0.50-0.75 CP 5%',
                                    '0.50-0.75 CP 10%',
                                    '0.50-0.75 CP 20%'), 
                         levels = c('Arbitrary Upstrap',
                                    'Arbitrary Upstrap (0.25 Stopping Point)',
                                    'Arbitrary Upstrap (0.50 Stopping Point)',
                                    'Arbitrary Upstrap (0.75 Stopping Point)',
                                    '0.50-0.75 Arbitrary Upstrap',
                                    'Arbitrary Upstrap M1',
                                    'Arbitrary Upstrap M2',
                                    'Arbitrary Upstrap M3',
                                    'Arbitrary Upstrap M4',
                                    'Arbitrary Upstrap M5',
                                    'Arbitrary Upstrap E1',
                                    'Arbitrary Upstrap E2',
                                    '0.50-0.75 Arbitrary Upstrap M1',
                                    '0.50-0.75 Arbitrary Upstrap M2',
                                    '0.50-0.75 Arbitrary Upstrap M3',
                                    '0.50-0.75 Arbitrary Upstrap M4',
                                    '0.50-0.75 Arbitrary Upstrap M5',
                                    '0.50-0.75 Arbitrary Upstrap E1',
                                    '0.50-0.75 Arbitrary Upstrap E2',
                                    'Calibrated Upstrap',
                                    'Calibrated Upstrap (0.25 Stopping Point)',
                                    'Calibrated Upstrap (0.50 Stopping Point)',
                                    'Calibrated Upstrap (0.75 Stopping Point)',
                                    '0.50-0.75 Calibrated Upstrap',
                                    'Calibrated Upstrap M1',
                                    'Calibrated Upstrap M2',
                                    'Calibrated Upstrap M3',
                                    'Calibrated Upstrap M4',
                                    'Calibrated Upstrap M5',
                                    'Calibrated Upstrap E1',
                                    'Calibrated Upstrap E2',
                                    '0.50-0.75 Calibrated Upstrap M1',
                                    '0.50-0.75 Calibrated Upstrap M2',
                                    '0.50-0.75 Calibrated Upstrap M3',
                                    '0.50-0.75 Calibrated Upstrap M4',
                                    '0.50-0.75 Calibrated Upstrap M5',
                                    '0.50-0.75 Calibrated Upstrap E1',
                                    '0.50-0.75 Calibrated Upstrap E2',
                                    'Alpha Beta Spending',
                                    'Alpha Beta Spending (0.25 Stopping Point)',
                                    'Alpha Beta Spending (0.50 Stopping Point)',
                                    'Alpha Beta Spending (0.75 Stopping Point)',
                                    '0.50-0.75 AB Calbrated Upstrap',
                                    'AB Upstrap M1',
                                    'AB Upstrap M2',
                                    'AB Upstrap M3',
                                    'AB Upstrap M4',
                                    'AB Upstrap M5',
                                    'AB Upstrap E1',
                                    'AB Upstrap E2',
                                    '0.50-0.75 AB Upstrap M1',
                                    '0.50-0.75 AB Upstrap M2',
                                    '0.50-0.75 AB Upstrap M3',
                                    '0.50-0.75 AB Upstrap M4',
                                    '0.50-0.75 AB Upstrap M5',
                                    '0.50-0.75 AB Upstrap E1',
                                    '0.50-0.75 AB Upstrap E2',
                                    'OBrien Fleming GSD',
                                    'Pocock GSD',
                                    'CP 1%',
                                    'CP 5%',
                                    'CP 10%',
                                    'CP 20%',
                                    '0.50-0.75 CP 1%',
                                    '0.50-0.75 CP 5%',
                                    '0.50-0.75 CP 10%',
                                    '0.50-0.75 CP 20%')))

###############################################################
### COMBINE SIMULATION RESULTS: TABULAR VERSION             ###
###############################################################

#reformat for tabular summary 
dat_final_output_tab <- dat_final_output %>% 
  filter(power %in% c(0.05,0.80)) %>%
  mutate(Method = factor(Method, 
                         labels = c('AU',
                                    '0.25 AU',
                                    '0.50 AU',
                                    '0.75 AU',
                                    '0.50-0.75 AU',
                                    'M1 AU',
                                    'M2 AU',
                                    'M3 AU',
                                    'M4 AU',
                                    'M5 AU',
                                    'E1 AU',
                                    'E2 AU',
                                    '0.50-0.75 M1 AU',
                                    '0.50-0.75 M2 AU',
                                    '0.50-0.75 M3 AU',
                                    '0.50-0.75 M4 AU',
                                    '0.50-0.75 M5 AU',
                                    '0.50-0.75 E1 AU',
                                    '0.50-0.75 E2 AU',
                                    'CU',
                                    '0.25 CU',
                                    '0.50 CU',
                                    '0.75 CU',
                                    '0.50-0.75 CU',
                                    'M1 CU',
                                    'M2 CU',
                                    'M3 CU',
                                    'M4 CU',
                                    'M5 CU',
                                    'E1 CU',
                                    'E2 CU',
                                    '0.50-0.75 M1 CU',
                                    '0.50-0.75 M2 CU',
                                    '0.50-0.75 M3 CU',
                                    '0.50-0.75 M4 CU',
                                    '0.50-0.75 M5 CU',
                                    '0.50-0.75 E1 CU',
                                    '0.50-0.75 E2 CU',
                                    'GU',
                                    '0.25 GU',
                                    '0.50 GU',
                                    '0.75 GU',
                                    '0.50-0.75 GU',
                                    'M1 GU',
                                    'M2 GU',
                                    'M3 GU',
                                    'M4 GU',
                                    'M5 GU',
                                    'E1 GU',
                                    'E2 GU',
                                    '0.50-0.75 M1 GU',
                                    '0.50-0.75 M2 GU',
                                    '0.50-0.75 M3 GU',
                                    '0.50-0.75 M4 GU',
                                    '0.50-0.75 M5 GU',
                                    '0.50-0.75 E1 GU',
                                    '0.50-0.75 E2 GU',
                                    'OBF',
                                    'PO',
                                    'CP 1%',
                                    'CP 5%',
                                    'CP 10%',
                                    'CP 20%',
                                    '0.50-0.75 CP 1%',
                                    '0.50-0.75 CP 5%',
                                    '0.50-0.75 CP 10%',
                                    '0.50-0.75 CP 20%',
                                    'FS'), 
                         levels = c('Arbitrary Upstrap',
                                    'Arbitrary Upstrap (0.25 Stopping Point)',
                                    'Arbitrary Upstrap (0.50 Stopping Point)',
                                    'Arbitrary Upstrap (0.75 Stopping Point)',
                                    '0.50-0.75 Arbitrary Upstrap',
                                    'Arbitrary Upstrap M1',
                                    'Arbitrary Upstrap M2',
                                    'Arbitrary Upstrap M3',
                                    'Arbitrary Upstrap M4',
                                    'Arbitrary Upstrap M5',
                                    'Arbitrary Upstrap E1',
                                    'Arbitrary Upstrap E2',
                                    '0.50-0.75 Arbitrary Upstrap M1',
                                    '0.50-0.75 Arbitrary Upstrap M2',
                                    '0.50-0.75 Arbitrary Upstrap M3',
                                    '0.50-0.75 Arbitrary Upstrap M4',
                                    '0.50-0.75 Arbitrary Upstrap M5',
                                    '0.50-0.75 Arbitrary Upstrap E1',
                                    '0.50-0.75 Arbitrary Upstrap E2',
                                    'Calibrated Upstrap',
                                    'Calibrated Upstrap (0.25 Stopping Point)',
                                    'Calibrated Upstrap (0.50 Stopping Point)',
                                    'Calibrated Upstrap (0.75 Stopping Point)',
                                    '0.50-0.75 Calibrated Upstrap',
                                    'Calibrated Upstrap M1',
                                    'Calibrated Upstrap M2',
                                    'Calibrated Upstrap M3',
                                    'Calibrated Upstrap M4',
                                    'Calibrated Upstrap M5',
                                    'Calibrated Upstrap E1',
                                    'Calibrated Upstrap E2',
                                    '0.50-0.75 Calibrated Upstrap M1',
                                    '0.50-0.75 Calibrated Upstrap M2',
                                    '0.50-0.75 Calibrated Upstrap M3',
                                    '0.50-0.75 Calibrated Upstrap M4',
                                    '0.50-0.75 Calibrated Upstrap M5',
                                    '0.50-0.75 Calibrated Upstrap E1',
                                    '0.50-0.75 Calibrated Upstrap E2',
                                    'Alpha Beta Spending',
                                    'Alpha Beta Spending (0.25 Stopping Point)',
                                    'Alpha Beta Spending (0.50 Stopping Point)',
                                    'Alpha Beta Spending (0.75 Stopping Point)',
                                    '0.50-0.75 AB Calbrated Upstrap',
                                    'AB Upstrap M1',
                                    'AB Upstrap M2',
                                    'AB Upstrap M3',
                                    'AB Upstrap M4',
                                    'AB Upstrap M5',
                                    'AB Upstrap E1',
                                    'AB Upstrap E2',
                                    '0.50-0.75 AB Upstrap M1',
                                    '0.50-0.75 AB Upstrap M2',
                                    '0.50-0.75 AB Upstrap M3',
                                    '0.50-0.75 AB Upstrap M4',
                                    '0.50-0.75 AB Upstrap M5',
                                    '0.50-0.75 AB Upstrap E1',
                                    '0.50-0.75 AB Upstrap E2',
                                    'OBrien Fleming GSD',
                                    'Pocock GSD',
                                    'CP 1%',
                                    'CP 5%',
                                    'CP 10%',
                                    'CP 20%',
                                    '0.50-0.75 CP 1%',
                                    '0.50-0.75 CP 5%',
                                    '0.50-0.75 CP 10%',
                                    '0.50-0.75 CP 20%',
                                    'FS'))) %>%
  mutate(ExpectedN_FO_MeanSD = sprintf("%0.00f (%0.02f)", ceiling(ExpectedN_FO_mean), round(ExpectedN_FO_sd,2))) %>%
  select(Method,n,power,ExpectedN_FO_MeanSD,DecisionFO,RejectionRate_FO_FS) %>%
  rbind(dat_pow_err %>%
          merge(dat_setting %>%
                  select(Setting,n,power),
                by = 'Setting') %>%
          select(-Setting) %>%
          rename(RejectionRate_FO_FS = Power) %>%
          filter(power %in% c(0.05,0.80)) %>%
          mutate(Method = as.factor('FS'),
                 ExpectedN_FO_MeanSD = as.character(n),
                 DecisionFO = 0) %>%
          distinct() %>%
          select(Method,n,power,ExpectedN_FO_MeanSD,DecisionFO,RejectionRate_FO_FS)) %>%
  arrange(n,Method) %>% 
  pivot_wider(
    id_cols = c('n','Method'),
    names_from = c('power'),
    values_from = c('ExpectedN_FO_MeanSD','DecisionFO','RejectionRate_FO_FS')
  ) %>%
  mutate(Method = as.factor(Method)) %>%
  select(Method,n,ExpectedN_FO_MeanSD_0.05,DecisionFO_0.05,RejectionRate_FO_FS_0.05,ExpectedN_FO_MeanSD_0.8,DecisionFO_0.8,RejectionRate_FO_FS_0.8)

###############################################################
### COMBINE SIMULATION RESULTS: RELATIVE TABULAR VERSION    ###
###############################################################

#reformat for tabular summary 
dat_final_output_tab_relative <- dat_final_output %>% 
  filter(power %in% c(0.05,0.80)) %>%
  mutate(Method = factor(Method, 
                         labels = c('AU',
                                    '0.25 AU',
                                    '0.50 AU',
                                    '0.75 AU',
                                    '0.50-0.75 AU',
                                    'M1 AU',
                                    'M2 AU',
                                    'M3 AU',
                                    'M4 AU',
                                    'M5 AU',
                                    'E1 AU',
                                    'E2 AU',
                                    '0.50-0.75 M1 AU',
                                    '0.50-0.75 M2 AU',
                                    '0.50-0.75 M3 AU',
                                    '0.50-0.75 M4 AU',
                                    '0.50-0.75 M5 AU',
                                    '0.50-0.75 E1 AU',
                                    '0.50-0.75 E2 AU',
                                    'CU',
                                    '0.25 CU',
                                    '0.50 CU',
                                    '0.75 CU',
                                    '0.50-0.75 CU',
                                    'M1 CU',
                                    'M2 CU',
                                    'M3 CU',
                                    'M4 CU',
                                    'M5 CU',
                                    'E1 CU',
                                    'E2 CU',
                                    '0.50-0.75 M1 CU',
                                    '0.50-0.75 M2 CU',
                                    '0.50-0.75 M3 CU',
                                    '0.50-0.75 M4 CU',
                                    '0.50-0.75 M5 CU',
                                    '0.50-0.75 E1 CU',
                                    '0.50-0.75 E2 CU',
                                    'GU',
                                    '0.25 GU',
                                    '0.50 GU',
                                    '0.75 GU',
                                    '0.50-0.75 GU',
                                    'M1 GU',
                                    'M2 GU',
                                    'M3 GU',
                                    'M4 GU',
                                    'M5 GU',
                                    'E1 GU',
                                    'E2 GU',
                                    '0.50-0.75 M1 GU',
                                    '0.50-0.75 M2 GU',
                                    '0.50-0.75 M3 GU',
                                    '0.50-0.75 M4 GU',
                                    '0.50-0.75 M5 GU',
                                    '0.50-0.75 E1 GU',
                                    '0.50-0.75 E2 GU',
                                    'OBF',
                                    'PO',
                                    'CP 1%',
                                    'CP 5%',
                                    'CP 10%',
                                    'CP 20%',
                                    '0.50-0.75 CP 1%',
                                    '0.50-0.75 CP 5%',
                                    '0.50-0.75 CP 10%',
                                    '0.50-0.75 CP 20%',
                                    'FS'), 
                         levels = c('Arbitrary Upstrap',
                                    'Arbitrary Upstrap (0.25 Stopping Point)',
                                    'Arbitrary Upstrap (0.50 Stopping Point)',
                                    'Arbitrary Upstrap (0.75 Stopping Point)',
                                    '0.50-0.75 Arbitrary Upstrap',
                                    'Arbitrary Upstrap M1',
                                    'Arbitrary Upstrap M2',
                                    'Arbitrary Upstrap M3',
                                    'Arbitrary Upstrap M4',
                                    'Arbitrary Upstrap M5',
                                    'Arbitrary Upstrap E1',
                                    'Arbitrary Upstrap E2',
                                    '0.50-0.75 Arbitrary Upstrap M1',
                                    '0.50-0.75 Arbitrary Upstrap M2',
                                    '0.50-0.75 Arbitrary Upstrap M3',
                                    '0.50-0.75 Arbitrary Upstrap M4',
                                    '0.50-0.75 Arbitrary Upstrap M5',
                                    '0.50-0.75 Arbitrary Upstrap E1',
                                    '0.50-0.75 Arbitrary Upstrap E2',
                                    'Calibrated Upstrap',
                                    'Calibrated Upstrap (0.25 Stopping Point)',
                                    'Calibrated Upstrap (0.50 Stopping Point)',
                                    'Calibrated Upstrap (0.75 Stopping Point)',
                                    '0.50-0.75 Calibrated Upstrap',
                                    'Calibrated Upstrap M1',
                                    'Calibrated Upstrap M2',
                                    'Calibrated Upstrap M3',
                                    'Calibrated Upstrap M4',
                                    'Calibrated Upstrap M5',
                                    'Calibrated Upstrap E1',
                                    'Calibrated Upstrap E2',
                                    '0.50-0.75 Calibrated Upstrap M1',
                                    '0.50-0.75 Calibrated Upstrap M2',
                                    '0.50-0.75 Calibrated Upstrap M3',
                                    '0.50-0.75 Calibrated Upstrap M4',
                                    '0.50-0.75 Calibrated Upstrap M5',
                                    '0.50-0.75 Calibrated Upstrap E1',
                                    '0.50-0.75 Calibrated Upstrap E2',
                                    'Alpha Beta Spending',
                                    'Alpha Beta Spending (0.25 Stopping Point)',
                                    'Alpha Beta Spending (0.50 Stopping Point)',
                                    'Alpha Beta Spending (0.75 Stopping Point)',
                                    '0.50-0.75 AB Calbrated Upstrap',
                                    'AB Upstrap M1',
                                    'AB Upstrap M2',
                                    'AB Upstrap M3',
                                    'AB Upstrap M4',
                                    'AB Upstrap M5',
                                    'AB Upstrap E1',
                                    'AB Upstrap E2',
                                    '0.50-0.75 AB Upstrap M1',
                                    '0.50-0.75 AB Upstrap M2',
                                    '0.50-0.75 AB Upstrap M3',
                                    '0.50-0.75 AB Upstrap M4',
                                    '0.50-0.75 AB Upstrap M5',
                                    '0.50-0.75 AB Upstrap E1',
                                    '0.50-0.75 AB Upstrap E2',
                                    'OBrien Fleming GSD',
                                    'Pocock GSD',
                                    'CP 1%',
                                    'CP 5%',
                                    'CP 10%',
                                    'CP 20%',
                                    '0.50-0.75 CP 1%',
                                    '0.50-0.75 CP 5%',
                                    '0.50-0.75 CP 10%',
                                    '0.50-0.75 CP 20%',
                                    'FS'))) %>%
  select(Method,n,power,ExpectedN_FO_mean,RejectionRate_FO_FS) %>%
  left_join(dat_pow_err %>%
              left_join(dat_setting %>%
                          select(Setting,n,power),
                        by = 'Setting') %>%
              ungroup() %>%
              select(-Setting) %>%
              distinct() %>%
              filter(power %in% c(0.05,0.80)),
            by = c('n','power')) %>%
  mutate(RelativeExpectedN_FO_mean = round(ceiling(ExpectedN_FO_mean)/n,2),
         RelativeRejectionRate_FO_FS = RejectionRate_FO_FS - Power) %>%
  select(-ExpectedN_FO_mean,-RejectionRate_FO_FS,-Power) %>%
  arrange(Method,n) %>% 
  pivot_wider(
    id_cols = c('n','Method'),
    names_from = c('power'),
    values_from = c('RelativeExpectedN_FO_mean','RelativeRejectionRate_FO_FS')
  ) %>%
  mutate(Method = as.factor(Method)) %>%
  select(Method,n,RelativeExpectedN_FO_mean_0.05,RelativeRejectionRate_FO_FS_0.05,RelativeExpectedN_FO_mean_0.8,RelativeRejectionRate_FO_FS_0.8)

###############################################################
### SAVE OUTPUT DATASETS                                    ###
###############################################################

save(dat_final_output, file = './01 Data/Data Processed/dat_final_output.RData')
#load('./01 Data/Data Processed/dat_final_output.RData')
save(dat_final_output_fig, file = './01 Data/Data Processed/dat_final_output_fig.RData')
#load('./01 Data/Data Processed/dat_final_output_fig.RData')
save(dat_final_output_tab, file = './01 Data/Data Processed/dat_final_output_tab.RData')
#load('./01 Data/Data Processed/dat_final_output_tab.RData')
save(dat_final_output_tab_relative, file = './01 Data/Data Processed/dat_final_output_tab_relative.RData')
#load('./01 Data/Data Processed/dat_final_output_tab_relative.RData')

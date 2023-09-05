###############################################################
### file: 01_PrimaryDataGeneration_Calibration              ###
### authors: Jess Wild                                      ###
### creation date:    07/18/2023                            ###
### latest edit date: 07/20/2023                            ###
### description: code to calibrate tuning parameters        ###
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
load('./01 Data/Data Processed/dat_heat.RData')
load('./01 Data/Data Processed/dat_setting.RData')

###############################################################
### PERFORM CALIBRATION OF P-VALUE/PROPORTION THRESHOLDS    ###
###############################################################

#consider three calibration strategies: 
#(1) set p-value <= 0.05 and find proportion such that 
#rejection rate <= 0.05 for the null case and is maximized for the alternative case
#(2) set proportion >= 0.80 and find p-value such that
#rejection rate <= 0.05 for the null case and is maximiezed for the alternative case 
#(3) find p-value/proportion combinations with rejection rate <= 0.05 for the null case 
#and maximized for the alternative case
#each calibration strategy should be considered separately for each planned alternative 
#power case (50%, 80%, 95%)
#each calibration should also be considered for both efficacy (set alpha<=0.05 and maximize beta)
#and futility (set beta>=0.80 and minimize alpha)

#first break the heatmap data up by simulation setting, keeping the null and alternative
#settings for the same sample size/interim fraction connected
dat_calibrate <- dat_heat %>% 
  gather(Setting, RejectionRate, Setting1:Setting48, factor_key=TRUE) %>%
  mutate(Setting = as.numeric(gsub('Setting','',Setting))) %>%
  merge(dat_setting, by = 'Setting') %>%
  mutate(RR = as.factor(ifelse(RR==1, yes = 'NullRejectionRate', no = 'AlternativeRejectionRate'))) %>%
  mutate(RR = str_c(RR,"",power)) %>%
  select(-accrual_type,-outcome_name, -Setting, -p1, -p2, -power) %>%
  spread(RR, RejectionRate)

#next create a table showing results for calibration with arbitrarily chosen proportion value thresholds and fixed p=0.05
dat_calibrate_arbitrary <- dat_calibrate %>%
  subset(P == 0.05 & as.character(Proportion) %in% c('0.05','0.10','0.80','0.95')|P == 0.05 & Proportion %in% c(0.05,0.10,0.80,0.95)) #allow P to be equal to or less than 0.05

#now identify all eligible p-value/proportion combinations for each setting and each 
#calibration strategy (eligibility based on type I error/power in the rejection rate)
dat_calibrate_efficacy_pfixed <- dat_calibrate %>%
  subset(P <= 0.05) %>% #allow P to be equal to or less than 0.05
  subset(NullRejectionRate0.05 <= 0.05) %>%
  group_by(n,interim_fraction) %>%
  mutate(maxAlternativeRejectionRate0.5 = max(AlternativeRejectionRate0.5),
         maxAlternativeRejectionRate0.8 = max(AlternativeRejectionRate0.8),
         maxAlternativeRejectionRate0.95 = max(AlternativeRejectionRate0.95)) %>% 
  ungroup() %>%
  filter(maxAlternativeRejectionRate0.5==AlternativeRejectionRate0.5|maxAlternativeRejectionRate0.8==AlternativeRejectionRate0.8|maxAlternativeRejectionRate0.95==AlternativeRejectionRate0.95) %>%
  mutate(MaxPower0.5 = ifelse(maxAlternativeRejectionRate0.5==AlternativeRejectionRate0.5, yes = 'Yes', no = 'No'),
         MaxPower0.8 = ifelse(maxAlternativeRejectionRate0.8==AlternativeRejectionRate0.8, yes = 'Yes', no = 'No'),
         MaxPower0.95 = ifelse(maxAlternativeRejectionRate0.95==AlternativeRejectionRate0.95, yes = 'Yes', no = 'No')) %>%
  select(-maxAlternativeRejectionRate0.5,-maxAlternativeRejectionRate0.8,-maxAlternativeRejectionRate0.95) %>%
  arrange(n,interim_fraction,MaxPower0.5,MaxPower0.8,MaxPower0.95)

dat_calibrate_futility_pfixed <- dat_calibrate %>%
  subset(P <= 0.05) %>% #allow P to be equal to or less than 0.05
  subset(AlternativeRejectionRate0.8 >= 0.80) %>%
  group_by(n,interim_fraction) %>%
  mutate(minNullRejectionRate0.05 = min(NullRejectionRate0.05)) %>% 
  ungroup() %>%
  filter(minNullRejectionRate0.05==NullRejectionRate0.05) %>%
  select(-minNullRejectionRate0.05,-AlternativeRejectionRate0.5,-AlternativeRejectionRate0.95) %>%
  arrange(n,interim_fraction,AlternativeRejectionRate0.8,NullRejectionRate0.05)

dat_calibrate_efficacy_pfixedstrict <- dat_calibrate %>%
  subset(P == 0.05) %>% #only allow P to be strictly equal to 0.05
  subset(NullRejectionRate0.05 <= 0.05) %>%
  group_by(n,interim_fraction) %>%
  mutate(maxAlternativeRejectionRate0.5 = max(AlternativeRejectionRate0.5),
         maxAlternativeRejectionRate0.8 = max(AlternativeRejectionRate0.8),
         maxAlternativeRejectionRate0.95 = max(AlternativeRejectionRate0.95)) %>% 
  ungroup() %>%
  filter(maxAlternativeRejectionRate0.5==AlternativeRejectionRate0.5|maxAlternativeRejectionRate0.8==AlternativeRejectionRate0.8|maxAlternativeRejectionRate0.95==AlternativeRejectionRate0.95) %>%
  mutate(MaxPower0.5 = ifelse(maxAlternativeRejectionRate0.5==AlternativeRejectionRate0.5, yes = 'Yes', no = 'No'),
         MaxPower0.8 = ifelse(maxAlternativeRejectionRate0.8==AlternativeRejectionRate0.8, yes = 'Yes', no = 'No'),
         MaxPower0.95 = ifelse(maxAlternativeRejectionRate0.95==AlternativeRejectionRate0.95, yes = 'Yes', no = 'No')) %>%
  select(-maxAlternativeRejectionRate0.5,-maxAlternativeRejectionRate0.8,-maxAlternativeRejectionRate0.95) %>%
  arrange(n,interim_fraction,MaxPower0.5,MaxPower0.8,MaxPower0.95)

dat_calibrate_futility_pfixedstrict <- dat_calibrate %>%
  subset(P == 0.05) %>% #allow P to be equal to or less than 0.05
  subset(AlternativeRejectionRate0.8 >= 0.80) %>%
  group_by(n,interim_fraction) %>%
  mutate(minNullRejectionRate0.05 = min(NullRejectionRate0.05)) %>% 
  ungroup() %>%
  filter(minNullRejectionRate0.05==NullRejectionRate0.05) %>%
  select(-minNullRejectionRate0.05,-AlternativeRejectionRate0.5,-AlternativeRejectionRate0.95) %>%
  arrange(n,interim_fraction,AlternativeRejectionRate0.8,NullRejectionRate0.05)

dat_calibrate_efficacy_propfixed <- dat_calibrate %>%
  subset(Proportion >= 0.80) %>% #allow proportion to be greater than or equal to 0.8
  subset(NullRejectionRate0.05 <= 0.05) %>%
  group_by(n,interim_fraction) %>%
  mutate(maxAlternativeRejectionRate0.5 = max(AlternativeRejectionRate0.5),
         maxAlternativeRejectionRate0.8 = max(AlternativeRejectionRate0.8),
         maxAlternativeRejectionRate0.95 = max(AlternativeRejectionRate0.95)) %>% 
  ungroup() %>%
  filter(maxAlternativeRejectionRate0.5==AlternativeRejectionRate0.5|maxAlternativeRejectionRate0.8==AlternativeRejectionRate0.8|maxAlternativeRejectionRate0.95==AlternativeRejectionRate0.95) %>%
  mutate(MaxPower0.5 = ifelse(maxAlternativeRejectionRate0.5==AlternativeRejectionRate0.5, yes = 'Yes', no = 'No'),
         MaxPower0.8 = ifelse(maxAlternativeRejectionRate0.8==AlternativeRejectionRate0.8, yes = 'Yes', no = 'No'),
         MaxPower0.95 = ifelse(maxAlternativeRejectionRate0.95==AlternativeRejectionRate0.95, yes = 'Yes', no = 'No')) %>%
  select(-maxAlternativeRejectionRate0.5,-maxAlternativeRejectionRate0.8,-maxAlternativeRejectionRate0.95) %>%
  arrange(n,interim_fraction,MaxPower0.5,MaxPower0.8,MaxPower0.95)

dat_calibrate_futility_propfixed <- dat_calibrate %>%
  subset(as.numeric(as.character(Proportion)) >= 0.80) %>% #allow P to be equal to or less than 0.05
  subset(AlternativeRejectionRate0.8 >= 0.50) %>% #should be 0.80 but couldn't find any rows meeting criteria
  group_by(n,interim_fraction) %>%
  mutate(minNullRejectionRate0.05 = min(NullRejectionRate0.05)) %>% 
  ungroup() %>%
  filter(minNullRejectionRate0.05==NullRejectionRate0.05) %>%
  select(-minNullRejectionRate0.05,-AlternativeRejectionRate0.5,-AlternativeRejectionRate0.95) %>%
  arrange(n,interim_fraction,AlternativeRejectionRate0.8,NullRejectionRate0.05)

dat_calibrate_efficacy_propfixedstrict <- dat_calibrate %>%
  subset(as.numeric(as.character(Proportion)) %in% c(0.95,0.80,0.10,0.05)) %>%
  subset(NullRejectionRate0.05 <= 0.05) %>%
  group_by(Proportion,n,interim_fraction) %>%
  mutate(maxAlternativeRejectionRate0.5 = max(AlternativeRejectionRate0.5),
         maxAlternativeRejectionRate0.8 = max(AlternativeRejectionRate0.8),
         maxAlternativeRejectionRate0.95 = max(AlternativeRejectionRate0.95)) %>% 
  ungroup() %>%
  filter(maxAlternativeRejectionRate0.5==AlternativeRejectionRate0.5|maxAlternativeRejectionRate0.8==AlternativeRejectionRate0.8|maxAlternativeRejectionRate0.95==AlternativeRejectionRate0.95) %>%
  mutate(MaxPower0.5 = ifelse(maxAlternativeRejectionRate0.5==AlternativeRejectionRate0.5, yes = 'Yes', no = 'No'),
         MaxPower0.8 = ifelse(maxAlternativeRejectionRate0.8==AlternativeRejectionRate0.8, yes = 'Yes', no = 'No'),
         MaxPower0.95 = ifelse(maxAlternativeRejectionRate0.95==AlternativeRejectionRate0.95, yes = 'Yes', no = 'No')) %>%
  select(-maxAlternativeRejectionRate0.5,-maxAlternativeRejectionRate0.8,-maxAlternativeRejectionRate0.95) %>%
  arrange(Proportion,n,interim_fraction,MaxPower0.5,MaxPower0.8,MaxPower0.95)

dat_calibrate_futility_propfixedstrict <- dat_calibrate %>%
  subset(Proportion == 0.80) %>% #allow P to be equal to or less than 0.05
  subset(AlternativeRejectionRate0.8 >= 0.50) %>%
  group_by(n,interim_fraction) %>%
  mutate(minNullRejectionRate0.05 = min(NullRejectionRate0.05)) %>% 
  ungroup() %>%
  filter(minNullRejectionRate0.05==NullRejectionRate0.05) %>%
  select(-minNullRejectionRate0.05,-AlternativeRejectionRate0.5,-AlternativeRejectionRate0.95) %>%
  arrange(n,interim_fraction,AlternativeRejectionRate0.8,NullRejectionRate0.05)

dat_calibrate_efficacy_variable <- dat_calibrate %>%
  subset(NullRejectionRate0.05 <= 0.05) %>%
  group_by(n,interim_fraction) %>%
  mutate(maxAlternativeRejectionRate0.5 = max(AlternativeRejectionRate0.5),
         maxAlternativeRejectionRate0.8 = max(AlternativeRejectionRate0.8),
         maxAlternativeRejectionRate0.95 = max(AlternativeRejectionRate0.95)) %>% 
  ungroup() %>%
  filter(maxAlternativeRejectionRate0.5==AlternativeRejectionRate0.5|maxAlternativeRejectionRate0.8==AlternativeRejectionRate0.8|maxAlternativeRejectionRate0.95==AlternativeRejectionRate0.95) %>%
  mutate(MaxPower0.5 = ifelse(maxAlternativeRejectionRate0.5==AlternativeRejectionRate0.5, yes = 'Yes', no = 'No'),
         MaxPower0.8 = ifelse(maxAlternativeRejectionRate0.8==AlternativeRejectionRate0.8, yes = 'Yes', no = 'No'),
         MaxPower0.95 = ifelse(maxAlternativeRejectionRate0.95==AlternativeRejectionRate0.95, yes = 'Yes', no = 'No')) %>%
  select(-maxAlternativeRejectionRate0.5,-maxAlternativeRejectionRate0.8,-maxAlternativeRejectionRate0.95) %>%
  arrange(n,interim_fraction,MaxPower0.5,MaxPower0.8,MaxPower0.95)

dat_calibrate_futility_variable <- dat_calibrate %>%
  subset(AlternativeRejectionRate0.8 >= 0.80) %>%
  group_by(n,interim_fraction) %>%
  mutate(minNullRejectionRate0.05 = min(NullRejectionRate0.05)) %>% 
  ungroup() %>%
  filter(minNullRejectionRate0.05==NullRejectionRate0.05) %>%
  select(-minNullRejectionRate0.05,-AlternativeRejectionRate0.5,-AlternativeRejectionRate0.95) %>%
  arrange(n,interim_fraction,AlternativeRejectionRate0.8,NullRejectionRate0.05)

###############################################################
### CONSOLIDATE ALL CALIBRATION DATASETS                    ###
###############################################################

#combine and simplify all calibration table results into a single dataframe, then export to excel for easy manual manipulation
dat_calibrate_results <- rbind(
  rbind(
    dat_calibrate_arbitrary %>% 
      mutate(Calibration = 'Arbitrary') %>%
      select(-AlternativeRejectionRate0.5,-AlternativeRejectionRate0.95),
    dat_calibrate_efficacy_pfixed %>% 
      mutate(Calibration = 'Fixed p value') %>%
      filter(MaxPower0.8=='Yes') %>%
      select(-AlternativeRejectionRate0.5,-AlternativeRejectionRate0.95,-MaxPower0.5,-MaxPower0.8,-MaxPower0.95),
    dat_calibrate_efficacy_pfixedstrict %>% 
      mutate(Calibration = 'Strictly fixed p value') %>%
      filter(MaxPower0.8=='Yes') %>%
      select(-AlternativeRejectionRate0.5,-AlternativeRejectionRate0.95,-MaxPower0.5,-MaxPower0.8,-MaxPower0.95),
    dat_calibrate_efficacy_propfixed %>% 
      mutate(Calibration = 'Fixed proportion') %>%
      filter(MaxPower0.8=='Yes') %>%
      select(-AlternativeRejectionRate0.5,-AlternativeRejectionRate0.95,-MaxPower0.5,-MaxPower0.8,-MaxPower0.95),
    dat_calibrate_efficacy_propfixedstrict %>% 
      mutate(Calibration = 'Strictly fixed proportion') %>%
      filter(MaxPower0.8=='Yes') %>%
      select(-AlternativeRejectionRate0.5,-AlternativeRejectionRate0.95,-MaxPower0.5,-MaxPower0.8,-MaxPower0.95),
    dat_calibrate_efficacy_variable %>% 
      mutate(Calibration = 'Variable') %>%
      filter(MaxPower0.8=='Yes') %>%
      select(-AlternativeRejectionRate0.5,-AlternativeRejectionRate0.95,-MaxPower0.5,-MaxPower0.8,-MaxPower0.95)) %>%
    mutate(MonitoringType = 'Efficacy'),
  rbind(
    dat_calibrate_futility_pfixed %>% 
      mutate(Calibration = 'Fixed p value'),
    dat_calibrate_futility_pfixedstrict %>% 
      mutate(Calibration = 'Strictly fixed p value'),
    dat_calibrate_futility_propfixed %>% 
      mutate(Calibration = 'Fixed proportion'),
    dat_calibrate_futility_propfixedstrict %>% 
      mutate(Calibration = 'Strictly fixed proportion'),
    dat_calibrate_futility_variable %>% 
      mutate(Calibration = 'Variable')) %>%
    mutate(MonitoringType = 'Futility')
)

###############################################################
### DEFINE CHOSEN CALIBRATION THRESHOLDS                    ###
###############################################################

#NOTE: upstrapping p value/proportion threshold decision rules are recorded in excel document for convenience

#add threshold decisions to simulation setting dataframe
dat_setting_calibration <- dat_setting %>%
  mutate(P_futility = ifelse((n == 40 & interim_fraction == 0.25), yes = 0.055, no = 0) + ifelse((n == 40 & interim_fraction == 0.50), yes = 0.015, no = 0) + 
           ifelse((n == 40 & interim_fraction == 0.75), yes = 0.09, no = 0) + ifelse((n == 160 & interim_fraction == 0.25), yes = 0.04, no = 0) +
           ifelse((n == 160 & interim_fraction == 0.50), yes = 0.03, no = 0) + ifelse((n == 160 & interim_fraction == 0.75), yes = 0.04, no = 0) + 
           ifelse((n == 600 & interim_fraction == 0.25), yes = 0.06, no = 0) + ifelse((n == 600 & interim_fraction == 0.50), yes = 0.035, no = 0) +
           ifelse((n == 600 & interim_fraction == 0.75), yes = 0.035, no = 0) + ifelse((n == 2000 & interim_fraction == 0.25), yes = 0.05, no = 0) + 
           ifelse((n == 2000 & interim_fraction == 0.50), yes = 0.02, no = 0) + ifelse((n == 2000 & interim_fraction == 0.75), yes = 0.04, no = 0),
         Prop_futility = ifelse((n == 40 & interim_fraction == 0.25), yes = 0.15, no = 0) + ifelse((n == 40 & interim_fraction == 0.50), yes = 0.05, no = 0) + 
           ifelse((n == 40 & interim_fraction == 0.75), yes = 0.4, no = 0) + ifelse((n == 160 & interim_fraction == 0.25), yes = 0.15, no = 0) +
           ifelse((n == 160 & interim_fraction == 0.50), yes = 0.15, no = 0) + ifelse((n == 160 & interim_fraction == 0.75), yes = 0.25, no = 0) + 
           ifelse((n == 600 & interim_fraction == 0.25), yes = 0.2, no = 0) + ifelse((n == 600 & interim_fraction == 0.50), yes = 0.2, no = 0) +
           ifelse((n == 600 & interim_fraction == 0.75), yes = 0.2, no = 0) + ifelse((n == 2000 & interim_fraction == 0.25), yes = 0.2, no = 0) + 
           ifelse((n == 2000 & interim_fraction == 0.50), yes = 0.15, no = 0) + ifelse((n == 2000 & interim_fraction == 0.75), yes = 0.35, no = 0),
         P_efficacy = ifelse((n == 40 & interim_fraction == 0.25), yes = 0, no = 0) + ifelse((n == 40 & interim_fraction == 0.50), yes = 0.01, no = 0) + 
           ifelse((n == 40 & interim_fraction == 0.75), yes = 0.06, no = 0) + ifelse((n == 160 & interim_fraction == 0.25), yes = 0.005, no = 0) +
           ifelse((n == 160 & interim_fraction == 0.50), yes = 0.08, no = 0) + ifelse((n == 160 & interim_fraction == 0.75), yes = 0.05, no = 0) + 
           ifelse((n == 600 & interim_fraction == 0.25), yes = 0.01, no = 0) + ifelse((n == 600 & interim_fraction == 0.50), yes = 0.05, no = 0) +
           ifelse((n == 600 & interim_fraction == 0.75), yes = 0.075, no = 0) + ifelse((n == 2000 & interim_fraction == 0.25), yes = 0.01, no = 0) + 
           ifelse((n == 2000 & interim_fraction == 0.50), yes = 0.025, no = 0) + ifelse((n == 2000 & interim_fraction == 0.75), yes = 0.035, no = 0),
         Prop_efficacy = ifelse((n == 40 & interim_fraction == 0.25), yes = 1, no = 0) + ifelse((n == 40 & interim_fraction == 0.50), yes = 0.4, no = 0) + 
           ifelse((n == 40 & interim_fraction == 0.75), yes = 0.7, no = 0) + ifelse((n == 160 & interim_fraction == 0.25), yes = 0.95, no = 0) +
           ifelse((n == 160 & interim_fraction == 0.50), yes = 0.9, no = 0) + ifelse((n == 160 & interim_fraction == 0.75), yes = 0.65, no = 0) + 
           ifelse((n == 600 & interim_fraction == 0.25), yes = 0.95, no = 0) + ifelse((n == 600 & interim_fraction == 0.50), yes = 0.8, no = 0) +
           ifelse((n == 600 & interim_fraction == 0.75), yes = 0.7, no = 0) + ifelse((n == 2000 & interim_fraction == 0.25), yes = 0.95, no = 0) + 
           ifelse((n == 2000 & interim_fraction == 0.50), yes = 0.7, no = 0) + ifelse((n == 2000 & interim_fraction == 0.75), yes = 0.5, no = 0))

###############################################################
### SAVE OUTPUT DATASETS                                    ###
###############################################################

save(dat_calibrate, file = './01 Data/Data Processed/dat_calibrate.RData')
#load('./01 Data/Data Processed/dat_calibrate.RData')
save(dat_calibrate_arbitrary, file = './01 Data/Data Processed/dat_calibrate_arbitrary.RData')
#load('./01 Data/Data Processed/dat_calibrate_arbitrary.RData')
save(dat_calibrate_efficacy_pfixed, file = './01 Data/Data Processed/dat_calibrate_efficacy_pfixed.RData')
#load('./01 Data/Data Processed/dat_calibrate_efficacy_pfixed.RData')
save(dat_calibrate_futility_pfixed, file = './01 Data/Data Processed/dat_calibrate_futility_pfixed.RData')
#load('./01 Data/Data Processed/dat_calibrate_futility_pfixed.RData')
save(dat_calibrate_efficacy_pfixedstrict, file = './01 Data/Data Processed/dat_calibrate_efficacy_pfixedstrict.RData')
#load('./01 Data/Data Processed/dat_calibrate_efficacy_pfixedstrict.RData')
save(dat_calibrate_futility_pfixedstrict, file = './01 Data/Data Processed/dat_calibrate_futility_pfixedstrict.RData')
#load('./01 Data/Data Processed/dat_calibrate_futility_pfixedstrict.RData')
save(dat_calibrate_efficacy_propfixed, file = './01 Data/Data Processed/dat_calibrate_efficacy_propfixed.RData')
#load('./01 Data/Data Processed/dat_calibrate_efficacy_propfixed.RData')
save(dat_calibrate_futility_propfixed, file = './01 Data/Data Processed/dat_calibrate_futility_propfixed.RData')
#load('./01 Data/Data Processed/dat_calibrate_futility_propfixed.RData')
save(dat_calibrate_efficacy_propfixedstrict, file = './01 Data/Data Processed/dat_calibrate_efficacy_propfixedstrict.RData')
#load('./01 Data/Data Processed/dat_calibrate_efficacy_propfixedstrict.RData')
save(dat_calibrate_futility_propfixedstrict, file = './01 Data/Data Processed/dat_calibrate_futility_propfixedstrict.RData')
#load('./01 Data/Data Processed/dat_calibrate_futility_propfixedstrict.RData')
save(dat_calibrate_efficacy_variable, file = './01 Data/Data Processed/dat_calibrate_efficacy_variable.RData')
#load('./01 Data/Data Processed/dat_calibrate_efficacy_variable.RData')
save(dat_calibrate_futility_variable, file = './01 Data/Data Processed/dat_calibrate_futility_variable.RData')
#load('./01 Data/Data Processed/dat_calibrate_futility_variable.RData')
save(dat_calibrate_results, file = './01 Data/Data Processed/dat_calibrate_results.RData')
#load('./01 Data/Data Processed/dat_calibrate_results.RData')
write.csv(dat_calibrate_results,file = './01 Data/Data Processed/CombinedCalibrationResults.csv')
save(dat_setting_calibration, file = './01 Data/Data Processed/dat_setting_calibration.RData')
#load('./01 Data/Data Processed/dat_setting_calibration.RData')

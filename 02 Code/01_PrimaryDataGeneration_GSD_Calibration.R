###############################################################
### file: 01_PrimaryDataGeneration_GSD_Calibration          ###
### authors: Jess Wild                                      ###
### creation date:    07/20/2023                            ###
### latest edit date: 07/21/2023                            ###
### description: code to perform calibrate GSD proportions  ###
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
load('./01 Data/Data Processed/dat_GSDbounds.RData')
load('./01 Data/Data Processed/dat_results.RData')
load('./01 Data/Data Processed/dat_setting.RData')

###############################################################
###  CALIBRATE ALPHA & BETA SPENDING PROPORTION THRESHOLDS  ###
###############################################################

#for this method we start with the GSD boundary p values and then calibrate the upstrapping proportion 
#threshold to get the best( in terms of max power with type I error <+ 5%) given those p values
#set up list of simulation scenario variables
p <- append(dat_GSDbounds$FutilityPValue,dat_GSDbounds$EfficacyPValue) %>%
  na.omit() %>%
  unique()
prop <- seq(0,1,by=0.05)

#set up results table to hold final heatmap proportions
dat_heat_GSD <- matrix(nrow = (length(p)*length(prop)), ncol = 50) %>%
  as.data.frame() 
colnames(dat_heat_GSD) <- c('P','Proportion','Setting1','Setting2','Setting3',
                            'Setting4','Setting5','Setting6','Setting7','Setting8',
                            'Setting9','Setting10','Setting11','Setting12',
                            'Setting13','Setting14','Setting15','Setting16',
                            'Setting17','Setting18','Setting19','Setting20',
                            'Setting21','Setting22','Setting23','Setting24',
                            'Setting25','Setting26','Setting27','Setting28',
                            'Setting29','Setting30','Setting31','Setting32',
                            'Setting33','Setting34','Setting35','Setting36',
                            'Setting37','Setting38','Setting39','Setting40',
                            'Setting41','Setting42','Setting43','Setting44',
                            'Setting45','Setting46','Setting47','Setting48')
dat_heat_GSD[,1:2] <- expand.grid(p, prop) 

#loop through each setting, p value, and proportion threshold to get proportions
for (setting in 1:48) {
  
  #create progress bar to measure completion of main loop (per simulation setting)
  pb = txtProgressBar(min = 0, max = nrow(dat_heat_GSD), initial = 0, style = 3)
  
  for (index in 1:nrow(dat_heat_GSD)) {
    #define the p value and proportion thresholds to reference
    p <- dat_heat_GSD$P[index]
    prop <- dat_heat_GSD$Proportion[index]
    
    #calculate proportion
    dat_loop <- subset(dat_results,Setting==setting) %>%
      select(P_Upstrap,Iteration)%>% 
      mutate(Reject = ifelse(as.numeric(P_Upstrap) < p,yes = 1, no = 0)) %>%
      group_by(Iteration) %>%
      mutate(Prop = mean(Reject)) %>%
      select(Iteration,Prop) %>%
      unique() %>%
      mutate(Comparison = ifelse(Prop > prop,yes = 1, no = 0))
    prop_loop <- mean(dat_loop$Comparison)
    
    #record proportion
    dat_heat_GSD[index,(setting+2)] <- prop_loop
    
    #update progress bar
    setTxtProgressBar(pb,index)
  }
  #progress update
  marker <- paste('Finished Setting: ',setting, sep = '')
  print(marker)
}

dat_calibrate_GSD <- dat_heat_GSD %>% 
  gather(Setting, RejectionRate, Setting1:Setting48, factor_key=TRUE) %>%
  mutate(Setting = as.numeric(gsub('Setting','',Setting))) %>%
  merge(dat_setting, by = 'Setting') %>%
  mutate(RR = as.factor(ifelse(RR==1, yes = 'NullRejectionRate', no = 'AlternativeRejectionRate'))) %>%
  mutate(RR = str_c(RR,"",power)) %>%
  select(-accrual_type,-outcome_name, -Setting, -p1, -p2, -power) %>%
  spread(RR, RejectionRate) %>%
  select(-AlternativeRejectionRate0.5,-AlternativeRejectionRate0.95) %>%
  subset(P > 0 & P <= 1) %>%
  merge(dat_GSDbounds %>% 
          subset(interim_fraction < 1.00 & BoundType == 'OBF') %>% 
          select(-FutilityCriticalValue,-EfficacyCriticalValue,-BoundType), by = 'interim_fraction') %>%
  mutate(TypeIIError = 1 - AlternativeRejectionRate0.8,
         LevelOfSignificance = 1 - NullRejectionRate0.05) %>%
  subset(P == FutilityPValue | P == EfficacyPValue) %>%
  select(n,interim_fraction,MonitoringType,P,Proportion,FutilityPValue,EfficacyPValue,AlternativeRejectionRate0.8,NullRejectionRate0.05,TypeIIError,LevelOfSignificance) %>%
  arrange(n,interim_fraction,MonitoringType)

dat_calibrate_GSD_futility <- dat_calibrate_GSD %>%
  subset(P == FutilityPValue) %>%
  select(-EfficacyPValue,-P) %>%
  mutate(TypeIIError = 1 - AlternativeRejectionRate0.8,
         LevelOfSignificance = 1 - NullRejectionRate0.05) %>%
  subset(TypeIIError <= 0.20) %>%
  subset(Proportion != 0) %>%
  group_by(n,interim_fraction,MonitoringType) %>%
  filter(LevelOfSignificance == max(LevelOfSignificance)) %>%
  select(n,interim_fraction,MonitoringType,FutilityPValue,Proportion,AlternativeRejectionRate0.8,NullRejectionRate0.05,TypeIIError,LevelOfSignificance) %>%
  arrange(n,interim_fraction,MonitoringType)

dat_calibrate_GSD_efficacy <- dat_calibrate_GSD %>%
  subset(P == EfficacyPValue) %>%
  select(-FutilityPValue,-P) %>%
  subset(NullRejectionRate0.05 <= 0.05) %>%
  subset(Proportion > 0 & Proportion < 1) %>%
  group_by(n,interim_fraction,MonitoringType) %>%
  filter(AlternativeRejectionRate0.8 == max(AlternativeRejectionRate0.8)) %>%
  select(n,interim_fraction,MonitoringType,EfficacyPValue,Proportion,AlternativeRejectionRate0.8,NullRejectionRate0.05) %>%
  arrange(n,interim_fraction,MonitoringType)

###############################################################
### DEFINE CHOSEN GSD CALIBRATION THRESHOLDS                ###
###############################################################

#combine futility and efficacy GSD bounds in reference table
dat_AB_settings <- merge(
  dat_calibrate_GSD_efficacy %>%
    mutate(EfficacyP = EfficacyPValue,
           EfficacyProp = Proportion) %>%
    select(n,interim_fraction,MonitoringType,EfficacyP,EfficacyProp),
  dat_calibrate_GSD_futility %>%
    mutate(FutilityP = FutilityPValue,
           FutilityProp = Proportion) %>%
    select(n,interim_fraction,MonitoringType,FutilityP,FutilityProp),
  by = c('n','interim_fraction','MonitoringType'),
  all = T
) %>%
  distinct(n,interim_fraction,MonitoringType, .keep_all = TRUE) %>%
  rbind(
    data.frame(
      n = c(40,160,600,2000),
      interim_fraction = c(0.25,0.25,0.25,0.25),
      MonitoringType = c('FO','FO','FO','FO'),
      EfficacyP = NA,
      EfficacyProp = NA,
      FutilityP = NA,
      FutilityProp = NA
    )
  ) %>%
  rbind(
    dat_GSDbounds %>% filter(interim_fraction == 1.00 & BoundType == 'OBF') %>% 
      mutate(n = 40, EfficacyProp = NA, FutilityProp = NA, FutilityP = NA) %>%
      mutate(EfficacyP = EfficacyPValue) %>% 
      select(-BoundType,-FutilityCriticalValue,-EfficacyCriticalValue,-FutilityPValue,-EfficacyPValue),
    dat_GSDbounds %>% filter(interim_fraction == 1.00 & BoundType == 'OBF') %>% 
      mutate(n = 160, EfficacyProp = NA, FutilityProp = NA, FutilityP = NA) %>%
      mutate(EfficacyP = EfficacyPValue) %>% 
      select(-BoundType,-FutilityCriticalValue,-EfficacyCriticalValue,-FutilityPValue,-EfficacyPValue),
    dat_GSDbounds %>% filter(interim_fraction == 1.00 & BoundType == 'OBF') %>% 
      mutate(n = 600, EfficacyProp = NA, FutilityProp = NA, FutilityP = NA) %>%
      mutate(EfficacyP = EfficacyPValue) %>% 
      select(-BoundType,-FutilityCriticalValue,-EfficacyCriticalValue,-FutilityPValue,-EfficacyPValue),
    dat_GSDbounds %>% filter(interim_fraction == 1.00 & BoundType == 'OBF') %>% 
      mutate(n = 2000, EfficacyProp = NA, FutilityProp = NA, FutilityP = NA) %>%
      mutate(EfficacyP = EfficacyPValue) %>% 
      select(-BoundType,-FutilityCriticalValue,-EfficacyCriticalValue,-FutilityPValue,-EfficacyPValue)
  )

###############################################################
### SAVE OUTPUT DATASETS                                    ###
###############################################################

save(dat_heat_GSD, file = './01 Data/Data Processed/dat_heat_GSD.RData')
#load('./01 Data/Data Processed/dat_heat_GSD.RData')
save(dat_calibrate_GSD, file = './01 Data/Data Processed/dat_calibrate_GSD.RData')
#load('./01 Data/Data Processed/dat_calibrate_GSD.RData')
write.csv(dat_calibrate_GSD,file = './01 Data/Data Processed/GSDCalibrationData.csv')
save(dat_calibrate_GSD_futility, file = './01 Data/Data Processed/dat_calibrate_GSD_futility.RData')
#load('./01 Data/Data Processed/dat_calibrate_GSD_futility.RData')  
write.csv(dat_calibrate_GSD_futility,file = './01 Data/Data Processed/GSDFutilityCalibrationData.csv')
save(dat_calibrate_GSD_efficacy, file = './01 Data/Data Processed/dat_calibrate_GSD_efficacy.RData')
#load('./01 Data/Data Processed/dat_calibrate_GSD_efficacy.RData')
write.csv(dat_calibrate_GSD_efficacy,file = './01 Data/Data Processed/GSDEfficacyCalibrationData.csv')
save(dat_AB_settings, file = './01 Data/Data Processed/dat_AB_settings.RData')
#load('./01 Data/Data Processed/dat_AB_settings.RData')

###############################################################
### file: 01_PrimaryDataGeneration_ServerVersion            ###
### authors: Jess Wild                                      ###
### creation date: 07/20/2022                               ###
### latest edit date:07/21/2023                             ###
### description: code to run upstrapping simulation study   ###
### on the server                                           ###
###############################################################

###############################################################
### SETUP WORKSPACE                                         ###
###############################################################

#set working directory
setwd('C:/Users/wildje/Desktop/Dissertation/Paper 1')

#import packages
library(tidyverse)
library(randomizeR)

#source helper function script
source(file = './02 Code/00_BasicFunctions.R')

###############################################################
### RUN SIMULATIONS                                         ###
###############################################################

#set up list of simulation scenario variables
n <- c(40,160,600,2000)
power <- c(0.05,0.5,0.8,0.95)
accrual_type <- 'unif'
outcome_name <- 'outcome_constant_high'
interim_fraction <- c(0.25,0.5,0.75)

#create cheatsheet of simulation scenarios
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
save(dat_setting, file = './01 Data/SimulationResultsUpstrap_Updated/SettingsData.RData')

#loop to run all possible simulation combinations
for (i in 1:nrow(dat_setting)) {
  #run each simulation scenario individually by referencing the cheatsheet
  runSimUpstrap(seed = 2022,n_sim = 1000,n_upstrap = 1000,T = 100,p_start = 0.6,p_end = 0.6, #non-variable function parameters
             filename_out = 'C:/Users/wildje/Desktop/Projects/TREATNOW/01 Data/SimulationResultsUpstrap_Updated', #directory path for output
             interim_fraction = dat_setting$interim_fraction[i],n = dat_setting$n[i],RR = dat_setting$RR[i], #variable function parameters
                         accrual_type = dat_setting$accrual_type[i],outcome_name = dat_setting$outcome_name[i],  
             sim_name = paste('Setting',dat_setting$Setting[i],sep = ''))  #identifiable simulation name for output
}

#create empty data frame and add in simulation results
dat_results <- data.frame()
for (i in 1:nrow(dat_setting)) {
  file = paste('./01 Data/SimulationResultsUpstrap_Updated/Setting',i,'_Output.RData',sep = '')
  load(file)
  dat_out <- dat_out %>% mutate(Setting = i)
  dat_results <- rbind(dat_results,dat_out)
}
save(dat_results, file = './01 Data/SimulationResultsUpstrap_Updated/ResultsTotal.RData')

###############################################################
### Calculate Power/Type I Error                            ###
###############################################################

#calculate power (proportion of iterations per setting with p<0.05) 
dat_pow_err <- dat_results %>%
  select(P_Observed_Full, Setting, Iteration) %>%
  unique() %>%
  mutate(Reject = ifelse(as.numeric(P_Observed_Full) < 0.05, yes = 1, no = 0)) %>%
  group_by(Setting) %>%
  mutate(Power = round(mean(Reject),digits=2)) %>%
  select(Setting, Power) %>%
  unique()
save(dat_pow_err, file = './01 Data/SimulationResultsUpstrap_Updated/Power.RData')

###############################################################
### Generate Heatmaps                                       ###
###############################################################

#set up list of simulation scenario variables
p <- seq(0,0.1,by=0.005)
prop <- seq(0,1,by=0.05)

#set up results table to hold final heatmap proportions
dat_heat <- matrix(nrow = (length(p)*length(prop)), ncol = 50) %>%
  as.data.frame() 
colnames(dat_heat) <- c('P','Proportion','Setting1','Setting2','Setting3',
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
dat_heat[,1:2] <- expand.grid(p, prop) 

#loop through each setting, p value, and proportion threshold to get proportions
for (setting in 1:48) {
  
  for (index in 1:nrow(dat_heat)) {
    #define the p value and proportion thresholds to reference
    p <- dat_heat$P[index]
    prop <- dat_heat$Proportion[index]
      
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
    dat_heat[index,(setting+2)] <- prop_loop
    
  }
}
save(dat_heat, file = './01 Data/SimulationResultsUpstrap_Updated/HeatmapData.RData')

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
save(dat_calibrate, file = './01 Data/SimulationResultsUpstrap_Updated/CalibrationData.RData')

#next create a table showing results for calibration with arbitrarily chosen proportion value thresholds and fixed p=0.05
dat_calibrate_arbitrary <- dat_calibrate %>%
  subset(P == 0.05 & as.character(Proportion) %in% c('0.05','0.10','0.80','0.95')|P == 0.05 & Proportion %in% c(0.05,0.10,0.80,0.95)) #allow P to be equal to or less than 0.05
save(dat_calibrate_arbitrary, file = './01 Data/SimulationResultsUpstrap_Updated/ArbitraryCalibrationData.RData')

#now identify all eligible p-value/proportion combinations for each setting and each 
 #calibration strategy (eligibility based on type I error/power in the rejection rate)
dat_calibrate_pfixed <- dat_calibrate %>%
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
save(dat_calibrate_pfixed, file = './01 Data/SimulationResultsUpstrap_Updated/FixedPValueCalibrationData.RData')
dat_calibrate_pfixedstrict <- dat_calibrate %>%
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
save(dat_calibrate_pfixedstrict, file = './01 Data/SimulationResultsUpstrap_Updated/StrictlyFixedPValueCalibrationData.RData')

dat_calibrate_propfixed <- dat_calibrate %>%
  subset(Proportion >= 0.80) %>%
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
save(dat_calibrate_propfixed, file = './01 Data/SimulationResultsUpstrap_Updated/FixedProportionCalibrationData.RData')
dat_calibrate_propfixedstrict <- dat_calibrate %>%
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
save(dat_calibrate_propfixedstrict, file = './01 Data/SimulationResultsUpstrap_Updated/StrictlyFixedProportionCalibrationData.RData')

dat_calibrate_variable <- dat_calibrate %>%
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
save(dat_calibrate_variable, file = './01 Data/SimulationResultsUpstrap_Updated/VariableCalibrationData.RData')
###############################################################
### file: 01_PrimaryDataGeneration_ConditionalPower         ###
### authors: Jess Wild                                      ###
### creation date: 06/26/2023                               ###
### latest edit date: 07/31/2023                            ###
### description: code to apply conditional power to upstrap ###
### data                                                    ###
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
library(randomizeR)
library(rpact) #for group sequential
library(LongCART) #for calculating conditional power

#source helper function script
source(file = './02 Code/TREATNOW_Simulations_Functions.R')

#load simulation settings
load('./01 Data/Data Processed/dat_setting.RData')

#load full sample rejection rate data
load('./01 Data/Data Processed/dat_fullsample_reject.RData')

###############################################################
### CALCULATE CONDITIONAL POWER                             ###
###############################################################

#create empty matrix for overall results
dat_condpow <- matrix(nrow = 0, ncol = 5) %>% as.data.frame()
colnames(dat_condpow) <- c('CondPow','PredProbSucc','Seed','Iteration','Setting')

#loop to run all possible simulation combinations
for (i in 1:nrow(dat_setting)) {
  
  #track loop progress by setting
  print(paste('Setting: ',i))
  
  #define parameters
  seed = 2022
  n_sim = 1000
  n_upstrap = 1000
  T = 100
  p_start = 0.6
  p_end = 0.6 #non-variable function parameters
  interim_fraction = dat_setting$interim_fraction[i]
  n = dat_setting$n[i]
  RR = dat_setting$RR[i] #variable function parameters
  accrual_type = dat_setting$accrual_type[i]
  outcome_name = dat_setting$outcome_name[i]  
  
  #set seed for reproducibility
  set.seed(seed)
  
  #create progress bar to measure completion of main loop
  pb = txtProgressBar(min = 0, max = n_sim, initial = 0, style = 3)
  
  #iterate over the number of simulations
  for (j in 1:n_sim) {
    
    #set seed for reproducibility
    seed_loop <- seed + j
    set.seed(seed_loop)
    
    #simulate data
    dat <- datSim(T,n,RR,p_start,p_end,accrual_type)
    
    #restrict data to desired interim fraction
    dat_interim <- dat[1:(interim_fraction*nrow(dat)),]
    
    #restrict data to desired outcome measure
    dat_interim <- dat_interim %>%
      select(id, outcome_name, randomized, time_accrual) %>%
      rename(outcome = outcome_name)
    
    #apply conditional power to interim data
    freq_table <- dat_interim %>% select(outcome,randomized) %>% table()
    prop_table <- prop.table(freq_table, margin = 1)
    p1 <- prop_table %>% data.frame() %>% filter(outcome == 1 & randomized == 'Treatment') %>% select(Freq) %>% as.numeric() #proportion with outcome = 1 in treatment group
    p2 <- prop_table %>% data.frame() %>% filter(outcome == 1 & randomized == 'Control') %>% select(Freq) %>% as.numeric() #proportion with outcome = 1 in control group
    n1 <- freq_table %>% data.frame() %>% filter(outcome == 1 & randomized == 'Treatment') %>% select(Freq) %>% as.numeric()
    n2 <- freq_table %>% data.frame() %>% filter(outcome == 1 & randomized == 'Control') %>% select(Freq) %>% as.numeric()
    p1 <- ifelse(p1==0|is.na(p1), yes = 0.0000000001, no = p1) #apply small correction to avoid function errors caused by zero inputs
    p2 <- ifelse(p2==0|is.na(p2), yes = 0.0000000001, no = p2)
    n1 <- ifelse(n1==0, yes = 0.0000000001, no = n1)
    n2 <- ifelse(n2==0, yes = 0.0000000001, no = n2)
    exp_diff <- (0.6-0.18046127)*(n==40)+(0.6-0.38017989)*(n==160)+(0.6-0.48629074)*(n==600)+(0.6-0.53797850)*(n==2000) #uses the 80% power scenario's proportion difference (by sample size): see settings dataset for reference values
    results_loop <- capture.output(
      succ_ia(type="bin", nsamples=2, null.value=0, alternative="greater",
            N=n, n=nrow(dat_interim),  a=1,
            propdiff.ia=p2-p1,
            alpha.final = 0.05,
            stderr.ia=sqrt(p1*(1-p1)/n1 + p2*(1-p2)/n2),
            succ.crit="trial", Z.crit.final=1.96,
            propdiff.exp=exp_diff,
            propdiff.prior=exp_diff, sd.prior=sqrt(0.06)
            )
    ) %>% parse_number() %>% t() %>% data.frame() %>% select(X1,X3)
    
    #record conditional power data
    dat_out_loop <- cbind(results_loop,seed,j,i)
    colnames(dat_out_loop) <- c('CondPow','PredProbSucc','Seed','Iteration','Setting')
    dat_condpow <- rbind(dat_condpow,dat_out_loop)
    
    #update progress bar
    setTxtProgressBar(pb,j)
  }
  
  #close progress bar
  close(pb)
}

###############################################################
### PROCESS RAW CONDITIONAL POWER DATA                      ###
###############################################################

#process raw conditional power data
#1%, 5%, 10%, 20% FO
dat_condpow_processed <- dat_condpow %>% 
  merge(dat_setting %>% select(Setting,n,interim_fraction,power), by = 'Setting') %>%
  crossing(CPThreshold = c(0.01,0.05,0.1,0.2)) %>%
  select(-Seed,-Setting,-PredProbSucc) %>%
  mutate(stop_futility = ifelse(CondPow < CPThreshold, yes = 1, no = 0)) %>%
  select(-CondPow) %>%
  mutate(interim_fraction = factor(interim_fraction)) %>%
  as.data.frame() %>%
  reshape(idvar = c("n","Iteration","power","CPThreshold"), 
          timevar = "interim_fraction", 
          direction = "wide") %>%
  select(n,power,Iteration,CPThreshold,stop_futility.0.25,stop_futility.0.5,stop_futility.0.75) %>%
  mutate(ExpectedN_FO = ifelse(stop_futility.0.25==1, yes = 0.25*n, no = 
                                 ifelse(stop_futility.0.5==1, yes = 0.5*n, no = 
                                          ifelse(stop_futility.0.75==1, yes = 0.75*n, no = n))), #case_when
         Decision_FO = ifelse(ExpectedN_FO<n,yes = 1, no = 0),
         Decision_FO_25 = ifelse(ExpectedN_FO==0.25*n,yes = 1, no = 0),
         Decision_FO_50 = ifelse(ExpectedN_FO==0.5*n,yes = 1, no = 0),
         Decision_FO_75 = ifelse(ExpectedN_FO==0.75*n,yes = 1, no = 0)) %>%
  merge(dat_fullsample_reject, by=c('n','power','Iteration')) %>%
  mutate(RejectionRate_FO_FS = ifelse(Decision_FO==0&Reject_FS==1,yes=1,no=0)) %>%
  select(n,power,Iteration,CPThreshold,ExpectedN_FO,Decision_FO,Decision_FO_25,Decision_FO_50,Decision_FO_75,RejectionRate_FO_FS) %>%
  group_by(n,power,CPThreshold) %>%
  summarise(ExpectedN_FO_mean = mean(ExpectedN_FO),
            ExpectedN_FO_sd = sd(ExpectedN_FO),
            DecisionFO = mean(Decision_FO),
            DecisionFO_25 = mean(Decision_FO_25),
            DecisionFO_50 = mean(Decision_FO_50),
            DecisionFO_75 = mean(Decision_FO_75),
            RejectionRate_FO_FS = mean(RejectionRate_FO_FS)) %>%
  mutate(Method = paste('CP',paste(CPThreshold*100,'%',sep = ''))) %>%
  select(-CPThreshold)

###############################################################
### PROCESS RAW CONDITIONAL POWER DATA: DROP 25% Look       ###
###############################################################

#process raw conditional power data
#1%, 5%, 10%, 20% FO
dat_condpow_processed_no25 <- dat_condpow %>% 
  merge(dat_setting %>% select(Setting,n,interim_fraction,power), by = 'Setting') %>%
  crossing(CPThreshold = c(0.01,0.05,0.1,0.2)) %>%
  select(-Seed,-Setting,-PredProbSucc) %>%
  mutate(stop_futility = ifelse(CondPow < CPThreshold, yes = 1, no = 0)) %>%
  select(-CondPow) %>%
  mutate(interim_fraction = factor(interim_fraction)) %>%
  as.data.frame() %>%
  reshape(idvar = c("n","Iteration","power","CPThreshold"), 
          timevar = "interim_fraction", 
          direction = "wide") %>%
  select(n,power,Iteration,CPThreshold,stop_futility.0.5,stop_futility.0.75) %>%
  mutate(ExpectedN_FO = ifelse(stop_futility.0.5==1, yes = 0.5*n, no = 
                                 ifelse(stop_futility.0.75==1, yes = 0.75*n, no = n)), #case_when
         Decision_FO = ifelse(ExpectedN_FO<n,yes = 1, no = 0),
         Decision_FO_50 = ifelse(ExpectedN_FO==0.5*n,yes = 1, no = 0),
         Decision_FO_75 = ifelse(ExpectedN_FO==0.75*n,yes = 1, no = 0)) %>%
  merge(dat_fullsample_reject, by=c('n','power','Iteration')) %>%
  mutate(RejectionRate_FO_FS = ifelse(Decision_FO==0&Reject_FS==1,yes=1,no=0)) %>%
  select(n,power,Iteration,CPThreshold,ExpectedN_FO,Decision_FO,Decision_FO_50,Decision_FO_75,RejectionRate_FO_FS) %>%
  group_by(n,power,CPThreshold) %>%
  summarise(ExpectedN_FO_mean = mean(ExpectedN_FO),
            ExpectedN_FO_sd = sd(ExpectedN_FO),
            DecisionFO = mean(Decision_FO),
            DecisionFO_25 = 0,
            DecisionFO_50 = mean(Decision_FO_50),
            DecisionFO_75 = mean(Decision_FO_75),
            RejectionRate_FO_FS = mean(RejectionRate_FO_FS)) %>%
  mutate(Method = paste('0.50-0.75 CP',paste(CPThreshold*100,'%',sep = ''))) %>%
  select(-CPThreshold)

###############################################################
### SAVE DATASETS                                           ###
###############################################################

#save data
save(dat_condpow, file = './01 Data/Data Processed/dat_condpow.RData')
#load('./01 Data/Data Processed/dat_condpow.RData')
save(dat_condpow_processed, file = './01 Data/Data Processed/dat_condpow_processed.RData')
#load('./01 Data/Data Processed/dat_condpow_processed.RData')
save(dat_condpow_processed_no25, file = './01 Data/Data Processed/dat_condpow_processed_no25.RData')
#load('./01 Data/Data Processed/dat_condpow_processed_no25.RData')

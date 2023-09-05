###############################################################
### file: 01_PrimaryDataGeneration_SimulationImplementation ###
### authors: Jess Wild                                      ###
### creation date:    07/18/2023                            ###
### latest edit date: 08/19/2023                            ###
### description: code to run upstrapping simulation study   ###
### locally                                                 ###
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

#source helper function script
source(file = './02 Code/00_BasicFunctions.R')

#load input datasets
load('./01 Data/Data Processed/dat_setting.RData')

###############################################################
### DEFINE HELPER FUNCTIONS                                 ###
###############################################################

#function to perform calibration simulation for upstrapping
runSimUpstrap <- function(seed,n_sim,n_upstrap,interim_fraction,T,n,RR,p_start,p_end,accrual_type,outcome_name,sim_name,filename_out){
  #seed: seed passed to each simulation for reproducibility
  #n_sim: number of total simulations to run for the other given parameters
  #n_upstrap: number of upstrapped datasets to generate per simulation
  #interim_fraction: proportion of accumulated data available for the interim analysis
  #T: duration of study (in days)
  #n: total sample size
  #RR: relative risk of treatment group compared to control group
  #p_start: probability of outcome in the control group (between 0,1) at lower end of probability range
  #p_end: probability of outcome in the control group (between 0,1) at higher end of probability range
  #accrual_type: choice of accrual pattern from (unif,beta_a,beta_b,beta_c)
  #beta_a: more accrual early in the trial
  #beta_b: more accrual later in the trial
  #beta_c: bimodal shape with two peaks of accrual
  #outcome_name: choice of efficacy time trend from (outcome_constant_low,outcome_constant_high,outcome_linear,outcome_linearNEG,outcome_quad,outcome_quadNEG)
  #outcome_constant_low: constant low probability of success 
  #constant high: constant high probability of success 
  #outcome_linear: linearly increasing probability of success
  #outcome_linearNEG: linearly decreasing probability of success
  #outcome_quad: quadratic concave probability of success
  #outcome_quadNEG: quadratic convex probability of success
  #sim_name: identifiable name for specific simulation scenario (used in naming output files)
  #filename_out: directory path to save results files
  
  #set seed for reproducibility
  set.seed(seed)
  
  #create empty matrix for overall results
  dat_out <- matrix(nrow = 0, ncol = 8) %>% as.data.frame()
  colnames(dat_out) <- c('P_Upstrap','P_Upstrap_TestUsed','P_Observed_Int','P_Observed_Int_TestUsed','P_Observed_Full','P_Observed_Full_TestUsed','Seed','Iteration')
  
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
    dat_full <- dat %>%
      select(id, outcome_name, randomized, time_accrual) %>%
      rename(outcome = outcome_name)
    
    #run complete data LR models
    mod_CSFE_Comp <- model_CSFE_simp(dat_current = dat_full, outcome_name = 'outcome')
    
    #create empty dataframe to hold simulation values
    dat_out_loop <- matrix(nrow = n_upstrap, ncol = 8) %>% as.data.frame()
    colnames(dat_out_loop) <- c('P_Upstrap','P_Upstrap_TestUsed','P_Observed_Int','P_Observed_Int_TestUsed','P_Observed_Full','P_Observed_Full_TestUsed','Seed','Iteration')
    
    #iterate over the number of upstrapped datasets
    for (i in 1:n_upstrap) {
      
      #run interim LR models
      mod_CSFE_Int <- model_CSFE_simp(dat_current = dat_interim, outcome_name = 'outcome')
      
      #run uppstrapped LR models
      mod_CSFE_UP <- model_CSFE_uppstrap_simp(dat_current = dat_interim, outcome_name = 'outcome', n_total = n)
      
      #collect and standardize output
      dat_out_loop[i,] <- c(mod_CSFE_UP$P_treatment,mod_CSFE_UP$TestUsed,
                            mod_CSFE_Int$P_treatment,mod_CSFE_Int$TestUsed,
                            mod_CSFE_Comp$P_treatment,mod_CSFE_Comp$TestUsed,
                            seed,j)
    }
    dat_out <- rbind(dat_out,dat_out_loop)
    
    #update progress bar
    setTxtProgressBar(pb,j)
  }
  
  #save output dataset in the appropriate file directory
  file_out_data <- paste(paste(paste(filename_out,sim_name,sep = '/'),'_Output',sep = ''),'.Rdata',sep = '')
  save(dat_out,file = file_out_data)
  
  #close progress bar
  close(pb)
  
}

###############################################################
### RUN SIMULATIONS                                         ###
###############################################################

#loop to run all possible simulation combinations
for (i in 1:nrow(dat_setting)) {
  print('Starting setting:')
  print(dat_setting$Setting[i])
  #run each simulation scenario individually by referencing the cheatsheet
  runSimUpstrap(seed = 2022,n_sim = 1000,n_upstrap = 1000,T = 100,p_start = 0.6,p_end = 0.6, #non-variable function parameters
                filename_out = 'C:/Users/wildje/Desktop/Projects/TREATNOW/01 Data/Upstrap 10-3-2022 Results', #directory path for output
                interim_fraction = dat_setting$interim_fraction[i],n = dat_setting$n[i],RR = dat_setting$RR[i], #variable function parameters
                accrual_type = dat_setting$accrual_type[i],outcome_name = dat_setting$outcome_name[i],  
                sim_name = paste('Setting',dat_setting$Setting[i],sep = ''))  #identifiable simulation name for output
}

#create empty data frame and add in simulation results
dat_results <- data.frame()

#create progress bar to measure completion of main loop
pb = txtProgressBar(min = 0, max = nrow(dat_setting), initial = 0, style = 3)

#loop through to process and combine datafiles
for (i in 1:nrow(dat_setting)) {
  file = paste('./01 Data/Upstrap 10-3-2022 Results/Setting',i,'_Output.RData',sep = '')
  load(file)
  dat_out <- dat_out %>% mutate(Setting = i)
  dat_results <- rbind(dat_results,dat_out)
  
  #update progress bar
  setTxtProgressBar(pb,i)
}

###############################################################
### SAVE OUTPUT DATASETS                                    ###
###############################################################

save(dat_results, file = './01 Data/Data Processed/dat_results.RData')
#load('./01 Data/Data Processed/dat_results.RData')
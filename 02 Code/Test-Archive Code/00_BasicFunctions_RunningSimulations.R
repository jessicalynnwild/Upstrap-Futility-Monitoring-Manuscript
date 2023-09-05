###############################################################
### file: 00_BasicFunctions _RunningSimulations             ###
### authors: Jess Wild                                      ###
### creation date:    07/18/2023                            ###
### latest edit date: 07/18/2023                            ###
### description: Code for implementing simulations/data     ### 
###   generation                                            ###
###############################################################

###############################################################
### FUNCTIONS FOR RUNNING SIMULATION SCENARIOS:             ###
###############################################################


runSimFreq <- function(seed,n_sim,interim_fraction,T,n,RR,p_start,p_end,accrual_type,outcome_name,sim_name,filename_out){
  #seed: seed passed to each simulation for reproducibility
  #n_sim: number of total simulations to run for the other given parameters
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
  
  #create empty dataframe to hold simulation values
  dat_out <- matrix(nrow = (n_sim*12), ncol = 6) %>% as.data.frame()
  colnames(dat_out) <- c('Value','ValueName','Covariates','Model','Simulation','Seed')
  
  #loop through all simulations
  for (i in 1:n_sim) {
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
    
    #run LR models
    mod_LR <- model_LR(dat_current = dat_interim, outcome_name = 'outcome')
    mod_LR_comp <- model_LR(dat_current = dat_full, outcome_name = 'outcome') %>% 
      mutate(Model = c('LR (Complete Data)','LR (Complete Data)','LR (Complete Data)'))
    
    #run uppstrapped LR models
    mod_LR_UP <- model_LR_uppstrap(dat_current = dat_interim, outcome_name = 'outcome', n_total = n)
    
    #run weighted upstrapped LR models
    mod_LR_UP_WT <- model_LR_uppstrap_wt(dat_current = dat_interim, outcome_name = 'outcome', n_total = n)
    
    #collect and standardize output
    dat_LR <- mod_LR %>%
      mutate(Value = P_treatment,
             ValueName = rep('P_treatment',nrow(mod_LR)),
             Simulation = rep(i,nrow(mod_LR)),
             Seed = rep(seed,nrow(mod_LR))) %>%
      select(Value,ValueName,Covariates,Model,Simulation,Seed)
    dat_LR_comp <- mod_LR_comp %>%
      mutate(Value = P_treatment,
             ValueName = rep('P_treatment',nrow(mod_LR_comp)),
             Simulation = rep(i,nrow(mod_LR_comp)),
             Seed = rep(seed,nrow(mod_LR_comp))) %>%
      select(Value,ValueName,Covariates,Model,Simulation,Seed)
    dat_LR_UP <- mod_LR_UP %>%
      mutate(Value = P_treatment,
             ValueName = rep('P_treatment',nrow(mod_LR_UP)),
             Simulation = rep(i,nrow(mod_LR_UP)),
             Seed = rep(seed,nrow(mod_LR_UP))) %>%
      select(Value,ValueName,Covariates,Model,Simulation,Seed)
    dat_LR_UP_WT <- mod_LR_UP_WT %>%
      mutate(Value = P_treatment,
             ValueName = rep('P_treatment',nrow(mod_LR_UP_WT)),
             Simulation = rep(i,nrow(mod_LR_UP_WT)),
             Seed = rep(seed,nrow(mod_LR_UP_WT))) %>%
      select(Value,ValueName,Covariates,Model,Simulation,Seed)
    dat_out[(1+(i-1)*12):(i*12),] <- rbind(dat_LR,dat_LR_comp,dat_LR_UP,dat_LR_UP_WT)
  }
  
  #save output dataset in the appropriate file directory
  file_out_data <- paste(paste(paste(filename_out,sim_name,sep = '/'),'_Output',sep = ''),'.Rdata',sep = '')
  save(dat_out,file = file_out_data)
}

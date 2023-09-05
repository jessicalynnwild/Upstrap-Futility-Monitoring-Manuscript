###############################################################
### file: 00_BasicFunctions_BayesianFunctions               ###
### authors: Jess Wild                                      ###
### creation date:    07/18/2023                            ###
### latest edit date: 07/18/2023                            ###
### description: All basic function code for specifically   ###
###   Bayesian analysis                                     ###
###############################################################

###############################################################
### FUNCTIONS TO PERFORM BAYESIAN MODELING:                 ###
###############################################################

#function to apply bayesian methods to calculate posterior probability and predicted probability of success given an interim dataset
model_bayes <- function(dat_current, outcome_name, n_total){
  #function arguments:
  #dat_current: simulated dataset
  #n_total: total desired sample size
  #outcome_name: selected outcome column from simulated dataset.  Must be passed as a string (i.e. 'outcome_linear')
  
  #define number of subjects left to enroll (calculated for treatment and control groups separately)
  n <- nrow(dat_current)
  n_trt <- nrow(subset(dat_current, randomized == 'Treatment'))
  n_ctr <- nrow(subset(dat_current, randomized == 'Control'))
  n_remaining <- n_total - n #CHECK IF WE NEED THIS
  n_remaining_trt <- n_total/2 - n_trt
  n_remaining_ctr <- n_total/2 - n_ctr
  
  #streamline the dataset
  dat_model <- dat_current %>%
    select(id, outcome_name, randomized, time_accrual) %>%
    mutate(time_accrual_sq = time_accrual^2) %>%
    rename(outcome = outcome_name)
  
  #calculate posterior probability for observed data
  mod1 <- brm(outcome ~ randomized, refresh=0, family='bernoulli', data=dat_model)
  asd1 <- posterior_samples(mod1, pars='b_randomized')
  postprob_trt1 <- mean(asd1$b_randomized > 0) #calculate posterior probability that the log-odds is >0 for treatment
  mod2 <- brm(outcome ~ randomized + time_accrual, refresh=0, family='bernoulli', data=dat_model)
  asd2 <- posterior_samples(mod2, pars='b_randomized')
  postprob_trt2 <- mean(asd2$b_randomized > 0) #calculate posterior probability that the log-odds is >0 for treatment
  mod3 <- brm(outcome ~ randomized + time_accrual + time_accrual_sq, refresh=0, family='bernoulli', data=dat_model)
  asd3 <- posterior_samples(mod3, pars='b_randomized')
  postprob_trt3 <- mean(asd3$b_randomized > 0) #calculate posterior probability that the log-odds is >0 for treatment
  
  #create data frame with remaining participants to enroll
  pred_dat <- data.frame(randomized = c(rep('Treatment',n_remaining_trt),rep('Control',n_remaining_ctr)), 
                         time_accrual = c(seq(max(dat_model$time_accrual),360,length.out=n_remaining_trt),
                                          seq(max(dat_model$time_accrual),360,length.out=n_remaining_ctr)) ) %>%
    mutate(time_accrual_sq = time_accrual^2)
  
  # generate remaining data based on the predictive posterior distribution
  mod1Pred <- posterior_predict(newdata=pred_dat, object=mod1, allow_new_levels=TRUE, re_formula=NULL, nsamples=1000)
  mod2Pred <- posterior_predict(newdata=pred_dat, object=mod2, allow_new_levels=TRUE, re_formula=NULL, nsamples=1000)
  mod3Pred <- posterior_predict(newdata=pred_dat, object=mod3, allow_new_levels=TRUE, re_formula=NULL, nsamples=1000)
  
  #set up a matrix to store posterior probability and frequentist p values for the 1000 predicted data sets (bayesian models commented out due to frequent crashes)
  #results_pred <- data.frame(PP_mod1 = rep(NA,1000),Pval_mod1 = rep(NA,1000),PP_mod2 = rep(NA,1000),Pval_mod2 = rep(NA,1000))
  results_pred <- data.frame(Pval_mod1 = rep(NA,1000),Pval_mod2 = rep(NA,1000),Pval_mod3 = rep(NA,1000))
  
  #calculate PP and frequentist p value for each of the 1000 predicted datasets
  for (i in 1:1000) {
    #set up the dataset
    dat_pp_mod1 <- rbind(dat_model[,2:5], cbind(outcome=mod1Pred[i,], pred_dat))
    dat_pp_mod2 <- rbind(dat_model[,2:5], cbind(outcome=mod2Pred[i,], pred_dat))
    
    #fit relevant models (bayesian models commented out due to frequent crashes)
    #mod1_loop <- brm(outcome ~ randomized, refresh=0, family='bernoulli', data=dat_pp_mod1)
    mod1_loop_freq <- glm(outcome ~ randomized, family='binomial', data=dat_pp_mod1)
    #mod2_loop <- brm(outcome ~ randomized + time_accrual, refresh=0, family='bernoulli', data=dat_pp_mod2)
    mod2_loop_freq <- glm(outcome ~ randomized + time_accrual, family='binomial', data=dat_pp_mod2)
    #mod3_loop <- brm(outcome ~ randomized + time_accrual + time_accrual_sq, refresh=0, family='bernoulli', data=dat_pp_mod2)
    mod3_loop_freq <- glm(outcome ~ randomized + time_accrual + time_accrual_sq, family='binomial', data=dat_pp_mod2)
    
    #calculate posterior probability (bayesian models commented out due to frequent crashes)
    #asd_loop1 <- posterior_samples(mod1_loop, pars='b_randomized')
    #postprob_trt_loop1 <- mean(asd_loop1$b_randomized > 0) # calculate posterior probability that the log-odds is >0 for treatment
    #asd_loop2 <- posterior_samples(mod2_loop, pars='b_randomized')
    #postprob_trt_loop2 <- mean(asd_loop2$b_randomized > 0)
    #asd_loop3 <- posterior_samples(mod3_loop, pars='b_randomized')
    #postprob_trt_loop3 <- mean(asd_loop3$b_randomized > 0)
    
    #record PP/frequentist p values (bayesian models commented out due to frequent crashes)
    #results_pred[i,] <- c(postprob_trt_loop1,coef(summary(mod1_loop_freq))[2,4],postprob_trt_loop2,coef(summary(mod2_loop_freq))[2,4])
    results_pred[i,] <- c(coef(summary(mod1_loop_freq))[2,4],coef(summary(mod2_loop_freq))[2,4],coef(summary(mod3_loop_freq))[2,4])
    
  }
 
  #calculate PPOS and frequentist PPOS (bayesian models commented out due to frequent crashes)
  #PPOS_trt_mod1 <- sum(results_pred$PP_mod1 > 0.975)/nrow(results_pred)
  PPOS_freq_trt_mod1 <- sum(results_pred$Pval_mod1 < 0.025)/nrow(results_pred)
  #PPOS_trt_mod2 <- sum(results_pred$PP_mod2 > 0.975)/nrow(results_pred)
  PPOS_freq_trt_mod2 <- sum(results_pred$Pval_mod2 < 0.025)/nrow(results_pred)
  #PPOS_trt_mod3 <- sum(results_pred$PP_mod3 > 0.975)/nrow(results_pred)
  PPOS_freq_trt_mod3 <- sum(results_pred$Pval_mod3 < 0.025)/nrow(results_pred)
  
  #consolidate results
  dat_out <- data.frame(
    PP = c(postprob_trt1, postprob_trt2, postprob_trt3),
    #PPOS = c(PPOS_trt_mod1, PPOS_trt_mod2, PPOS_trt_mod3),
    PPOS_freq = c(PPOS_freq_trt_mod1, PPOS_freq_trt_mod2, PPOS_freq_trt_mod3),
    Covariates = c('treatment','treatment,time','treatment,time,time^2'),
    Model = c('Bayesian','Bayesian','Bayesian')
  )
  
  #return PP, pvalue, PPOS and frequentist PPOS for all three models
  return(dat_out)
}

#function to apply bayesian methods to calculate posterior probability given a complete dataset
model_bayes_comp <- function(dat_current, outcome_name){
  #function arguments:
  #dat_current: simulated dataset
  #outcome_name: selected outcome column from simulated dataset.  Must be passed as a string (i.e. 'outcome_linear')
  
  #streamline the dataset
  dat_model <- dat_current %>%
    select(id, outcome_name, randomized, time_accrual) %>%
    mutate(time_accrual_sq = time_accrual^2) %>%
    rename(outcome = outcome_name)
  
  #calculate posterior probability for observed data
  mod1 <- brm(outcome ~ randomized, refresh=0, family='bernoulli', data=dat_model)
  asd1 <- posterior_samples(mod1, pars='b_randomized')
  postprob_trt1 <- mean(asd1$b_randomized > 0) #calculate posterior probability that the log-odds is >0 for treatment
  mod2 <- brm(outcome ~ randomized + time_accrual, refresh=0, family='bernoulli', data=dat_model)
  asd2 <- posterior_samples(mod2, pars='b_randomized')
  postprob_trt2 <- mean(asd2$b_randomized > 0) #calculate posterior probability that the log-odds is >0 for treatment
  mod3 <- brm(outcome ~ randomized + time_accrual + time_accrual_sq, refresh=0, family='bernoulli', data=dat_model)
  asd3 <- posterior_samples(mod3, pars='b_randomized')
  postprob_trt3 <- mean(asd3$b_randomized > 0) #calculate posterior probability that the log-odds is >0 for treatment
  
  #consolidate results
  dat_out <- data.frame(
    PP = c(postprob_trt1, postprob_trt2, postprob_trt3),
    Covariates = c('treatment','treatment,time','treatment,time,time^2'),
    Model = c('Bayesian (Complete Data)','Bayesian (Complete Data)','Bayesian (Complete Data)')
  )
  
  #return PP results
  return(dat_out)
}

###############################################################
### FUNCTIONS FOR RUNNING SIMULATION SCENARIOS:             ###
###############################################################

runSim <- function(seed,n_sim,interim_fraction,T,n,RR,p_start,p_end,accrual_type,outcome_name,sim_name,filename_out){
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
  dat_out <- matrix(nrow = (n_sim*21), ncol = 6) %>% as.data.frame()
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
    
    #run Bayesian models
    mod_Bayes <- model_bayes(dat_current = dat_interim, outcome_name = 'outcome', n_total = n)
    mod_Bayes_comp <- model_bayes_comp(dat_current = dat_full, outcome_name = 'outcome')
    
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
    dat_Bayes <- rbind(
      mod_Bayes %>% 
        mutate(Value = PP,
               ValueName = rep('PP',nrow(mod_Bayes)),
               Simulation = rep(i,nrow(mod_Bayes)),
               Seed = rep(seed,nrow(mod_Bayes))) %>% 
        select(Value,ValueName,Covariates,Model,Simulation,Seed),
      mod_Bayes %>% 
        mutate(Value = PPOS_freq,
               ValueName = rep('PPOS_freq',nrow(mod_Bayes)),
               Simulation = rep(i,nrow(mod_Bayes)),
               Seed = rep(seed,nrow(mod_Bayes))) %>% 
        select(Value,ValueName,Covariates,Model,Simulation,Seed)
    )
    dat_Bayes_comp <- mod_Bayes_comp %>%
      mutate(Value = PP,
             ValueName = rep('PP',nrow(mod_Bayes_comp)),
             Simulation = rep(i,nrow(mod_Bayes_comp)),
             Seed = rep(seed,nrow(mod_Bayes_comp))) %>%
      select(Value,ValueName,Covariates,Model,Simulation,Seed)
    dat_out[(1+(i-1)*21):(i*21),] <- rbind(dat_LR,dat_LR_comp,dat_LR_UP,dat_LR_UP_WT,dat_Bayes,dat_Bayes_comp)
    
    #update output and dataset in the appropriate file directory
    file_out_data <- paste(paste(paste(paste(filename_out,sim_name,sep = '/'),'_Dataset_Simulation',sep = ''),i,sep = ''),'.Rdata',sep = '')
    file_out_results_txt <- paste(paste(paste(paste(filename_out,sim_name,sep = '/'),'_CurrentResults_Simulation',sep = ''),i,sep = ''),'.txt',sep = '')
    save(dat,file = file_out_data)
    write.table(dat_out,file = file_out_results_txt, row.names=FALSE, col.names=FALSE)
  }
}

###############################################################
### file: 00_BasicFunctions                                 ###
### authors: Jess Wild                                      ###
### creation date:    07/18/2023                            ###
### latest edit date: 07/18/2023                            ###
### description: Code defining basic functions for data     ### 
###   generation and model fitting                          ###
###############################################################

###############################################################
### FUNCTIONS TO GENERATE DATA:                             ###
###############################################################

#construct the overall data simulation function
datSim <- function(T,n,RR,p_start,p_end,accrual_type){
  #function arguments:
    #T: duration of study (in days)
    #n: total sample size
    #RR: relative risk of treatment group compared to control group
    #p_start: probability of outcome in the control group (between 0,1) at lower end of probability range
    #p_end: probability of outcome in the control group (between 0,1) at higher end of probability range
    #accrual_type: choice of accrual pattern from (unif,beta_a,beta_b,beta_c) 
          #beta_a: more accrual early in the trial
          #beta_b: more accrual later in the trial
          #beta_c: bimodal shape with two peaks of accrual
  
  #define vectors for subject ids and time interval from study beginning to end
  id <- seq(1:n) #arbitrary patient ids
  time <- seq(1:T) #sequential vector of study days from start to finish (for defining/checking time trends)
  
  #define time trend functions for specified time interval 
  slope <- (p_end-p_start)/(T-1) #m=(y_2-y_1)/(x_2-x_1) line formula solved for (p_start,1),(p_end,T)
  p_ctr_constant_low <- rep(p_start,T) #constant, low probability
  p_trt_constant_low <- RR*p_ctr_constant_low
  p_ctr_constant_high <- rep(p_end,T) #constant, high probability
  p_trt_constant_high <- RR*p_ctr_constant_high
  p_ctr_linear <- p_start + slope*(time-1) #linear
  p_trt_linear <- RR*p_ctr_linear
  p_ctr_linearNEG <- p_end - slope*(time-1) #linear
  p_trt_linearNEG <- RR*p_ctr_linearNEG
  p_ctr_quad <- ((p_start-p_end)/(1-T)^2)*(time-T)^2+p_end #quadratic y=a(x-h)^2+k where (h,k) is the vertex
  p_trt_quad <- RR*p_ctr_quad
  p_ctr_quadNEG <- ((p_end-p_start)/(T-1)^2)*(time-1)^2+p_start #quadratic
  p_trt_quadNEG <- RR*p_ctr_quadNEG
  #plot the trends over the time vector for Treatment & Control
  #par(mfrow=c(2,3))
  #plot(x = time, y = p_ctr_constant_low, col="pink", ylim=c(0,2), main='Constant')
  #lines(x = time, y = p_trt_constant_low, col="magenta")
  #lines(x = time, y = p_trt_constant_high, col="blue")
  #plot(x = time, y = p_ctr_linear, col="purple", main='Linear')
  #lines(x = time, y = p_trt_linear, col="plum")
  #plot(x = time, y = p_ctr_linearNEG, col="purple", main='Linear (decrease)')
  #lines(x = time, y = p_trt_linearNEG, col="plum")
  #plot(x = time, y = p_ctr_quad, col="blue",main='Quadratic')
  #lines(x = time, y = p_trt_quad, col="lightblue")
  #plot(x = time, y = p_ctr_quadNEG, col="blue",main='Quadratic (inverted)')
  #lines(x = time, y = p_trt_quadNEG, col="lightblue")
  
  #randomly simulate accrual times for each patient within the study duration, based on user input choice
  time_accrual_unif <- runif(n=n, min = 1, max = T) %>% #assumes uniform accrual between first and last study days
    floor() %>% #round down to nearest day (fractional days don't make sense here)
    sort() #reorder times from first to last
  time_accrual_beta_a  <- (rbeta(n=n, 1, 10, ncp = 1)*T+1) %>% #beta distributed accrual with more accrual early in the trial
    floor() %>% #round down to nearest day (fractional days don't make sense here)
    sort() #reorder times from first to last
  time_accrual_beta_b  <- (rbeta(n=n, 10, 1, ncp = 1)*T+1) %>% #beta distributed accrual with more accrual later in the trial
    floor() %>% #round down to nearest day (fractional days don't make sense here)
    sort() #reorder times from first to last
  time_accrual_beta_c <- c(rbeta(n=n/2, 1, 10, ncp = 1)*(T/2)+1,rbeta(n=n/2, 1, 10, ncp = 1)*(T/2)+(T/2)+1) %>% #beta distributed accrual with two peaks of accrual
    floor() %>% #round down to nearest day (fractional days don't make sense here)
    sort() #reorder times from first to last
  #specify particular accrual time
  time_accrual <- (accrual_type=='unif')*time_accrual_unif +
    (accrual_type=='beta_a')*time_accrual_beta_a +
    (accrual_type=='beta_b')*time_accrual_beta_b +
    (accrual_type=='beta_c')*time_accrual_beta_c
  #plot accrual pattern over time
  #plot(ecdf(time_accrual), main='Accrual CDF')
  
  
  #permuted block randomization
  rand1 <- rpbrPar(N=n,
                   rb=c(2,4,6,8,10), #potential block sizes
                   K=2, #number of study arms
                   groups = c('Treatment','Control') #labels for the study arms
  )
  rand2 <- genSeq(rand1)
  randomized <- t(getRandList(rand2))
  
  #probability calculated based on accrual time and treatment group
  prob_constant_low <- (randomized=='Control')*p_ctr_constant_low[time_accrual] + 
    (randomized=='Treatment')*p_trt_constant_low[time_accrual]
  prob_constant_high <- (randomized=='Control')*p_ctr_constant_high[time_accrual] + 
    (randomized=='Treatment')*p_trt_constant_high[time_accrual]
  prob_linear <- (randomized=='Control')*p_ctr_linear[time_accrual] + 
    (randomized=='Treatment')*p_trt_linear[time_accrual]
  prob_linearNEG <- (randomized=='Control')*p_ctr_linearNEG [time_accrual] + 
    (randomized=='Treatment')*p_trt_linearNEG [time_accrual]
  prob_quad <- (randomized=='Control')*p_ctr_quad[time_accrual] + 
    (randomized=='Treatment')*p_trt_quad[time_accrual]
  prob_quadNEG <- (randomized=='Control')*p_ctr_quadNEG[time_accrual] + 
    (randomized=='Treatment')*p_trt_quadNEG[time_accrual]
  
  #simulate outcomes with binomial distribution
  outcome_constant_low <- rbinom(n = n, #simulate outcomes for each subject
                             p = prob_constant_low, #use the individual specific probabilities for constant time trend
                             size = 1 #0,1 outcome
  ) %>% as.factor()
  outcome_constant_high <- rbinom(n = n, #simulate outcomes for each subject
                                 p = prob_constant_high, #use the individual specific probabilities for constant time trend
                                 size = 1 #0,1 outcome
  ) %>% as.factor()
  outcome_linear <- rbinom(n = n, #simulate outcomes for each subject
                           p = prob_linear, #use the individual specific probabilities for increasing linear time trend
                           size = 1 #0,1 outcome
  ) %>% as.factor()
  outcome_linearNEG <- rbinom(n = n, #simulate outcomes for each subject
                           p = prob_linearNEG, #use the individual specific probabilities for decreasing linear time trend
                           size = 1 #0,1 outcome
  ) %>% as.factor()
  outcome_quad <- rbinom(n = n, #simulate outcomes for each subject
                         p = prob_quad, #use the individual specific probabilities for quadratic time trend
                         size = 1 #0,1 outcome
  ) %>% as.factor()
  outcome_quadNEG <- rbinom(n = n, #simulate outcomes for each subject
                         p = prob_quadNEG, #use the individual specific probabilities for inverted quadratic time trend
                         size = 1 #0,1 outcome
  ) %>% as.factor()
  
  #combine into single dataset
  dat_sim <- data.frame(id, 
                        time_accrual, 
                        randomized,
                        prob_constant_low,
                        prob_constant_high,
                        prob_linear,
                        prob_linearNEG,
                        prob_quad,
                        prob_quadNEG,
                        outcome_constant_low,
                        outcome_constant_high,
                        outcome_linear,
                        outcome_linearNEG,
                        outcome_quad,
                        outcome_quadNEG
  )
  return(dat_sim)
}

###############################################################
### FUNCTIONS TO PERFORM MODELING:                          ###
###############################################################

#function to apply logistic regression using only current data
model_LR <- function(dat_current, outcome_name){
  #function arguments:
    #dat_current: simulated dataset
    #outcome_name: selected outcome column from simulated dataset.  Must be passed as a string (i.e. 'outcome_linear')
  
  #streamline the dataset
  dat_model <- dat_current %>%
    select(id, outcome_name, randomized, time_accrual) %>%
    mutate(time_accrual_sq = time_accrual^2) %>%
    rename(outcome = outcome_name)
  
  #fit logistic regression models
  mod1 <- glm(outcome ~ randomized, data = dat_model, family = "binomial")
  mod2 <- glm(outcome ~ randomized + time_accrual, data = dat_model, family = "binomial")
  mod3 <- glm(outcome ~ randomized + time_accrual + time_accrual_sq, data = dat_model, family = "binomial")
  
  #report p values
  dat_pval <- data.frame(
    P_treatment = c(coef(summary(mod1))[2,4],coef(summary(mod2))[2,4],coef(summary(mod3))[2,4]),
    P_time = c(NA,coef(summary(mod2))[3,4],coef(summary(mod3))[3,4]),
    P_time_sq = c(NA,NA,coef(summary(mod3))[4,4]),
    Covariates = c('treatment','treatment,time','treatment,time,time_sq'),
    Model = c('LR','LR','LR')
  )
  return(dat_pval)
}

#function to apply logistic regression with upstrapped complete data
model_LR_uppstrap <- function(dat_current, outcome_name, n_total){
  #function arguments:
    #dat_current: simulated dataset
    #outcome_name: selected outcome column from simulated dataset.  Must be passed as a string (i.e. 'outcome_linear')
    #n_total: total desired sample size for uppstrap resampling purposes
  
  #streamline the dataset
  dat_model <- dat_current %>%
    select(id, outcome_name, randomized, time_accrual) %>%
    mutate(time_accrual_sq = time_accrual^2) %>%
    rename(outcome = outcome_name)
  
  #split into treatment and control datasets
  dat_trt <- subset(dat_model, randomized == 'Treatment')
  dat_ctr <- subset(dat_model, randomized == 'Control')
  
  #define number of samples to upstrapp (calculated for treatment and control groups separately)
  n_uppstrap_trt <- n_total/2 - nrow(dat_trt)
  n_uppstrap_ctr <- n_total/2 - nrow(dat_ctr)
  
  #uppstrap resampling of the data (separated for treatment and control groups to preserve allocation ratio)
  dat_uppstrap_trt <- sample_n(dat_trt, size = n_uppstrap_trt, replace = T)
  dat_uppstrap_ctr <- sample_n(dat_ctr, size = n_uppstrap_ctr, replace = T)
  dat_complete <- rbind(dat_model,dat_uppstrap_trt,dat_uppstrap_ctr)
  
  #fit logistic regression models
  mod1 <- glm(outcome ~ randomized, data = dat_complete, family = "binomial")
  mod2 <- glm(outcome ~ randomized + time_accrual, data = dat_complete, family = "binomial")
  mod3 <- glm(outcome ~ randomized + time_accrual + time_accrual_sq, data = dat_complete, family = "binomial")
  
  #report p values 
  dat_pval <- data.frame(
    P_treatment = c(coef(summary(mod1))[2,4],coef(summary(mod2))[2,4],coef(summary(mod3))[2,4]),
    P_time = c(NA,coef(summary(mod2))[3,4],coef(summary(mod3))[3,4]),
    P_time_sq = c(NA,NA,coef(summary(mod3))[4,4]),
    Covariates = c('treatment','treatment,time','treatment,time,time_sq'),
    Model = c('Uppstrapped LR','Uppstrapped LR','Uppstrapped LR')
  )
  return(dat_pval)
}

#function to apply logistic regression with weighted upstrapping used to get complete data
model_LR_uppstrap_wt <- function(dat_current, outcome_name, n_total){
  #function arguments:
  #dat_current: simulated dataset
  #outcome_name: selected outcome column from simulated dataset.  Must be passed as a string (i.e. 'outcome_linear')
  #n_total: total desired sample size for uppstrap resampling purposes
  
  #define current study day and current sample size (needed to calculate weights)
  t_star <- max(dat_current$time_accrual)
  n <- nrow(dat_current)
  
  #streamline the dataset
  dat_model <- dat_current %>%
    select(id, outcome_name, randomized, time_accrual) %>%
    mutate(time_accrual_sq = time_accrual^2) %>%
    mutate(weight = ((1/n)^(time_accrual/t_star))/(sum((1/n)^(time_accrual/t_star)))) %>% #calculate weight for each individual based on accrual time
    rename(outcome = outcome_name)
  
  #split into treatment and control datasets
  dat_trt <- subset(dat_model, randomized == 'Treatment')
  dat_ctr <- subset(dat_model, randomized == 'Control')
  
  #define number of samples to upstrapp (calculated for treatment and control groups separately)
  n_uppstrap_trt <- n_total/2 - nrow(dat_trt)
  n_uppstrap_ctr <- n_total/2 - nrow(dat_ctr)
  
  #uppstrap resampling of the data (separated for treatment and control groups to preserve allocation ratio)
  dat_uppstrap_trt <- sample_n(dat_trt, size = n_uppstrap_trt, replace = T, weight = dat_trt$weight)
  dat_uppstrap_ctr <- sample_n(dat_ctr, size = n_uppstrap_ctr, replace = T, weight = dat_ctr$weight)
  dat_complete <- rbind(dat_model,dat_uppstrap_trt,dat_uppstrap_ctr)
  
  #fit logistic regression models
  mod1 <- glm(outcome ~ randomized, data = dat_complete, family = "binomial")
  mod2 <- glm(outcome ~ randomized + time_accrual, data = dat_complete, family = "binomial")
  mod3 <- glm(outcome ~ randomized + time_accrual + time_accrual_sq, data = dat_complete, family = "binomial")
  
  #report p values
  dat_pval <- data.frame(
    P_treatment = c(coef(summary(mod1))[2,4],coef(summary(mod2))[2,4],coef(summary(mod3))[2,4]),
    P_time = c(NA,coef(summary(mod2))[3,4],coef(summary(mod3))[3,4]),
    P_time_sq = c(NA,NA,coef(summary(mod3))[4,4]),
    Covariates = c('treatment','treatment,time','treatment,time,time_sq'),
    Model = c('Weighted Uppstrapped LR','Weighted Uppstrapped LR','Weighted Uppstrapped LR')
  )
  return(dat_pval)
}

#function to apply logistic regression using only current data: only basic model (no time covariates)
model_LR_simp <- function(dat_current, outcome_name){
  #function arguments:
  #dat_current: simulated dataset
  #outcome_name: selected outcome column from simulated dataset.  Must be passed as a string (i.e. 'outcome_linear')
  
  #streamline the dataset
  dat_model <- dat_current %>%
    select(id, outcome_name, randomized, time_accrual) %>%
    mutate(time_accrual_sq = time_accrual^2) %>%
    rename(outcome = outcome_name)
  
  #fit logistic regression models
  mod1 <- glm(outcome ~ randomized, data = dat_model, family = "binomial")
  
  #report p values
  dat_pval <- data.frame(
    P_treatment = c(coef(summary(mod1))[2,4])
  )
  return(dat_pval)
}

#function to apply logistic regression with upstrapped complete data: only basic model (no time covariates)
model_LR_uppstrap_simp <- function(dat_current, outcome_name, n_total){
  #function arguments:
  #dat_current: simulated dataset
  #outcome_name: selected outcome column from simulated dataset.  Must be passed as a string (i.e. 'outcome_linear')
  #n_total: total desired sample size for uppstrap resampling purposes
  
  #streamline the dataset
  dat_model <- dat_current %>%
    select(id, outcome_name, randomized, time_accrual) %>%
    mutate(time_accrual_sq = time_accrual^2) %>%
    rename(outcome = outcome_name)
  
  #split into treatment and control datasets
  dat_trt <- subset(dat_model, randomized == 'Treatment')
  dat_ctr <- subset(dat_model, randomized == 'Control')
  
  #define number of samples to upstrapp (calculated for treatment and control groups separately)
  n_uppstrap_trt <- n_total/2 - nrow(dat_trt)
  n_uppstrap_ctr <- n_total/2 - nrow(dat_ctr)
  
  #uppstrap resampling of the data (separated for treatment and control groups to preserve allocation ratio)
  dat_uppstrap_trt <- sample_n(dat_trt, size = n_uppstrap_trt, replace = T)
  dat_uppstrap_ctr <- sample_n(dat_ctr, size = n_uppstrap_ctr, replace = T)
  dat_complete <- rbind(dat_model,dat_uppstrap_trt,dat_uppstrap_ctr)
  
  #fit logistic regression models
  mod1 <- glm(outcome ~ randomized, data = dat_complete, family = "binomial")
  
  #report p values 
  dat_pval <- data.frame(
    P_treatment = c(coef(summary(mod1))[2,4])
  )
  return(dat_pval)
}

#function to apply chi-squared OR fisher's exact test (chi-squared preffered when useable): only basic model (no time covariates)
model_CSFE_simp <- function(dat_current, outcome_name){
    #function arguments:
    #dat_current: simulated dataset
    #outcome_name: selected outcome column from simulated dataset.  Must be passed as a string (i.e. 'outcome_linear')
    
    #streamline the dataset and create contingengy table
    dat_model <- dat_current %>%
      select(outcome_name, randomized) %>%
      rename(outcome = outcome_name) %>%
      ftable()
    
    #fit chi squared OR fishers exact test
    test <- 'CS'
    catk <- chisq.test(dat_model) #CS if possible
    if( sum(catk$expected < 5) > 0 ){ test <- 'FE' }
    if( sum(catk$expected < 5) > 0 ){ catk <- fisher.test(dat_model) } #FE if needed
    
    #report p values
    dat_pval <- data.frame(
      P_treatment = catk$p.value,
      TestUsed = test
    )
    return(dat_pval)
}

#function to apply chi-squared OR fisher's exact test (chi-squared preferred when usable) with upstrapped complete data: only basic model (no time covariates)
model_CSFE_uppstrap_simp <- function(dat_current, outcome_name, n_total){
  #function arguments:
  #dat_current: simulated dataset
  #outcome_name: selected outcome column from simulated dataset.  Must be passed as a string (i.e. 'outcome_linear')
  #n_total: total desired sample size for uppstrap resampling purposes
  
  #streamline the dataset
  dat_model <- dat_current %>%
    select(id, outcome_name, randomized, time_accrual) %>%
    mutate(time_accrual_sq = time_accrual^2) %>%
    rename(outcome = outcome_name)
  
  #split into treatment and control datasets
  dat_trt <- subset(dat_model, randomized == 'Treatment')
  dat_ctr <- subset(dat_model, randomized == 'Control')
  
  #define number of samples to upstrapp (calculated for treatment and control groups separately)
  n_uppstrap_trt <- n_total/2 - nrow(dat_trt)
  n_uppstrap_ctr <- n_total/2 - nrow(dat_ctr)
  
  #uppstrap resampling of the data (separated for treatment and control groups to preserve allocation ratio)
  dat_uppstrap_trt <- sample_n(dat_trt, size = n_uppstrap_trt, replace = T)
  dat_uppstrap_ctr <- sample_n(dat_ctr, size = n_uppstrap_ctr, replace = T)
  dat_complete <- rbind(dat_model,dat_uppstrap_trt,dat_uppstrap_ctr)
  
  #create contingency table
  dat_tab <- dat_complete %>%
    select(outcome, randomized) %>%
    ftable()
  
  #fit chi squared OR fishers exact test
  test <- 'CS'
  catk <- chisq.test(dat_tab) #CS if possible
  if( sum(catk$expected < 5) > 0 ){ test <- 'FE' }
  if( sum(catk$expected < 5) > 0 ){ catk <- fisher.test(dat_tab) } #FE if needed
  
  #report p values
  dat_pval <- data.frame(
    P_treatment = catk$p.value,
    TestUsed = test
  )
  return(dat_pval)
}

# Simulate example data set
set.seed(030322)
trt <- c( rep(0,250), rep(1,250) )
t <- runif(n=500)
xb <- -5 + 1*trt + 0*t
p <- 1/(1 + exp(-xb))
y <- rbinom(n=500, size=1, prob=p)
dat <- data.frame( y=y, trt=trt, t=t )

####################################
### Using brms

library(brms)

mod1 <- brm(y ~ trt + t, refresh=0, family='bernoulli', data=dat )

# calculate posterior probability for observed data
asd <- posterior_samples(mod1, pars='b_trt')
postprob_trt <- mean(asd$b_trt > 0) # calculate posterior probability that the log-odds is >0 for treatment

# create data frame with remaining participants to enroll
pred_dat <- data.frame(trt=c(rep(0,50),rep(1,50)), t=rep(seq(1,2,length.out=50), 2) )  

# generate remaining data based on the predictive posterior distribution
# NOTE: produces results that are 1000 (nsamples) x 100 (the number of overall new observations)
mod1Pred <- posterior_predict(newdata=pred_dat, object=mod1, allow_new_levels=TRUE, re_formula=NULL, nsamples=1000)

# need to calculate the posterior probability for each of the 1000 datasets where we add the predictive posterior generated observations
# NOTE: might as well also calculate a frequentist p-value since we could compare this to an upstrap as a way of generating future data with a Bayesian approach but frequentist analysis...might be analogous to conditional power but with a Bayesian twist
# NOTE: probably easiest to wrap this in a for loop from 1:1000 and save the posterior probability (PP) and p-value for treatment

dat_pp <- rbind( dat, cbind(y=mod1Pred[1,], pred_dat))

mod1_final <- brm(y ~ trt + t, refresh=0, family='bernoulli', data=dat_pp)
mod1_final_freq <- glm(y ~ trt + t, family='binomial', data=dat_pp)

asd_final <- posterior_samples(mod1_final, pars='b_trt')
postprob_trt_final <- mean(asd_final$b_trt > 0) # calculate posterior probability that the log-odds is >0 for treatment

## NOTES ON POSTERIOR PROBABILITY AND P-VALUE THRESHOLDS:
# If we are considering P(OR>1) as a one-sided test, we can probably use PP>0.975 and p<0.025 as our thresholds
# For the predictive probability of success we'd ultimately summarize the proportion of the 1000 simulated future observation sets where PP>0.975
# For the frequentist twist on the approach we'd similarly summarize the proportion of the 1000 simulated future observation sets where p<0.025

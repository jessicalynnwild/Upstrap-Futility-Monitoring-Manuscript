###############################################################
### file: 01_PrimaryDataGeneration_HeatmapData              ###
### authors: Jess Wild                                      ###
### creation date:    07/18/2023                            ###
### latest edit date: 07/18/2023                            ###
### description: code to generate heatmap style dataset     ###
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
load('./01 Data/Data Processed/dat_results.RData')

###############################################################
### PROCESS RESULTS DATA INTO HEATMAP FORMAT                ###
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
  
  #create progress bar to measure completion of main loop (per simulation setting)
  pb = txtProgressBar(min = 0, max = nrow(dat_heat), initial = 0, style = 3)
  
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
    
    #update progress bar
    setTxtProgressBar(pb,index)
  }
  #progress update
  marker <- paste('Finished Setting: ',setting, sep = '')
  print(marker)
}

###############################################################
### SAVE OUTPUT DATASETS                                    ###
###############################################################

save(dat_heat, file = './01 Data/Data Processed/dat_heat.RData')
#load('./01 Data/Data Processed/dat_heat.RData')

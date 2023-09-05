###############################################################
### file: 04_PrimaryDataAnalysis_Radarplot                  ###
### authors: Jess Wild                                      ###
### creation date:    07/31/2022                            ###
### latest edit date: 08/01/2023                            ###
### description: code to generate radar plot visualizations ###
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
library(plotly)
remotes::install_github("ricardo-bion/ggradar")

#load input datasets
load('./01 Data/Data Processed/dat_final_output_fig.RData') #change to generate from dat_final_output_fig

###############################################################
### Radar Plots: DATA PROCESSING                            ###
###############################################################

#modify raw dataset
dat_radar <- dat_final_output_fig %>% 
  ungroup() %>%
  filter(Method %in% c('AU','CU','GU','OBF','PO','CP 1%','CP 5%','CP 10%','CP 20%')) %>%
  filter(power %in% c(0.05,0.80)) %>%
  select(n,power,ExpectedN_FO_mean,DecisionFO,RejectionRate_FO_FS,Method) %>%
  group_by(n) %>%
  pivot_wider(
    id_cols = c('n','Method'),
    names_from = c('power'),
    values_from = c('ExpectedN_FO_mean','DecisionFO','RejectionRate_FO_FS')
  ) %>%
  mutate(ExpectedN_FO_mean_0.05 = ExpectedN_FO_mean_0.05/n,
         ExpectedN_FO_mean_0.8 = ExpectedN_FO_mean_0.8/n,
         Method = as.factor(Method)) %>%
  select(n,Method,ExpectedN_FO_mean_0.05,ExpectedN_FO_mean_0.8,DecisionFO_0.05,DecisionFO_0.8,RejectionRate_FO_FS_0.05,
         RejectionRate_FO_FS_0.8) %>%
  rename(ESS_H0 = ExpectedN_FO_mean_0.05,ESS_H1 = ExpectedN_FO_mean_0.8,Prop_H0 = DecisionFO_0.05,
         Prop_H1 = DecisionFO_0.8,TIE = RejectionRate_FO_FS_0.05,Power = RejectionRate_FO_FS_0.8) %>%
  ungroup() %>%
  mutate(ESS_H0 = as.numeric(ESS_H0),ESS_H1 = as.numeric(ESS_H1),Prop_H0 = as.numeric(Prop_H0),Prop_H1 = as.numeric(Prop_H1),Power = as.numeric(Power),TIE = as.numeric(TIE), n = as.character(n)) %>%
  rename(group = Method,`ESS (H0)` = ESS_H0, `ESS (H1)` = ESS_H1, `Proportion (H0)` = Prop_H0, `Proportion (H1)` = Prop_H1)

#split by sample size
dat_radar_40 <- dat_radar %>% filter(n == 40) %>% select(-n) 
dat_radar_160 <- dat_radar %>% filter(n == 160) %>% select(-n)
dat_radar_600 <- dat_radar %>% filter(n == 600) %>% select(-n)
dat_radar_2000 <- dat_radar %>% filter(n == 2000) %>% select(-n)

###############################################################
### Radar Plots: COMBINED PLOT                              ###
###############################################################

ggradar::ggradar(dat_radar) +
  facet_wrap(~ n, ncol = 2)

#######################PLOTLY VERSION
#n = 40
fig_40 <- plot_ly(
  type = 'scatterpolar',
  fill = 'toself'
) 
for (i in 1:nrow(dat_radar_40)) {
  fig_40 <- fig_40 %>%
    add_trace(
      r = dat_radar_40[i,2:ncol(dat_radar_40)] %>% as.numeric(),
      theta = colnames(dat_radar_40)[2:ncol(dat_radar_40)],
      name = dat_radar_40$group[i]
    ) 
}
fig_40 <- fig_40 %>%
  layout(
    polar = list(
      radialaxis = list(
        visible = T,
        range = c(0,1)
      )
    )
  )

#n = 160
fig_160 <- plot_ly(
  type = 'scatterpolar',
  fill = 'toself'
) 
for (i in 1:nrow(dat_radar_160)) {
  fig_160 <- fig_160 %>%
    add_trace(
      r = dat_radar_160[i,3:ncol(dat_radar_160)] %>% as.numeric(),
      theta = colnames(dat_radar_160)[3:ncol(dat_radar_160)],
      name = dat_radar_160$Method[i]
    ) 
}
fig_160 <- fig_160 %>%
  layout(
    polar = list(
      radialaxis = list(
        visible = T,
        range = c(0,1)
      )
    )
  )

#n = 600
fig_600 <- plot_ly(
  type = 'scatterpolar',
  fill = 'toself'
) 
for (i in 1:nrow(dat_radar_600)) {
  fig_600 <- fig_600 %>%
    add_trace(
      r = dat_radar_600[i,3:ncol(dat_radar_600)] %>% as.numeric(),
      theta = colnames(dat_radar_600)[3:ncol(dat_radar_600)],
      name = dat_radar_600$Method[i]
    ) 
}
fig_600 <- fig_600 %>%
  layout(
    polar = list(
      radialaxis = list(
        visible = T,
        range = c(0,1)
      )
    )
  )

#n = 2000

fig_2000 <- plot_ly(
  type = 'scatterpolar',
  fill = 'toself'
) 
for (i in 1:nrow(dat_radar_2000)) {
  fig_2000 <- fig_2000 %>%
    add_trace(
      r = dat_radar_2000[i,3:ncol(dat_radar_2000)] %>% as.numeric(),
      theta = colnames(dat_radar_2000)[3:ncol(dat_radar_2000)],
      name = dat_radar_2000$Method[i]
    ) 
}
fig_2000 <- fig_2000 %>%
  layout(
    polar = list(
      radialaxis = list(
        visible = T,
        range = c(0,1)
      )
    )
  )

#facet wrap plots
title <- list(
  font = list(size = 18),
  xref = 'paper',
  yref = 'paper',
  yanchor = 'bottom',
  xanchor = 'center',
  align = 'center',
  x = .5,
  y = 1,
  showarrow = FALSE)
subplot(b1,b2, shareY = T, titleX = T, titleY = T, 
        nrows = 1)
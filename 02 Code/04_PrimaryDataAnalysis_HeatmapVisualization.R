###############################################################
### file: 05_PrimaryDataAnalysis_HeatmapVisualization       ###
### authors: Jess Wild                                      ###
### creation date:    07/18/2023                            ###
### latest edit date: 07/18/2023                            ###
### description: code to visualize heatmap dataset          ###
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
library(lattice) #for heatmaps
library(latticeExtra) #for heatmaps
library(gridExtra) #for heatmaps

#load input datasets
load('./01 Data/Data Processed/dat_heat.RData')

###############################################################
### VISUALIZE HEATMAP DATA                                  ###
###############################################################

#generate heatmaps for each simulation setting
plot_S1 <- levelplot(Setting1 ~ Proportion*P, data=dat_heat  ,xlab="Proportion", ylab='P-Value',
                     main="Simulation Setting: 1",panel = panel.levelplot.points, cex = 0) + 
  layer_(panel.2dsmoother(..., n = 200))
plot_S2 <- levelplot(Setting2 ~ Proportion*P, data=dat_heat  ,xlab="Proportion", ylab='P-Value',
                     main="Simulation Setting: 2",panel = panel.levelplot.points, cex = 0) + 
  layer_(panel.2dsmoother(..., n = 200))
plot_S3 <- levelplot(Setting3 ~ Proportion*P, data=dat_heat  ,xlab="Proportion", ylab='P-Value',
                     main="Simulation Setting: 3",panel = panel.levelplot.points, cex = 0) + 
  layer_(panel.2dsmoother(..., n = 200))
plot_S4 <- levelplot(Setting4 ~ Proportion*P, data=dat_heat  ,xlab="Proportion", ylab='P-Value',
                     main="Simulation Setting: 4",panel = panel.levelplot.points, cex = 0) + 
  layer_(panel.2dsmoother(..., n = 200))
plot_S5 <- levelplot(Setting5 ~ Proportion*P, data=dat_heat  ,xlab="Proportion", ylab='P-Value',
                     main="Simulation Setting: 5",panel = panel.levelplot.points, cex = 0) + 
  layer_(panel.2dsmoother(..., n = 200))
plot_S6 <- levelplot(Setting6 ~ Proportion*P, data=dat_heat  ,xlab="Proportion", ylab='P-Value',
                     main="Simulation Setting: 6",panel = panel.levelplot.points, cex = 0) + 
  layer_(panel.2dsmoother(..., n = 200))
plot_S7 <- levelplot(Setting7 ~ Proportion*P, data=dat_heat  ,xlab="Proportion", ylab='P-Value',
                     main="Simulation Setting: 7",panel = panel.levelplot.points, cex = 0) + 
  layer_(panel.2dsmoother(..., n = 200))
plot_S8 <- levelplot(Setting8 ~ Proportion*P, data=dat_heat  ,xlab="Proportion", ylab='P-Value',
                     main="Simulation Setting: 8",panel = panel.levelplot.points, cex = 0) + 
  layer_(panel.2dsmoother(..., n = 200))
plot_S9 <- levelplot(Setting9 ~ Proportion*P, data=dat_heat  ,xlab="Proportion", ylab='P-Value',
                     main="Simulation Setting: 9",panel = panel.levelplot.points, cex = 0) + 
  layer_(panel.2dsmoother(..., n = 200))
plot_S10 <- levelplot(Setting10 ~ Proportion*P, data=dat_heat  ,xlab="Proportion", ylab='P-Value',
                      main="Simulation Setting: 10",panel = panel.levelplot.points, cex = 0) + 
  layer_(panel.2dsmoother(..., n = 200))
plot_S11 <- levelplot(Setting11 ~ Proportion*P, data=dat_heat  ,xlab="Proportion", ylab='P-Value',
                      main="Simulation Setting: 11",panel = panel.levelplot.points, cex = 0) + 
  layer_(panel.2dsmoother(..., n = 200))
plot_S12 <- levelplot(Setting12 ~ Proportion*P, data=dat_heat  ,xlab="Proportion", ylab='P-Value',
                      main="Simulation Setting: 12",panel = panel.levelplot.points, cex = 0) + 
  layer_(panel.2dsmoother(..., n = 200))
plot_S13 <- levelplot(Setting13 ~ Proportion*P, data=dat_heat  ,xlab="Proportion", ylab='P-Value',
                      main="Simulation Setting: 13",panel = panel.levelplot.points, cex = 0) + 
  layer_(panel.2dsmoother(..., n = 200))
plot_S14 <- levelplot(Setting14 ~ Proportion*P, data=dat_heat  ,xlab="Proportion", ylab='P-Value',
                      main="Simulation Setting: 14",panel = panel.levelplot.points, cex = 0) + 
  layer_(panel.2dsmoother(..., n = 200))
plot_S15 <- levelplot(Setting15 ~ Proportion*P, data=dat_heat  ,xlab="Proportion", ylab='P-Value',
                      main="Simulation Setting: 15",panel = panel.levelplot.points, cex = 0) + 
  layer_(panel.2dsmoother(..., n = 200))
plot_S16 <- levelplot(Setting16 ~ Proportion*P, data=dat_heat  ,xlab="Proportion", ylab='P-Value',
                      main="Simulation Setting: 16",panel = panel.levelplot.points, cex = 0) + 
  layer_(panel.2dsmoother(..., n = 200))
plot_S17 <- levelplot(Setting17 ~ Proportion*P, data=dat_heat  ,xlab="Proportion", ylab='P-Value',
                      main="Simulation Setting: 17",panel = panel.levelplot.points, cex = 0) + 
  layer_(panel.2dsmoother(..., n = 200))
plot_S18 <- levelplot(Setting18 ~ Proportion*P, data=dat_heat  ,xlab="Proportion", ylab='P-Value',
                      main="Simulation Setting: 18",panel = panel.levelplot.points, cex = 0) + 
  layer_(panel.2dsmoother(..., n = 200))
plot_S19 <- levelplot(Setting19 ~ Proportion*P, data=dat_heat  ,xlab="Proportion", ylab='P-Value',
                      main="Simulation Setting: 19",panel = panel.levelplot.points, cex = 0) + 
  layer_(panel.2dsmoother(..., n = 200))
plot_S20 <- levelplot(Setting20 ~ Proportion*P, data=dat_heat  ,xlab="Proportion", ylab='P-Value',
                      main="Simulation Setting: 20",panel = panel.levelplot.points, cex = 0) + 
  layer_(panel.2dsmoother(..., n = 200))
plot_S21 <- levelplot(Setting21 ~ Proportion*P, data=dat_heat  ,xlab="Proportion", ylab='P-Value',
                      main="Simulation Setting: 21",panel = panel.levelplot.points, cex = 0) + 
  layer_(panel.2dsmoother(..., n = 200))
plot_S22 <- levelplot(Setting22 ~ Proportion*P, data=dat_heat  ,xlab="Proportion", ylab='P-Value',
                      main="Simulation Setting: 22",panel = panel.levelplot.points, cex = 0) + 
  layer_(panel.2dsmoother(..., n = 200))
plot_S23 <- levelplot(Setting23 ~ Proportion*P, data=dat_heat  ,xlab="Proportion", ylab='P-Value',
                      main="Simulation Setting: 23",panel = panel.levelplot.points, cex = 0) + 
  layer_(panel.2dsmoother(..., n = 200))
plot_S24 <- levelplot(Setting24 ~ Proportion*P, data=dat_heat  ,xlab="Proportion", ylab='P-Value',
                      main="Simulation Setting: 24",panel = panel.levelplot.points, cex = 0) + 
  layer_(panel.2dsmoother(..., n = 200))
plot_S25 <- levelplot(Setting25 ~ Proportion*P, data=dat_heat  ,xlab="Proportion", ylab='P-Value',
                      main="Simulation Setting: 25",panel = panel.levelplot.points, cex = 0) + 
  layer_(panel.2dsmoother(..., n = 200))
plot_S26 <- levelplot(Setting26 ~ Proportion*P, data=dat_heat  ,xlab="Proportion", ylab='P-Value',
                      main="Simulation Setting: 26",panel = panel.levelplot.points, cex = 0) + 
  layer_(panel.2dsmoother(..., n = 200))
plot_S27 <- levelplot(Setting27 ~ Proportion*P, data=dat_heat  ,xlab="Proportion", ylab='P-Value',
                      main="Simulation Setting: 27",panel = panel.levelplot.points, cex = 0) + 
  layer_(panel.2dsmoother(..., n = 200))
plot_S28 <- levelplot(Setting28 ~ Proportion*P, data=dat_heat  ,xlab="Proportion", ylab='P-Value',
                      main="Simulation Setting: 28",panel = panel.levelplot.points, cex = 0) + 
  layer_(panel.2dsmoother(..., n = 200))
plot_S29 <- levelplot(Setting29 ~ Proportion*P, data=dat_heat  ,xlab="Proportion", ylab='P-Value',
                      main="Simulation Setting: 29",panel = panel.levelplot.points, cex = 0) + 
  layer_(panel.2dsmoother(..., n = 200))
plot_S30 <- levelplot(Setting30 ~ Proportion*P, data=dat_heat  ,xlab="Proportion", ylab='P-Value',
                      main="Simulation Setting: 30",panel = panel.levelplot.points, cex = 0) + 
  layer_(panel.2dsmoother(..., n = 200))
plot_S31 <- levelplot(Setting31 ~ Proportion*P, data=dat_heat  ,xlab="Proportion", ylab='P-Value',
                      main="Simulation Setting: 31",panel = panel.levelplot.points, cex = 0) + 
  layer_(panel.2dsmoother(..., n = 200))
plot_S32 <- levelplot(Setting32 ~ Proportion*P, data=dat_heat  ,xlab="Proportion", ylab='P-Value',
                      main="Simulation Setting: 32",panel = panel.levelplot.points, cex = 0) + 
  layer_(panel.2dsmoother(..., n = 200))
plot_S33 <- levelplot(Setting33 ~ Proportion*P, data=dat_heat  ,xlab="Proportion", ylab='P-Value',
                      main="Simulation Setting: 33",panel = panel.levelplot.points, cex = 0) + 
  layer_(panel.2dsmoother(..., n = 200))
plot_S34 <- levelplot(Setting34 ~ Proportion*P, data=dat_heat  ,xlab="Proportion", ylab='P-Value',
                      main="Simulation Setting: 34",panel = panel.levelplot.points, cex = 0) + 
  layer_(panel.2dsmoother(..., n = 200))
plot_S35 <- levelplot(Setting35 ~ Proportion*P, data=dat_heat  ,xlab="Proportion", ylab='P-Value',
                      main="Simulation Setting: 35",panel = panel.levelplot.points, cex = 0) + 
  layer_(panel.2dsmoother(..., n = 200))
plot_S36 <- levelplot(Setting36 ~ Proportion*P, data=dat_heat  ,xlab="Proportion", ylab='P-Value',
                      main="Simulation Setting: 36",panel = panel.levelplot.points, cex = 0) + 
  layer_(panel.2dsmoother(..., n = 200))
plot_S37 <- levelplot(Setting37 ~ Proportion*P, data=dat_heat  ,xlab="Proportion", ylab='P-Value',
                      main="Simulation Setting: 37",panel = panel.levelplot.points, cex = 0) + 
  layer_(panel.2dsmoother(..., n = 200))
plot_S38 <- levelplot(Setting38 ~ Proportion*P, data=dat_heat  ,xlab="Proportion", ylab='P-Value',
                      main="Simulation Setting: 38",panel = panel.levelplot.points, cex = 0) + 
  layer_(panel.2dsmoother(..., n = 200))
plot_S39 <- levelplot(Setting39 ~ Proportion*P, data=dat_heat  ,xlab="Proportion", ylab='P-Value',
                      main="Simulation Setting: 39",panel = panel.levelplot.points, cex = 0) + 
  layer_(panel.2dsmoother(..., n = 200))
plot_S40 <- levelplot(Setting40 ~ Proportion*P, data=dat_heat  ,xlab="Proportion", ylab='P-Value',
                      main="Simulation Setting: 40",panel = panel.levelplot.points, cex = 0) + 
  layer_(panel.2dsmoother(..., n = 200))
plot_S41 <- levelplot(Setting41 ~ Proportion*P, data=dat_heat  ,xlab="Proportion", ylab='P-Value',
                      main="Simulation Setting: 41",panel = panel.levelplot.points, cex = 0) + 
  layer_(panel.2dsmoother(..., n = 200))
plot_S42 <- levelplot(Setting42 ~ Proportion*P, data=dat_heat  ,xlab="Proportion", ylab='P-Value',
                      main="Simulation Setting: 42",panel = panel.levelplot.points, cex = 0) + 
  layer_(panel.2dsmoother(..., n = 200))
plot_S43 <- levelplot(Setting43 ~ Proportion*P, data=dat_heat  ,xlab="Proportion", ylab='P-Value',
                      main="Simulation Setting: 43",panel = panel.levelplot.points, cex = 0) + 
  layer_(panel.2dsmoother(..., n = 200))
plot_S44 <- levelplot(Setting44 ~ Proportion*P, data=dat_heat  ,xlab="Proportion", ylab='P-Value',
                      main="Simulation Setting: 44",panel = panel.levelplot.points, cex = 0) + 
  layer_(panel.2dsmoother(..., n = 200))
plot_S45 <- levelplot(Setting45 ~ Proportion*P, data=dat_heat  ,xlab="Proportion", ylab='P-Value',
                      main="Simulation Setting: 45",panel = panel.levelplot.points, cex = 0) + 
  layer_(panel.2dsmoother(..., n = 200))
plot_S46 <- levelplot(Setting46 ~ Proportion*P, data=dat_heat  ,xlab="Proportion", ylab='P-Value',
                      main="Simulation Setting: 46",panel = panel.levelplot.points, cex = 0) + 
  layer_(panel.2dsmoother(..., n = 200))
plot_S47 <- levelplot(Setting47 ~ Proportion*P, data=dat_heat  ,xlab="Proportion", ylab='P-Value',
                      main="Simulation Setting: 47",panel = panel.levelplot.points, cex = 0) + 
  layer_(panel.2dsmoother(..., n = 200))
plot_S48 <- levelplot(Setting48 ~ Proportion*P, data=dat_heat  ,xlab="Proportion", ylab='P-Value',
                      main="Simulation Setting: 48",panel = panel.levelplot.points, cex = 0) + 
  layer_(panel.2dsmoother(..., n = 200))

#grid figures: 5% type I error (null scenario)
grid.arrange(plot_S4,plot_S20,plot_S36,
             plot_S3,plot_S19,plot_S35,
             plot_S2,plot_S18,plot_S34,
             plot_S1,plot_S17,plot_S33,
             ncol=3)

#grid figures: 50% power (alternative scenario)
grid.arrange(plot_S8,plot_S24,plot_S40,
             plot_S7,plot_S23,plot_S39,
             plot_S6,plot_S22,plot_S38,
             plot_S5,plot_S21,plot_S37,
             ncol=3)

#grid figures: 80% power (alternative scenario)
grid.arrange(plot_S12,plot_S28,plot_S44,
             plot_S11,plot_S27,plot_S43,
             plot_S10,plot_S26,plot_S42,
             plot_S9,plot_S25,plot_S41,
             ncol=3)

#grid figures: 95% power (alternative scenario)
grid.arrange(plot_S16,plot_S32,plot_S48,
             plot_S15,plot_S31,plot_S47,
             plot_S14,plot_S30,plot_S46,
             plot_S13,plot_S29,plot_S45,
             ncol=3)

#additional plotting options for reference
plot_S6a <- levelplot(Setting6 ~ Proportion*P, data=dat_heat  ,xlab="Proportion", ylab='P-Value',
                      main="Simulation Setting: 6",panel = panel.levelplot.points, cex = 0.5) + 
  layer_(panel.2dsmoother(..., n = 200))
plot_S6b <- levelplot(Setting6 ~ Proportion*P, data=dat_heat  ,xlab="Proportion", ylab='P-Value',
                      main="Simulation Setting: 6",panel = panel.levelplot.points, cex = 0) + 
  layer_(panel.2dsmoother(..., n = 200))
plot_S6c <- levelplot(Setting6 ~ Proportion*P, data=dat_heat  ,xlab="Proportion", ylab='P-Value',
                      main="Simulation Setting: 6") + 
  layer_(panel.2dsmoother(..., n = 200))
grid.arrange(plot_S6a,plot_S6b,plot_S6c,ncol=3)

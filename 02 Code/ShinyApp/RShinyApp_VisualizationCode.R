###############################################################
### file: RShinyApp_VisualizationCode                       ###
### authors: Jess Wild                                      ###
### creation date: 04/28/2023                               ###
### latest edit date:05/01/2023                             ###
### description: code to implement all shiny app visuals    ###
###############################################################

###############################################################
### SETUP WORKSPACE                                         ###
###############################################################

#clear workspace
remove(list = ls())

#set working directory
setwd('C:/Users/wildje/Desktop/Dissertation/Paper 1')

#import packages
library(tidyverse)
library(plotly)
library(fmsb)
library(RColorBrewer)

#load raw dataset
load('./01 Data/Upstrap 10-3-2022 Results/ConsolidatedFinalOutputData.RData')

#wide to long format NOTE: MOVE TO DATA CONSOLIDATION SCRIPT
dat_final_output <- dat_final_output %>% 
  pivot_longer(cols = c('ExpectedN_FO_mean','ExpectedN_EO_mean','ExpectedN_FE_mean'),
               names_to='monitor_type',values_to='ExpectedN_mean') %>%
  mutate(monitor_type = recode_factor(monitor_type, ExpectedN_FO_mean='FO',
                                      ExpectedN_EO_mean='EO',ExpectedN_FE_mean='FE')) %>% 
  pivot_longer(cols = c('ExpectedN_FO_sd','ExpectedN_EO_sd','ExpectedN_FE_sd'),
               names_to='monitor_type2',values_to='ExpectedN_sd') %>%
  mutate(monitor_type2 = recode_factor(monitor_type2, ExpectedN_FO_sd='FO',
                                       ExpectedN_EO_sd='EO',ExpectedN_FE_sd='FE')) %>%
  filter(monitor_type == monitor_type2) %>% 
  pivot_longer(cols = c('DecisionFO','DecisionEO','DecisionFE'),
               names_to='monitor_type3',values_to='Decision') %>%
  mutate(monitor_type3 = recode_factor(monitor_type3, DecisionFO='FO',
                                       DecisionEO='EO',DecisionFE='FE'))%>%
  filter(monitor_type == monitor_type3) %>% 
  pivot_longer(cols = c('DecisionFO_25','DecisionEO_25','DecisionFE_25'),
               names_to='monitor_type4',values_to='Decision_25') %>%
  mutate(monitor_type4 = recode_factor(monitor_type4, DecisionFO_25='FO',
                                       DecisionEO_25='EO',DecisionFE_25='FE'))%>%
  filter(monitor_type == monitor_type4) %>% 
  pivot_longer(cols = c('DecisionFO_50','DecisionEO_50','DecisionFE_50'),
               names_to='monitor_type5',values_to='Decision_50') %>%
  mutate(monitor_type5 = recode_factor(monitor_type5, DecisionFO_50='FO',
                                       DecisionEO_50='EO',DecisionFE_50='FE'))%>%
  filter(monitor_type == monitor_type5) %>%
  pivot_longer(cols = c('DecisionFO_75','DecisionEO_75','DecisionFE_75'),
               names_to='monitor_type6',values_to='Decision_75') %>%
  mutate(monitor_type6 = recode_factor(monitor_type6, DecisionFO_75='FO',
                                       DecisionEO_75='EO',DecisionFE_75='FE'))%>%
  filter(monitor_type == monitor_type6) %>%
  pivot_longer(cols = c('RejectionRate_FO_FS','RejectionRate_EO_FS','RejectionRate_FE_FS'),
               names_to='monitor_type7',values_to='RejectionRate') %>%
  mutate(monitor_type7 = recode_factor(monitor_type7, RejectionRate_FO_FS='FO',
                                       RejectionRate_EO_FS='EO',RejectionRate_FE_FS='FE'))%>%
  filter(monitor_type == monitor_type7) %>%
  mutate(method = recode_factor(Method,`Arbitrary Upstrap`='AU',
                                `Arbitrary Upstrap (0.25 Stopping Point)`='AU0.25',
                                `Arbitrary Upstrap (0.50 Stopping Point)`='AU0.50',
                                `Arbitrary Upstrap (0.75 Stopping Point)`='AU0.75',
                                `Calibrated Upstrap`='CU',
                                `Calibrated Upstrap (0.25 Stopping Point)`='CU0.25',
                                `Calibrated Upstrap (0.50 Stopping Point)`='CU0.50',
                                `Calibrated Upstrap (0.75 Stopping Point)`='CU0.75',
                                `OBrien Fleming GSD`='OBF',
                                `Pocock GSD`='PO',
                                `Alpha Beta Spending`='GU',
                                `Alpha Beta Spending (0.25 Stopping Point)`='GU0.25',
                                `Alpha Beta Spending (0.50 Stopping Point)`='GU0.50',
                                `Alpha Beta Spending (0.75 Stopping Point)`='GU0.75',
                                `0.50-0.75 Arbitrary Upstrap`='AU0.50-0.75',
                                `0.50-0.75 Calibrated Upstrap`='CU0.50-0.75',
                                `0.50-0.75 AB Calbrated Upstrap`='GU0.50-0.75')) %>%
  select(n,power,method,monitor_type,ExpectedN_mean,ExpectedN_sd,Decision,Decision_25,
         Decision_50,Decision_75,RejectionRate) %>%
  filter(power %in% c(0.05,0.80)) %>%
  pivot_wider(names_from='power',values_from=c('ExpectedN_mean','ExpectedN_sd','Decision',
                                               'Decision_25','Decision_50','Decision_75',
                                               'RejectionRate'))

###############################################################
### DEFINE USER INPUTS                                      ###
###############################################################

#sample size
n_set <- 600

#monitoring type
mt_set <- 'FO'

#power/type I error boundaries
max_tie <- 0.05
min_power <- 0.75

#radar chart view preference
 #??? ifelse statement???

###############################################################
### POWER/TIE SCATTERPLOTS                                  ###
###############################################################

#modify raw dataset
dat_scatter <- dat_final_output %>% 
  filter(method != 'OBF') %>%
  filter(n == n_set & monitor_type == mt_set)

#define plotly helper functions for adding lines
vline <- function(x = 0, color = "green") {
  list(
    type = "line",
    y0 = 0,
    y1 = 1,
    yref = "paper",
    x0 = x,
    x1 = x,
    line = list(color = color, dash="dot")
  )
}
hline <- function(y = 0, color = "blue") {
  list(
    type = "line",
    x0 = 0,
    x1 = 1,
    xref = "paper",
    y0 = y,
    y1 = y,
    line = list(color = color, dash="dot")
  )
}

#ployly scatterplot
dat_scatter %>% 
  plot_ly(y = ~RejectionRate_0.8, 
          x = ~RejectionRate_0.05,
          type = 'scatter', mode = 'markers',
          color = ~method) %>%
  layout(title='Power vs Type I Error',
         annotations = list(
           list(text = paste("Minimum Power = ",min_power,sep = ''), 
                x=0,y=0,
                showarrow = F, xref='paper', yref='paper', 
                xanchor='left', yanchor='auto', xshift=0, yshift=0,
                font=list(size=8, color="blue")),
           list(text = paste("Maximum Type I Error = ",max_tie,sep = ''), 
                x=0,y=0.05,
                showarrow = F, xref='paper', yref='paper', 
                xanchor='left', yanchor='auto', xshift=0, yshift=0,
                font=list(size=8, color="green"))
           ),
         yaxis = list(title = 'Power',
                      range = c(min(min_power-0.05,dat_scatter$RejectionRate_0.8-0.05),
                                max(min_power+0.05,dat_scatter$RejectionRate_0.8+0.05))),
         xaxis = list(title = 'Type I Error',
                      range = c(min(max_tie-0.005,dat_scatter$RejectionRate_0.05-0.005),
                                max(max_tie+0.005,dat_scatter$RejectionRate_0.05+0.005))),
         shapes = list(hline(y = min_power),vline(x = max_tie),
                       list(type = "rect",
                            fillcolor = "red", line = list(color = "red"), opacity = 0.2,
                            y0 = min_power, y1 = max(min_power+0.05,dat_scatter$RejectionRate_0.8+0.05), 
                            x0 = min(max_tie-0.005,dat_scatter$RejectionRate_0.05-0.005), x1 = max_tie))
         ) #%>%
  add_text(showlegend = F, 
           x = min(max_tie+0.05,dat_scatter$RejectionRate_0.05+0.05), 
           y = min(min_power-0.044,dat_scatter$RejectionRate_0.8-0.044),
           text = paste("Minimum Power = ",min_power,sep = ''),
           textposition = "middle",
           textfont = list(color = 'blue', 
                           size = 8
                           )
           ) %>%
  
  add_text(showlegend = F, 
           x = min(max_tie+0.063,dat_scatter$RejectionRate_0.05+0.063), 
           y = min(min_power-0.037,dat_scatter$RejectionRate_0.8-0.037),
           text = paste("Maximum Type I Error = ",max_tie,sep = ''),
           textposition = "middle",
           textfont = list(color = 'green', 
                           size = 8
           )
  )

###############################################################
### EXPECTED SAMPLE SIZE BAR PLOTS                          ###
###############################################################

#modify raw dataset
dat_bar <- dat_final_output %>% 
  filter(method != 'OBF') %>%
  filter(n == n_set & monitor_type == mt_set) %>%
  #filter(RejectionRate_0.05<=max_tie & RejectionRate_0.8>=min_power) %>% #boundary conditions
  mutate(method = droplevels(method))

#ployly bar plot
barPlot1 <- dat_bar %>%
  plot_ly(
    x = ~method,
    y = ~ExpectedN_mean_0.05,
    color = ~method,
    legendgroup=~method,
    type = "bar"
  ) %>%
  layout(title='Expected Sample Size H0',
         yaxis = list(title = 'N'),
         xaxis = list(title = 'Method'),
         shapes = list(hline(filter(dat_final_output,n == n_set & monitor_type == mt_set & 
                                      method=='OBF')$ExpectedN_mean_0.05))
  )
barPlot2 <- dat_bar %>%
  plot_ly(
    x = ~method,
    y = ~ExpectedN_mean_0.8,
    color = ~method,
    legendgroup=~method,
    showlegend = F,
    type = "bar"
  ) %>%
  layout(title='Expected Sample Size H1',
         yaxis = list(title = 'N'),
         xaxis = list(title = 'Method'),
         shapes = list(hline(filter(dat_final_output,n == n_set & monitor_type == mt_set & 
                                      method=='OBF')$ExpectedN_mean_0.8))
  )
subplot(barPlot1,barPlot2)

###############################################################
### METHOD ENROLLMENT CURVES                                ###
###############################################################

#modify raw dataset
dat_enroll_all <- dat_final_output %>% 
  mutate(method = droplevels(method)) %>%
  mutate(ExpectedN_0_H0 = 0,
         ExpectedN_0.25_H0 = 0.25*n,
         ExpectedN_0.50_H0 = ceiling(Decision_25_0.05*(0.25*n)+(1-Decision_25_0.05)*(0.5*n)),
         ExpectedN_0.75_H0 = ceiling(Decision_25_0.05*(0.25*n)+Decision_50_0.05*(0.5*n)+(1-Decision_25_0.05-Decision_50_0.05)*(0.75*n)),
         ExpectedN_1.00_H0 = ceiling(Decision_25_0.05*(0.25*n)+Decision_50_0.05*(0.5*n)+Decision_75_0.05*(0.75*n)+(1-Decision_25_0.05-Decision_50_0.05-Decision_75_0.05)*n),
         ExpectedN_0_H1 = 0,
         ExpectedN_0.25_H1 = 0.25*n,
         ExpectedN_0.50_H1 = ceiling(Decision_25_0.8*(0.25*n)+(1-Decision_25_0.8)*(0.5*n)),
         ExpectedN_0.75_H1 = ceiling(Decision_25_0.8*(0.25*n)+Decision_50_0.8*(0.5*n)+(1-Decision_25_0.8-Decision_50_0.8)*(0.75*n)),
         ExpectedN_1.00_H1 = ceiling(Decision_25_0.8*(0.25*n)+Decision_50_0.8*(0.5*n)+Decision_75_0.8*(0.75*n)+(1-Decision_25_0.8-Decision_50_0.8-Decision_75_0.8)*n)) %>%
  group_by(n,method,monitor_type) %>%
  pivot_longer(cols = c('ExpectedN_0_H0','ExpectedN_0.25_H0','ExpectedN_0.50_H0','ExpectedN_0.75_H0','ExpectedN_1.00_H0'),
               names_to='info_fract',values_to='ExpectedN_H0') %>%
  mutate(info_fract = recode_factor(info_fract, ExpectedN_0_H0=0,ExpectedN_0.25_H0=0.25,
                                    ExpectedN_0.50_H0=0.50,ExpectedN_0.75_H0=0.75,ExpectedN_1.00_H0=1)) %>% 
  pivot_longer(cols = c('ExpectedN_0_H1','ExpectedN_0.25_H1','ExpectedN_0.50_H1','ExpectedN_0.75_H1','ExpectedN_1.00_H1'),
               names_to='info_fract2',values_to='ExpectedN_H1') %>%
  mutate(info_fract2 = recode_factor(info_fract2, ExpectedN_0_H1=0,ExpectedN_0.25_H1=0.25,
                                    ExpectedN_0.50_H1=0.50,ExpectedN_0.75_H1=0.75,ExpectedN_1.00_H1=1)) %>%
  filter(info_fract == info_fract2) %>% 
  filter(method != 'OBF') %>%
  filter(n == n_set & monitor_type == mt_set) %>%
  filter(RejectionRate_0.05<=max_tie & RejectionRate_0.8>=min_power | method == 'FS') %>% #boundary conditions
  select(n,method,monitor_type,info_fract,ExpectedN_H0,ExpectedN_H1)

dat_enroll_FS <- rbind(dat_final_output %>% 
  mutate(method = droplevels(method)) %>%
  mutate(ExpectedN_0_H0 = 0,
         ExpectedN_0.25_H0 = 0.25*n,
         ExpectedN_0.50_H0 = ceiling(Decision_25_0.05*(0.25*n)+(1-Decision_25_0.05)*(0.5*n)),
         ExpectedN_0.75_H0 = ceiling(Decision_25_0.05*(0.25*n)+Decision_50_0.05*(0.5*n)+(1-Decision_25_0.05-Decision_50_0.05)*(0.75*n)),
         ExpectedN_1.00_H0 = ceiling(Decision_25_0.05*(0.25*n)+Decision_50_0.05*(0.5*n)+Decision_75_0.05*(0.75*n)+(1-Decision_25_0.05-Decision_50_0.05-Decision_75_0.05)*n),
         ExpectedN_0_H1 = 0,
         ExpectedN_0.25_H1 = 0.25*n,
         ExpectedN_0.50_H1 = ceiling(Decision_25_0.8*(0.25*n)+(1-Decision_25_0.8)*(0.5*n)),
         ExpectedN_0.75_H1 = ceiling(Decision_25_0.8*(0.25*n)+Decision_50_0.8*(0.5*n)+(1-Decision_25_0.8-Decision_50_0.8)*(0.75*n)),
         ExpectedN_1.00_H1 = ceiling(Decision_25_0.8*(0.25*n)+Decision_50_0.8*(0.5*n)+Decision_75_0.8*(0.75*n)+(1-Decision_25_0.8-Decision_50_0.8-Decision_75_0.8)*n)) %>%
  group_by(n,method,monitor_type) %>%
  pivot_longer(cols = c('ExpectedN_0_H0','ExpectedN_0.25_H0','ExpectedN_0.50_H0','ExpectedN_0.75_H0','ExpectedN_1.00_H0'),
               names_to='info_fract',values_to='ExpectedN_H0') %>%
  mutate(info_fract = recode_factor(info_fract, ExpectedN_0_H0=0,ExpectedN_0.25_H0=0.25,
                                    ExpectedN_0.50_H0=0.50,ExpectedN_0.75_H0=0.75,ExpectedN_1.00_H0=1)) %>% 
  pivot_longer(cols = c('ExpectedN_0_H1','ExpectedN_0.25_H1','ExpectedN_0.50_H1','ExpectedN_0.75_H1','ExpectedN_1.00_H1'),
               names_to='info_fract2',values_to='ExpectedN_H1') %>%
  mutate(info_fract2 = recode_factor(info_fract2, ExpectedN_0_H1=0,ExpectedN_0.25_H1=0.25,
                                     ExpectedN_0.50_H1=0.50,ExpectedN_0.75_H1=0.75,ExpectedN_1.00_H1=1)) %>%
  filter(info_fract == info_fract2) %>% 
  filter(method != 'OBF') %>%
  filter(RejectionRate_0.05<=max_tie & RejectionRate_0.8>=min_power) %>% #boundary conditions
  select(n,method,monitor_type,info_fract,ExpectedN_H0,ExpectedN_H1),
  expand.grid(n = c(40,160,600,2000),method = 'FS',monitor_type = c('FO','EO','FE'),
              info_fract = c(0,0.25,0.50,0.75,1)) %>%
  mutate(ExpectedN_H0 = n*info_fract,ExpectedN_H1 = n*info_fract,
         info_fract = as.factor(info_fract))
  ) %>% filter(n == n_set & monitor_type == mt_set)


dat_enroll_OBF <- dat_final_output %>% 
  mutate(method = droplevels(method)) %>%
  mutate(ExpectedN_0_H0 = 0,
         ExpectedN_0.25_H0 = 0.25*n,
         ExpectedN_0.50_H0 = ceiling(Decision_25_0.05*(0.25*n)+(1-Decision_25_0.05)*(0.5*n)),
         ExpectedN_0.75_H0 = ceiling(Decision_25_0.05*(0.25*n)+Decision_50_0.05*(0.5*n)+(1-Decision_25_0.05-Decision_50_0.05)*(0.75*n)),
         ExpectedN_1.00_H0 = ceiling(Decision_25_0.05*(0.25*n)+Decision_50_0.05*(0.5*n)+Decision_75_0.05*(0.75*n)+(1-Decision_25_0.05-Decision_50_0.05-Decision_75_0.05)*n),
         ExpectedN_0_H1 = 0,
         ExpectedN_0.25_H1 = 0.25*n,
         ExpectedN_0.50_H1 = ceiling(Decision_25_0.8*(0.25*n)+(1-Decision_25_0.8)*(0.5*n)),
         ExpectedN_0.75_H1 = ceiling(Decision_25_0.8*(0.25*n)+Decision_50_0.8*(0.5*n)+(1-Decision_25_0.8-Decision_50_0.8)*(0.75*n)),
         ExpectedN_1.00_H1 = ceiling(Decision_25_0.8*(0.25*n)+Decision_50_0.8*(0.5*n)+Decision_75_0.8*(0.75*n)+(1-Decision_25_0.8-Decision_50_0.8-Decision_75_0.8)*n)) %>%
  group_by(n,method,monitor_type) %>%
  pivot_longer(cols = c('ExpectedN_0_H0','ExpectedN_0.25_H0','ExpectedN_0.50_H0','ExpectedN_0.75_H0','ExpectedN_1.00_H0'),
               names_to='info_fract',values_to='ExpectedN_H0') %>%
  mutate(info_fract = recode_factor(info_fract, ExpectedN_0_H0=0,ExpectedN_0.25_H0=0.25,
                                    ExpectedN_0.50_H0=0.50,ExpectedN_0.75_H0=0.75,ExpectedN_1.00_H0=1)) %>% 
  pivot_longer(cols = c('ExpectedN_0_H1','ExpectedN_0.25_H1','ExpectedN_0.50_H1','ExpectedN_0.75_H1','ExpectedN_1.00_H1'),
               names_to='info_fract2',values_to='ExpectedN_H1') %>%
  mutate(info_fract2 = recode_factor(info_fract2, ExpectedN_0_H1=0,ExpectedN_0.25_H1=0.25,
                                     ExpectedN_0.50_H1=0.50,ExpectedN_0.75_H1=0.75,ExpectedN_1.00_H1=1)) %>%
  filter(info_fract == info_fract2) %>% 
  filter(n == n_set & monitor_type == mt_set) %>%
  filter(RejectionRate_0.05<=max_tie & RejectionRate_0.8>=min_power | method == 'OBF') %>% #boundary conditions
  select(n,method,monitor_type,info_fract,ExpectedN_H0,ExpectedN_H1)

#all methods enrollment curve
palette <- brewer.pal(length(unique(dat_enroll_all$method)), "Dark2")
plot1 <- plot_ly(data = dat_enroll_all,legendgroup=~method,color=~method)%>%
  layout(xaxis = list(title = "Information Fraction"),
         yaxis = list (title = "Expected Sample Size (H0)",range = c(0,n_set))) %>%
  add_trace(data = dat_enroll_all, x = ~info_fract, y = ~ExpectedN_H0, name = ~method,
            type = 'scatter',
            mode = 'line+markers',
            line = list(width = 4)
            )
plot2 <- plot_ly(data = dat_enroll_all,legendgroup=~method,showlegend=F,color=~method)%>%
  layout(xaxis = list(title = "Information Fraction"),
         yaxis = list (title = "Expected Sample Size (H1)",range = c(0,n_set))) %>% 
  add_trace(data = dat_enroll_all, x = ~info_fract, y = ~ExpectedN_H1, name = ~method,
            type = 'scatter',
            mode = 'line+markers',
            line = list(width = 4)
            )
p = subplot(plot1,plot2,titleY = T,titleX=T,margin = 0.08) %>% 
  layout(title = 'H0 and H1 Cumulative Enrollment Curves')
p
  
#FS relative enrollment curve

#OBF relative enrollment curve


###############################################################
### METHOD PERFORMANCE RADAR PLOTS                          ###
###############################################################

#modify raw dataset
dat_radar <- dat_final_output %>% 
    ungroup() %>%
    filter(method != 'OBF') %>%
    filter(n == n_set & monitor_type == mt_set) %>%
    filter(RejectionRate_0.05<=max_tie & RejectionRate_0.8>=min_power) %>% #boundary conditions
    mutate(ExpectedN_mean_0.05 = ExpectedN_mean_0.05/n,
           ExpectedN_mean_0.8 = ExpectedN_mean_0.8/n) %>%
    select(method,ExpectedN_mean_0.05,ExpectedN_mean_0.8,Decision_0.05,Decision_0.8,RejectionRate_0.05,
           RejectionRate_0.8) %>%
    rename(ESS_H0 = ExpectedN_mean_0.05,ESS_H1 = ExpectedN_mean_0.8,Prop_H0 = Decision_0.05,
           Prop_H1 = Decision_0.8,TIE = RejectionRate_0.05,Power = RejectionRate_0.8)
dat_radar_OBF <- dat_final_output %>% 
  ungroup() %>%
  filter(method == 'OBF') %>%
  filter(n == n_set & monitor_type == mt_set) %>%
  filter(RejectionRate_0.05<=max_tie & RejectionRate_0.8>=min_power) %>% #boundary conditions
  mutate(ExpectedN_mean_0.05 = ExpectedN_mean_0.05/n,
         ExpectedN_mean_0.8 = ExpectedN_mean_0.8/n) %>%
  select(method,ExpectedN_mean_0.05,ExpectedN_mean_0.8,Decision_0.05,Decision_0.8,RejectionRate_0.05,
         RejectionRate_0.8) %>%
  rename(ESS_H0 = ExpectedN_mean_0.05,ESS_H1 = ExpectedN_mean_0.8,Prop_H0 = Decision_0.05,
         Prop_H1 = Decision_0.8,TIE = RejectionRate_0.05,Power = RejectionRate_0.8)

#all methods radar plot
radarPlot_single <- function(dat_radar){
  fig <- plot_ly(
    type = 'scatterpolar',
    fill = 'toself'
  ) 
  for (i in 1:nrow(dat_radar)) {
    fig <- fig %>%
      add_trace(
        r = dat_radar[i,2:ncol(dat_radar)] %>% as.numeric(),
        theta = colnames(dat_radar)[2:ncol(dat_radar)],
        name = dat_radar$method[i]
      ) 
  }
  fig <- fig %>%
    layout(
      polar = list(
        radialaxis = list(
          visible = T,
          range = c(0,1)
        )
      )
    )
  return(fig)
}

#small multiples radar plots
radarPlot_multiples <- function(dat_radar,dat_radar_OBF){
  fig_list <- list()
  for (i in 1:nrow(dat_radar)) {
    fig <- plot_ly(
      type = 'scatterpolar',
      fill = 'toself'
    ) 
    fig <- fig %>%
      add_trace(
        r = dat_radar[i,2:ncol(dat_radar)] %>% as.numeric(),
        theta = colnames(dat_radar)[2:ncol(dat_radar)],
        name = dat_radar$method[i]
      ) 
    fig <- fig %>%
      add_trace(
        r = dat_radar_OBF[,2:ncol(dat_radar)] %>% as.numeric(),
        theta = colnames(dat_radar)[2:ncol(dat_radar)],
        name = dat_radar_OBF$method
      ) 
    fig <- fig %>%
      layout(
        polar = list(
          radialaxis = list(
            visible = T,
            range = c(0,1)
          )
        )
      )
    fig_list[[i]] <- fig
  }
  return(fig_list)
}

#small multiple radar plots
radarPlot_multiple <- function(dat_radar,dat_radar_OBF,method_input){
  fig <- plot_ly(
    type = 'scatterpolar',
    fill = 'toself'
  ) 
  fig <- fig %>%
    add_trace(
      r = dat_radar %>% filter(method == method_input) %>% select(-method) %>% as.numeric(),
      theta = colnames(dat_radar)[2:ncol(dat_radar)],
      name = method_input
    ) 
  fig <- fig %>%
    add_trace(
      r = dat_radar_OBF %>% select(-method) %>% as.numeric(),
      theta = colnames(dat_radar %>% select(-method)),
      name = 'OBF'
    ) 
  fig <- fig %>%
    layout(
      polar = list(
        radialaxis = list(
          visible = T,
          range = c(0,1)
        )
      )
    )
  return(fig)
}

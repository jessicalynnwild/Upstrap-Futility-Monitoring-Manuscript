###############################################################
### file: app.R                                             ###
### authors: Jess Wild                                      ###
### creation date: 04/28/2023                               ###
### latest edit date:04/29/2023                             ###
### description: code to implement shiny app                ###
###############################################################

###############################################################
### SETUP WORKSPACE                                         ###
###############################################################

#clear workspace
remove(list = ls())

#set working directory
setwd('C:/Users/wildje/Desktop/Dissertation/Paper 1')

#read in and process the data
load('./01 Data/Upstrap 10-3-2022 Results/ConsolidatedFinalOutputData.RData')

#load required packages
library(shiny)
library(shinythemes)
library(tidyverse)
library(plotly)
library(fmsb)
library(RColorBrewer)

###############################################################
### PROCESS DATA FOR VISUALIZATION                          ###
###############################################################

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
                                                 'RejectionRate')) %>% 
    filter(method %in% c('AU','CU','GU','OBF','PO','AU0.50-0.75','CU0.50-0.75','GU0.50-0.75','AU0.25','AU0.50','AU0.75','CU0.25','CU0.50','CU0.75','GU0.25','GU0.50','GU0.75')) %>% #EDIT THIS!!!
    #filter(method %in% c('AU','CU','GU','OBF','PO')) %>% #EDIT THIS!!!
    droplevels()

#enrollment data
dat_enroll <- dat_final_output %>% 
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
    filter(info_fract == info_fract2)

#definition data
defDat <- function(dat_final_output,n_set,mt_set,max_tie,min_power){
    dat_def <- dat_final_output %>% 
        filter(method != 'OBF' & method != 'PO') %>%
        filter(n == n_set & monitor_type == mt_set) %>%
        filter(RejectionRate_0.05<=max_tie & RejectionRate_0.8>=min_power) %>% #boundary conditions
        mutate(method = droplevels(method))
    return(dat_def)
}

#scatter plot data
scatterDat <- function(dat_final_output,n_set,mt_set){
    dat_scatter <- dat_final_output %>% 
        filter(method != 'OBF' & method != 'PO') %>%
        filter(n == n_set & monitor_type == mt_set)
    return(dat_scatter)
}

#bar plot data
barDat <- function(dat_final_output,n_set,mt_set,max_tie,min_power){
    dat_bar <- dat_final_output %>% 
        mutate(FSExpectedN_mean_0.05 = ExpectedN_mean_0.05/n,
               FSExpectedN_mean_0.8 = ExpectedN_mean_0.8/n) %>%
        filter(method != 'OBF' & method != 'PO') %>%
        filter(n == n_set & monitor_type == mt_set) %>%
        mutate(OBFExpectedN_mean_0.05 = ExpectedN_mean_0.05/(
            dat_final_output %>% ungroup() %>% 
                filter(method == 'OBF' & n == n_set & monitor_type == mt_set) %>% 
                select(ExpectedN_mean_0.05) %>% as.numeric()),
            OBFExpectedN_mean_0.8 = ExpectedN_mean_0.8/(
                dat_final_output %>% ungroup() %>% 
                    filter(method == 'OBF' & n == n_set & monitor_type == mt_set) %>% 
                    select(ExpectedN_mean_0.8) %>% as.numeric())) %>%
        mutate(POExpectedN_mean_0.05 = ExpectedN_mean_0.05/(
            dat_final_output %>% ungroup() %>% 
                filter(method == 'PO' & n == n_set & monitor_type == mt_set) %>% 
                select(ExpectedN_mean_0.05) %>% as.numeric()),
            POExpectedN_mean_0.8 = ExpectedN_mean_0.8/(
                dat_final_output %>% ungroup() %>% 
                    filter(method == 'PO' & n == n_set & monitor_type == mt_set) %>% 
                    select(ExpectedN_mean_0.8) %>% as.numeric())) %>%
        filter(RejectionRate_0.05<=max_tie & RejectionRate_0.8>=min_power) %>% #boundary conditions
        mutate(method = droplevels(method))
    return(dat_bar)
}

#enrollment curve data
enrollDat <- function(dat_enroll,n_set,mt_set,max_tie,min_power){
    dat_enroll_all <- dat_enroll %>% 
        filter(method != 'OBF' & method != 'PO') %>%
        filter(n == n_set & monitor_type == mt_set) %>%
        filter(RejectionRate_0.05<=max_tie & RejectionRate_0.8>=min_power | method == 'FS') %>% #boundary conditions
        mutate(method = droplevels(method)) %>%
        select(n,method,monitor_type,info_fract,ExpectedN_H0,ExpectedN_H1)
    return(dat_enroll_all)
}

enrollFSDat <- function(dat_final_output,n_set,mt_set,max_tie,min_power,method_input){
    dat_enroll_FS <- rbind(dat_enroll %>% 
                               filter(method != 'OBF' & method != 'PO') %>%
                               filter(RejectionRate_0.05<=max_tie & RejectionRate_0.8>=min_power) %>% #boundary conditions
                               select(n,method,monitor_type,info_fract,ExpectedN_H0,ExpectedN_H1),
                           expand.grid(n = c(40,160,600,2000),method = 'FS',monitor_type = c('FO','EO','FE'),
                                       info_fract = c(0,0.25,0.50,0.75,1)) %>%
                               mutate(ExpectedN_H0 = n*info_fract,ExpectedN_H1 = n*info_fract,
                                      info_fract = as.factor(info_fract))
    ) %>% filter(n == n_set & monitor_type == mt_set) %>%
        filter(method == method_input | method == 'FS')
    return(dat_enroll_FS)
}

enrollOBFDat <- function(dat_final_output,n_set,mt_set,max_tie,min_power,method_input){
    dat_enroll_OBF <- dat_enroll %>% 
        filter(n == n_set & monitor_type == mt_set) %>%
        filter(RejectionRate_0.05<=max_tie & RejectionRate_0.8>=min_power | method == 'OBF') %>% #boundary conditions
        select(n,method,monitor_type,info_fract,ExpectedN_H0,ExpectedN_H1) %>%
        filter(method == method_input | method == 'OBF')
    return(dat_enroll_OBF)
}

enrollPODat <- function(dat_final_output,n_set,mt_set,max_tie,min_power,method_input){
    dat_enroll_PO <- dat_enroll %>% 
        filter(n == n_set & monitor_type == mt_set) %>%
        filter(RejectionRate_0.05<=max_tie & RejectionRate_0.8>=min_power | method == 'PO') %>% #boundary conditions
        select(n,method,monitor_type,info_fract,ExpectedN_H0,ExpectedN_H1) %>%
        filter(method == method_input | method == 'PO')
    return(dat_enroll_PO)
}

#radar chart data
radarDat <- function(dat_final_output,n_set,mt_set,max_tie,min_power){
    dat_radar <- dat_final_output %>% 
        ungroup() %>%
        filter(method != 'OBF' & method != 'PO') %>%
        filter(n == n_set & monitor_type == mt_set) %>%
        filter(RejectionRate_0.05<=max_tie & RejectionRate_0.8>=min_power) %>% #boundary conditions
        mutate(ExpectedN_mean_0.05 = ExpectedN_mean_0.05/n,
               ExpectedN_mean_0.8 = ExpectedN_mean_0.8/n) %>%
        select(method,ExpectedN_mean_0.05,ExpectedN_mean_0.8,Decision_0.05,Decision_0.8,RejectionRate_0.05,
               RejectionRate_0.8) %>%
        dplyr::rename(ESS_H0 = ExpectedN_mean_0.05,ESS_H1 = ExpectedN_mean_0.8,Prop_H0 = Decision_0.05,
               Prop_H1 = Decision_0.8,TIE = RejectionRate_0.05,Power = RejectionRate_0.8) %>%
        mutate(TIE_inverse = 1-TIE,
               Prop_H0_inverse = 1-Prop_H0,
               Prop_H1_inverse = 1-Prop_H1) %>%
        select(-TIE, -Prop_H0_inverse, -Prop_H1_inverse) #EDIT
    return(dat_radar)
}

radarOBFDat <- function(dat_final_output,n_set,mt_set){
    dat_radar_OBF <- dat_final_output %>% 
        ungroup() %>%
        filter(method == 'OBF') %>%
        filter(n == n_set & monitor_type == mt_set) %>%
        mutate(ExpectedN_mean_0.05 = ExpectedN_mean_0.05/n,
               ExpectedN_mean_0.8 = ExpectedN_mean_0.8/n) %>%
        select(method,ExpectedN_mean_0.05,ExpectedN_mean_0.8,Decision_0.05,Decision_0.8,RejectionRate_0.05,
               RejectionRate_0.8) %>%
        dplyr::rename(ESS_H0 = ExpectedN_mean_0.05,ESS_H1 = ExpectedN_mean_0.8,Prop_H0 = Decision_0.05,
               Prop_H1 = Decision_0.8,TIE = RejectionRate_0.05,Power = RejectionRate_0.8) %>%
        mutate(TIE_inverse = 1-TIE,
               Prop_H0_inverse = 1-Prop_H0,
               Prop_H1_inverse = 1-Prop_H1) %>%
        select(-TIE, -Prop_H0_inverse, -Prop_H1_inverse) #EDIT
    return(dat_radar_OBF)
}

radarPODat <- function(dat_final_output,n_set,mt_set){
    dat_radar_PO <- dat_final_output %>% 
        ungroup() %>%
        filter(method == 'PO') %>%
        filter(n == n_set & monitor_type == mt_set) %>%
        mutate(ExpectedN_mean_0.05 = ExpectedN_mean_0.05/n,
               ExpectedN_mean_0.8 = ExpectedN_mean_0.8/n) %>%
        select(method,ExpectedN_mean_0.05,ExpectedN_mean_0.8,Decision_0.05,Decision_0.8,RejectionRate_0.05,
               RejectionRate_0.8) %>%
        dplyr::rename(ESS_H0 = ExpectedN_mean_0.05,ESS_H1 = ExpectedN_mean_0.8,Prop_H0 = Decision_0.05,
                      Prop_H1 = Decision_0.8,TIE = RejectionRate_0.05,Power = RejectionRate_0.8) %>%
        mutate(TIE_inverse = 1-TIE,
               Prop_H0_inverse = 1-Prop_H0,
               Prop_H1_inverse = 1-Prop_H1) %>%
        select(-TIE, -Prop_H0_inverse, -Prop_H1_inverse) #EDIT
    return(dat_radar_PO)
}

###############################################################
### VISUALIZATION HELPER FUNCTIONS                          ###
###############################################################

#line functions
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

#scatterplot functions
scatterPlot <- function(dat_scatter,max_tie,min_power){
    plot <- dat_scatter %>% 
        plot_ly(y = ~RejectionRate_0.8, 
                x = ~RejectionRate_0.05,
                type = 'scatter', mode = 'markers',
                color = ~method) %>%
        layout(title='<b>Power vs Type I Error</b>',
               annotations = list(
                   list(text = paste("Minimum Power Threshold = ",min_power,sep = ''), 
                        x=0,y=0,
                        showarrow = F, xref='paper', yref='paper', 
                        xanchor='left', yanchor='auto', xshift=0, yshift=0,
                        font=list(size=8, color="blue")),
                   list(text = paste("Maximum Type I Error Threshold = ",max_tie,sep = ''), 
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
               shapes = list(hline(min_power),vline(max_tie),
                             list(type = "rect",
                                  fillcolor = "green", line = list(color = "green"), opacity = 0.05,
                                  y0 = min_power, y1 = max(min_power+0.05,dat_scatter$RejectionRate_0.8+0.05), 
                                  x0 = min(max_tie-0.005,dat_scatter$RejectionRate_0.05-0.005), x1 = max_tie))
        ) 
    return(plot)
}

#barplot functions
barPlot <- function(dat_bar,n_set,mt_set){
    plot1 <- dat_bar %>%
        plot_ly(
            x = ~method,
            y = ~ExpectedN_mean_0.8,
            color = ~method,
            legendgroup=~method,
            type = "bar",
            error_y = ~list(array = ExpectedN_sd_0.8,color = '#000000')
        ) %>%
        layout(showlegend = FALSE,
               yaxis = list(title = 'Expected Sample Size (H1)',range = c(0,n_set)),
               xaxis = list(title = 'Method'),
               shapes = list(hline(filter(dat_final_output,n == n_set & monitor_type == mt_set & 
                                              method=='OBF')$ExpectedN_mean_0.8)),
               annotations = list(
                   list(text = paste("OBF = ",filter(dat_final_output,n == n_set & monitor_type == mt_set & 
                                                         method=='OBF')$ExpectedN_mean_0.8,sep = ''), 
                        x=0,y=1.0001,
                        showarrow = F, xref='paper', yref='paper', 
                        xanchor='left', yanchor='auto', xshift=0, yshift=0,
                        font=list(size=8, color="blue")))
        )
    plot2 <- dat_bar %>%
        plot_ly(
            x = ~method,
            y = ~ExpectedN_mean_0.05,
            color = ~method,
            legendgroup=~method,
            showlegend = F,
            type = "bar",
            error_y = ~list(array = ExpectedN_sd_0.05,color = '#000000')
        ) %>%
        layout(yaxis = list(title = 'Expected Sample Size (H0)',range = c(0,n_set)),
               xaxis = list(title = 'Method'),
               shapes = list(hline(filter(dat_final_output,n == n_set & monitor_type == mt_set & 
                                              method=='OBF')$ExpectedN_mean_0.05)),
               annotations = list(
                   list(text = paste("OBF = ",filter(dat_final_output,n == n_set & monitor_type == mt_set & 
                                                         method=='OBF')$ExpectedN_mean_0.05,sep = ''), 
                        x=0,y=1.0001,
                        showarrow = F, xref='paper', yref='paper', 
                        xanchor='left', yanchor='auto', xshift=0, yshift=0,
                        font=list(size=8, color="blue")))
        )
    plot = subplot(plot2,plot1,titleY = TRUE,margin = 0.08) %>% 
        layout(title = '<b>Expected Sample Size</b>')
    return(plot)
}

FSbarPlot <- function(dat_bar,n_set,mt_set){
    plot1 <- dat_bar %>%
        plot_ly(
            x = ~method,
            y = ~FSExpectedN_mean_0.8,
            color = ~method,
            legendgroup=~method,
            type = "bar"
        ) %>%
        layout(showlegend = FALSE,
               yaxis = list(title = 'Expected Sample Size (H1)',range = c(0,1)),
               xaxis = list(title = 'Method')
        )
    plot2 <- dat_bar %>%
        plot_ly(
            x = ~method,
            y = ~FSExpectedN_mean_0.05,
            color = ~method,
            legendgroup=~method,
            showlegend = F,
            type = "bar"
        ) %>%
        layout(yaxis = list(title = 'Expected Sample Size (H0)',range = c(0,1)),
               xaxis = list(title = 'Method')
        )
    plot = subplot(plot2,plot1,titleY = TRUE,margin = 0.08) %>% 
        layout(title = '<b>Expected Sample Size Proportional to Fixed Sample Design</b>')
    return(plot)
}

OBFbarPlot <- function(dat_bar,n_set,mt_set){
    plot1 <- dat_bar %>%
        plot_ly(
            x = ~method,
            y = ~OBFExpectedN_mean_0.8,
            color = ~method,
            legendgroup=~method,
            type = "bar"
        ) %>%
        layout(showlegend = FALSE,
               yaxis = list(title = 'Expected Sample Size (H1)',range = c(0,1)),
               xaxis = list(title = 'Method')
        )
    plot2 <- dat_bar %>%
        plot_ly(
            x = ~method,
            y = ~OBFExpectedN_mean_0.05,
            color = ~method,
            legendgroup=~method,
            showlegend = F,
            type = "bar"
        ) %>%
        layout(yaxis = list(title = 'Expected Sample Size (H0)',range = c(0,1)),
               xaxis = list(title = 'Method')
        )
    plot = subplot(plot2,plot1,titleY = TRUE,margin = 0.08) %>% 
        layout(title = '<b>Expected Sample Size Proportional to OBF</b>')
    return(plot)
}

PObarPlot <- function(dat_bar,n_set,mt_set){
    plot1 <- dat_bar %>%
        plot_ly(
            x = ~method,
            y = ~POExpectedN_mean_0.8,
            color = ~method,
            legendgroup=~method,
            type = "bar"
        ) %>%
        layout(showlegend = FALSE,
               yaxis = list(title = 'Expected Sample Size (H1)',range = c(0,max(rbind(dat_bar$POExpectedN_mean_0.05,dat_bar$POExpectedN_mean_0.8)))),
               xaxis = list(title = 'Method')
        )
    plot2 <- dat_bar %>%
        plot_ly(
            x = ~method,
            y = ~POExpectedN_mean_0.05,
            color = ~method,
            legendgroup=~method,
            showlegend = F,
            type = "bar"
        ) %>%
        layout(yaxis = list(title = 'Expected Sample Size (H0)',range = c(0,max(rbind(dat_bar$POExpectedN_mean_0.05,dat_bar$POExpectedN_mean_0.8)))),
               xaxis = list(title = 'Method')
        )
    plot = subplot(plot2,plot1,titleY = TRUE,margin = 0.08) %>% 
        layout(title = '<b>Expected Sample Size Proportional to PO</b>')
    return(plot)
}

#enrollment curve functions
enrollmentPlot <- function(dat_enroll_all,n_set){
    plot1 <- plot_ly(data = dat_enroll_all,legendgroup=~method,color=~method)%>%
        layout(xaxis = list(title = "Information Fraction"),
               yaxis = list (title = "Expected Sample Size (H0)",range = c(0,n_set))) %>%
        add_trace(data = dat_enroll_all, x = ~info_fract, y = ~ExpectedN_H0, name = ~method,
                  type = 'scatter',
                  #error_y = ~list(array = ExpectedN_sd_0.05,color = '#000000'),
                  mode = 'line+markers',
                  line = list(width = 4)
        )
    plot2 <- plot_ly(data = dat_enroll_all,legendgroup=~method,color=~method,showlegend=F)%>%
        layout(xaxis = list(title = "Information Fraction"),
               yaxis = list (title = "Expected Sample Size (H1)",range = c(0,n_set))) %>% 
        add_trace(data = dat_enroll_all, x = ~info_fract, y = ~ExpectedN_H1, name = ~method,
                  type = 'scatter',
                  #error_y = ~list(array = ExpectedN_sd_0.8,color = '#000000')
                  mode = 'line+markers',
                  line = list(width = 4)
        )
    plot = subplot(plot1,plot2,titleY = T,titleX = T,margin = 0.08) %>% 
        layout(title = '<b>Cumulative Enrollment Curves</b>')
    return(plot)
}

FSenrollmentPlot <- function(dat_enroll_FS,n_set,subplot_set){
    plot1 <- plot_ly(data = dat_enroll_FS %>% filter(method == 'FS'),x = ~info_fract, y = ~ExpectedN_H0, name = ~method,
                     type = 'scatter',
                     mode = 'line+markers',
                     line = list(color = ~method, width = 4),
                     legendgroup=~method,color=~method)%>%
        layout(xaxis = list(title = "Information Fraction"),
               yaxis = list (title = "Expected Sample Size (H0)",range = c(0,n_set))) %>%
        add_trace(data = dat_enroll_FS %>% filter(method == subplot_set), x = ~info_fract, y = ~ExpectedN_H0, name = ~method,
                  type = 'scatter',
                  mode = 'line+markers',
                  fill = 'tonexty', 
                  fillcolor='rgba(0,100,80,0.2)',
                  line = list(width = 4)
        )
    plot2 <- plot_ly(data = dat_enroll_FS %>% filter(method == 'FS'),x = ~info_fract, y = ~ExpectedN_H1, name = ~method,
                     type = 'scatter',
                     mode = 'line+markers',
                     line = list(width = 4),
                     legendgroup=~method,color=~method,showlegend=F)%>%
        layout(xaxis = list(title = "Information Fraction"),
               yaxis = list (title = "Expected Sample Size (H1)",range = c(0,n_set))) %>%
        add_trace(data = dat_enroll_FS %>% filter(method == subplot_set), x = ~info_fract, y = ~ExpectedN_H1, name = ~method,
                  type = 'scatter',
                  mode = 'line+markers',
                  fill = 'tonexty', 
                  fillcolor='rgba(0,100,80,0.2)',
                  line = list(width = 4)
        )
    plot = subplot(plot1,plot2,titleY = T,titleX = T,margin = 0.08) %>% 
        layout(title = '<b>Cumulative Enrollment Curve Compared to Fixed Sample Design</b>')
    return(plot)
}

OBFenrollmentPlot <- function(dat_enroll_OBF,n_set,subplot_set){
    plot1 <- plot_ly(data = dat_enroll_OBF %>% filter(method == 'OBF'),x = ~info_fract, y = ~ExpectedN_H0, name = ~method,
                     type = 'scatter',
                     mode = 'line+markers',
                     line = list(width = 4),
                     legendgroup=~method,color=~method)%>%
        layout(xaxis = list(title = "Information Fraction"),
               yaxis = list (title = "Expected Sample Size (H0)",range = c(0,n_set))) %>%
        add_trace(data = dat_enroll_OBF %>% filter(method == subplot_set), x = ~info_fract, y = ~ExpectedN_H0, name = ~method,
                  type = 'scatter',
                  mode = 'line+markers',
                  fill = 'tonexty', 
                  fillcolor='rgba(0,100,80,0.2)',
                  line = list(width = 4)
        )
    plot2 <- plot_ly(data = dat_enroll_OBF %>% filter(method == 'OBF'),x = ~info_fract, y = ~ExpectedN_H1, name = ~method,
                     type = 'scatter',
                     mode = 'line+markers',
                     line = list(width = 4),
                     legendgroup=~method,color=~method,showlegend=F)%>%
        layout(xaxis = list(title = "Information Fraction"),
               yaxis = list (title = "Expected Sample Size (H1)",range = c(0,n_set))) %>%
        add_trace(data = dat_enroll_OBF %>% filter(method == subplot_set), x = ~info_fract, y = ~ExpectedN_H1, name = ~method,
                  type = 'scatter',
                  mode = 'line+markers',
                  fill = 'tonexty', 
                  fillcolor='rgba(0,100,80,0.2)',
                  line = list(width = 4)
        )
    plot = subplot(plot1,plot2,titleY = T,titleX = T,margin = 0.08) %>% 
        layout(title = '<b>Cumulative Enrollment Curve Compared to OBF</b>')
    return(plot)
}

POenrollmentPlot <- function(dat_enroll_PO,n_set,subplot_set){
    plot1 <- plot_ly(data = dat_enroll_PO %>% filter(method == 'PO'),x = ~info_fract, y = ~ExpectedN_H0, name = ~method,
                     type = 'scatter',
                     mode = 'line+markers',
                     line = list(width = 4),
                     legendgroup=~method,color=~method)%>%
        layout(xaxis = list(title = "Information Fraction"),
               yaxis = list (title = "Expected Sample Size (H0)",range = c(0,n_set))) %>%
        add_trace(data = dat_enroll_PO %>% filter(method == subplot_set), x = ~info_fract, y = ~ExpectedN_H0, name = ~method,
                  type = 'scatter',
                  mode = 'line+markers',
                  fill = 'tonexty', 
                  fillcolor='rgba(0,100,80,0.2)',
                  line = list(width = 4)
        )
    plot2 <- plot_ly(data = dat_enroll_PO %>% filter(method == 'PO'),x = ~info_fract, y = ~ExpectedN_H1, name = ~method,
                     type = 'scatter',
                     mode = 'line+markers',
                     line = list(width = 4),
                     legendgroup=~method,color=~method,showlegend=F)%>%
        layout(xaxis = list(title = "Information Fraction"),
               yaxis = list (title = "Expected Sample Size (H1)",range = c(0,n_set))) %>%
        add_trace(data = dat_enroll_PO %>% filter(method == subplot_set), x = ~info_fract, y = ~ExpectedN_H1, name = ~method,
                  type = 'scatter',
                  mode = 'line+markers',
                  fill = 'tonexty', 
                  fillcolor='rgba(0,100,80,0.2)',
                  line = list(width = 4)
        )
    plot = subplot(plot1,plot2,titleY = T,titleX = T,margin = 0.08) %>% 
        layout(title = '<b>Cumulative Enrollment Curve Compared to PO</b>')
    return(plot)
}

#radar plot functions
radarPlot <- function(dat_radar){
    fig <- plot_ly(
        type = 'scatterpolar',
        fill = 'toself'
    ) 
    for (i in 1:nrow(dat_radar)) {
        fig <- fig %>%
            add_trace(
                r = dat_radar[i,2:ncol(dat_radar)] %>% as.numeric(),
                theta = c('Expected Sample Size (H0)','Expected Sample Size (H1)',
                          'Interim Stopping Rate (H0)','Interim Stopping Rate (H1)','1 - Type I Error','Power'),#colnames(dat_radar)[2:ncol(dat_radar)],
                name = dat_radar$method[i]
            ) 
    }
    fig <- fig %>%
        layout(
            title = '<b>Combined Method Performance Radar Plot</b>',
            polar = list(
                radialaxis = list(
                    visible = T,
                    range = c(0,1)
                )
            )
        )
    return(fig)
}

OBFradarPlot <- function(dat_radar,dat_radar_OBF,method_input){
    fig <- plot_ly(
        type = 'scatterpolar',
        fill = 'toself'
    ) 
    fig <- fig %>%
        add_trace(
            r = dat_radar %>% filter(method == method_input) %>% select(-method) %>% as.numeric(),
            theta = c('Expected Sample Size (H0)','Expected Sample Size (H1)',
                              'Interim Stopping Rate (H0)','Interim Stopping Rate (H1)','1 - Type I Error','Power'),
            name = method_input
        ) 
    fig <- fig %>%
        add_trace(
            r = dat_radar_OBF %>% select(-method) %>% as.numeric(),
            theta = c('Expected Sample Size (H0)','Expected Sample Size (H1)',
                              'Interim Stopping Rate (H0)','Interim Stopping Rate (H1)','1 - Type I Error','Power'),
            name = 'OBF'
        ) 
    fig <- fig %>%
        layout(
            title = '<b>Method Performance Compared to OBF Radar Plot</b>',
            polar = list(
                radialaxis = list(
                    visible = T,
                    range = c(0,1)
                )
            )
        )
    return(fig)
}

POradarPlot <- function(dat_radar,dat_radar_PO,method_input){
    fig <- plot_ly(
        type = 'scatterpolar',
        fill = 'toself'
    ) 
    fig <- fig %>%
        add_trace(
            r = dat_radar %>% filter(method == method_input) %>% select(-method) %>% as.numeric(),
            theta = c('Expected Sample Size (H0)','Expected Sample Size (H1)',
                              'Interim Stopping Rate (H0)','Interim Stopping Rate (H1)','1 - Type I Error','Power'),
            name = method_input
        ) 
    fig <- fig %>%
        add_trace(
            r = dat_radar_PO %>% select(-method) %>% as.numeric(),
            theta = c('Expected Sample Size (H0)','Expected Sample Size (H1)',
                              'Interim Stopping Rate (H0)','Interim Stopping Rate (H1)','1 - Type I Error','Power'),
            name = 'PO'
        ) 
    fig <- fig %>%
        layout(
            title = '<b>Method Performance Compared to PO Radar Plot</b>',
            polar = list(
                radialaxis = list(
                    visible = T,
                    range = c(0,1)
                )
            )
        )
    return(fig)
}

###############################################################
### User Interface Code                                     ###
###############################################################

#define UI for application
ui <- fluidPage(

    #set theme
    theme = shinytheme("cerulean"),
    
    #set application title
    titlePanel("Interim Monitoring Method Performance Comparison"),

    #define sidebar with user inputs 
    sidebarLayout(
        sidebarPanel(
            sliderInput("min_power",
                        "Minimum Power:",
                        min = 0.5,
                        max = 1,
                        value = 0.6),
            sliderInput("max_tie",
                        "Maximum Type I Error",
                        min = 0,
                        max = 0.5,
                        value = 0.20),
            radioButtons("n_set", label = "Total Sample Size",
                        choices = c(40,160,600,2000), 
                        selected = 600),
            radioButtons("mt_set", label = "Monitoring Type",
                        choices = c('FO','EO','FE'), 
                        selected = 'FO'),
            uiOutput("subplot_set")
        ),
        

        #define tab structure
        mainPanel(
            navbarPage('',
                       #tabPanel("Project Background"),
                       tabPanel("Power & Type I Error Thresholds", 
                                plotlyOutput("scatterPlot")),
                       tabPanel("Expected Sample Size", 
                                plotlyOutput("barPlot"),
                                plotlyOutput("FSbarPlot"),
                                plotlyOutput("OBFbarPlot"),
                                plotlyOutput("PObarPlot")),
                       tabPanel("Cumulative Enrollment Curves", 
                                plotlyOutput("enrollPlot"),
                                plotlyOutput("FSenrollPlot"),
                                plotlyOutput("OBFenrollPlot"),
                                plotlyOutput("POenrollPlot")),
                       tabPanel("Performance Metric Radar Plots", 
                                plotlyOutput("radarPlot"),
                                plotlyOutput("OBFradarPlot"),
                                plotlyOutput("POradarPlot"))#,
                       #tabPanel("User Generated Simulation Study",
                       #         img(src='underConstruction.jpg', height="100%", width="100%"))
           )
        )
    )
)


###############################################################
### SERVER CODE                                             ###
###############################################################

# Define server logic required to draw plots
server <- function(input, output) {
    
    output$subplot_set <- renderUI({
        
        #define user inputs 
        max_tie <- input$max_tie
        min_power <- input$min_power
        n_set <- input$n_set
        mt_set <- input$mt_set
        
        #define dataset (just to get admissible methods)
        dat_def <- defDat(dat_final_output,n_set,mt_set,max_tie,min_power)
        
        #error message
        validate(
            need(nrow(dat_def)>0, "ERROR: Please select threshold values with at lease one admissible method")
        )
        
        
        selectInput("subplot_set", label = "Individual Plot Method",
                    choices = unique(dat_def$method),
                    selected = character(0))
    })
    
    output$scatterPlot <- renderPlotly({
        
        #define user inputs 
        max_tie <- input$max_tie
        min_power <- input$min_power
        n_set <- input$n_set
        mt_set <- input$mt_set
        
        #define dataset
        dat_scatter <- scatterDat(dat_final_output,n_set,mt_set)
        
        #plotly scatterplot
        scatterPlot(dat_scatter,max_tie,min_power)
    })
    
    output$barPlot <- renderPlotly({
        
        #define user inputs 
        max_tie <- input$max_tie
        min_power <- input$min_power
        n_set <- input$n_set
        mt_set <- input$mt_set
        
        #define dataset
        dat_bar <- barDat(dat_final_output,n_set,mt_set,max_tie,min_power)
        
        #error message
        validate(
            need(nrow(dat_bar)>0, "ERROR: Please select threshold values with at lease one admissible method")
        )
        
        #plotly barplot
        barPlot(dat_bar,n_set,mt_set)
    })
    
    output$FSbarPlot <- renderPlotly({
        
        #define user inputs 
        max_tie <- input$max_tie
        min_power <- input$min_power
        n_set <- input$n_set
        mt_set <- input$mt_set
        
        #define dataset
        dat_bar <- barDat(dat_final_output,n_set,mt_set,max_tie,min_power)
        
        #error message
        validate(
            need(nrow(dat_bar)>0, "ERROR: Please select threshold values with at lease one admissible method")
        )
        
        #plotly barplot
        FSbarPlot(dat_bar,n_set,mt_set)
    })
    
    output$OBFbarPlot <- renderPlotly({
        
        #define user inputs 
        max_tie <- input$max_tie
        min_power <- input$min_power
        n_set <- input$n_set
        mt_set <- input$mt_set
        
        #define dataset
        dat_bar <- barDat(dat_final_output,n_set,mt_set,max_tie,min_power)
        
        #error message
        validate(
            need(nrow(dat_bar)>0, "ERROR: Please select threshold values with at lease one admissible method")
        )
        
        #plotly barplot
        OBFbarPlot(dat_bar,n_set,mt_set)
    })
    
    output$PObarPlot <- renderPlotly({
        
        #define user inputs 
        max_tie <- input$max_tie
        min_power <- input$min_power
        n_set <- input$n_set
        mt_set <- input$mt_set
        
        #define dataset
        dat_bar <- barDat(dat_final_output,n_set,mt_set,max_tie,min_power)
        
        #error message
        validate(
            need(nrow(dat_bar)>0, "ERROR: Please select threshold values with at lease one admissible method")
        )
        
        #plotly barplot
        PObarPlot(dat_bar,n_set,mt_set)
    })
    
    output$enrollPlot <- renderPlotly({
        
        #define user inputs 
        max_tie <- input$max_tie
        min_power <- input$min_power
        n_set <- input$n_set
        mt_set <- input$mt_set
        
        #define dataset
        dat_enroll_all <- enrollDat(dat_enroll,n_set,mt_set,max_tie,min_power)
        
        #error message
        validate(
            need(nrow(dat_enroll_all)>0, "ERROR: Please select threshold values with at lease one admissible method")
        )
        
        #plotly enrollment curve
        enrollmentPlot(dat_enroll_all,n_set)
    })
    
    output$FSenrollPlot <- renderPlotly({
        
        #define user inputs 
        max_tie <- input$max_tie
        min_power <- input$min_power
        n_set <- input$n_set
        mt_set <- input$mt_set
        subplot_set <- input$subplot_set
        
        #define dataset
        dat_enroll_FS <- enrollFSDat(dat_enroll,n_set,mt_set,max_tie,min_power,subplot_set)
        
        #error message
        validate(
            need(length(unique(dat_enroll_FS$method))>1, "ERROR: Please select threshold values with at lease one admissible method")
        )
        
        #plotly enrollment curve
        FSenrollmentPlot(dat_enroll_FS,n_set,subplot_set)
    })
    
    output$OBFenrollPlot <- renderPlotly({
        
        #define user inputs 
        max_tie <- input$max_tie
        min_power <- input$min_power
        n_set <- input$n_set
        mt_set <- input$mt_set
        subplot_set <- input$subplot_set
        
        #define dataset
        dat_enroll_OBF <- enrollOBFDat(dat_enroll,n_set,mt_set,max_tie,min_power,subplot_set)
        
        #error message
        validate(
            need(length(unique(dat_enroll_OBF$method))>1, "ERROR: Please select threshold values with at lease one admissible method")
        )
        
        #plotly enrollment curve
        OBFenrollmentPlot(dat_enroll_OBF,n_set,subplot_set)
    })
    
    output$POenrollPlot <- renderPlotly({
        
        #define user inputs 
        max_tie <- input$max_tie
        min_power <- input$min_power
        n_set <- input$n_set
        mt_set <- input$mt_set
        subplot_set <- input$subplot_set
        
        #define dataset
        dat_enroll_PO <- enrollPODat(dat_enroll,n_set,mt_set,max_tie,min_power,subplot_set)
        
        #error message
        validate(
            need(length(unique(dat_enroll_PO$method))>1, "ERROR: Please select threshold values with at lease one admissible method")
        )
        
        #plotly enrollment curve
        POenrollmentPlot(dat_enroll_PO,n_set,subplot_set)
    })
    
    output$radarPlot <- renderPlotly({
        
        #define user inputs 
        max_tie <- input$max_tie
        min_power <- input$min_power
        n_set <- input$n_set
        mt_set <- input$mt_set
        
        #define dataset
        dat_radar <- radarDat(dat_final_output,n_set,mt_set,max_tie,min_power)
        
        #error message
        validate(
            need(nrow(dat_radar)>0, "ERROR: Please select threshold values with at lease one admissible method")
        )
        
        #plotly radarplot
        radarPlot(dat_radar)
    })
    
    output$OBFradarPlot <- renderPlotly({
        
        #define user inputs 
        max_tie <- input$max_tie
        min_power <- input$min_power
        n_set <- input$n_set
        mt_set <- input$mt_set
        subplot_set <- input$subplot_set
        
        #define dataset
        dat_radar <- radarDat(dat_final_output,n_set,mt_set,max_tie,min_power)
        dat_radar_OBF <- radarOBFDat(dat_final_output,n_set,mt_set)
        
        #error message
        validate(
            need(nrow(dat_radar)>0&nrow(dat_radar_OBF)>0, "ERROR: Please select threshold values with at lease one admissible method")
        )
        
        #plotly radarplot
        OBFradarPlot(dat_radar,dat_radar_OBF,method_input = subplot_set)
    })
    
    output$POradarPlot <- renderPlotly({
        
        #define user inputs 
        max_tie <- input$max_tie
        min_power <- input$min_power
        n_set <- input$n_set
        mt_set <- input$mt_set
        subplot_set <- input$subplot_set
        
        #define dataset
        dat_radar <- radarDat(dat_final_output,n_set,mt_set,max_tie,min_power)
        dat_radar_PO <- radarPODat(dat_final_output,n_set,mt_set)
        
        #error message
        validate(
            need(nrow(dat_radar)>0&nrow(dat_radar_PO)>0, "ERROR: Please select threshold values with at lease one admissible method")
        )
        
        #plotly radarplot
        POradarPlot(dat_radar,dat_radar_PO,method_input = subplot_set)
    })
}

#run the application 
shinyApp(ui = ui, server = server)



library(tidyverse)
library(lme4)
library(lmerTest)
library(kableExtra)
repoPath = "G:/My Drive/GitHub/imageAB_submission"
setwd(repoPath) 
source('AB_pic_helperFuncs.R')

#data file path for experiment 1 data: 
datFilePath = paste(getwd(), '/Exp1_rawData', sep = '')


#import the data as a list of data.frames
allDat = readCCdat(datFilePath)
df = allDat[[2]]
allDat = allDat[[1]]

setwd(paste0(repoPath, '/Figs')) 


##############    Figure 2b        #######################################
#function plots the T2|T1 proportion across lags averaged over participants
#and also adds summary variables to the dataframe at the participant level
df = plotBasicAB(allDat, df,3)

#extract the T2|T1 data from the data frame to ready for plotting at the 
#participant level
plotDat <- df[,grepl( 'T2.T1',colnames(df))] 
plotDat$subID <- seq(1,length(plotDat[,1]), 1)
plotDat <- plotDat %>%
            pivot_longer(cols = starts_with('lag'),
                        names_to = "lag", 
                        values_to = "value")

meanVals <- plotDat %>% 
            group_by(lag) %>% 
            summarise(mean_lag = mean(value), 
                      ste_lag = sd(value)/sqrt(13) * qt(.975, 13))
meanVals$subID = 99


xVals = unique(plotDat$lag)
png("pilot_T2.png",         # File name
    width=500, height=500)
plotDat %>% ggplot(aes(x=lag, y=value, group = subID)) +
  geom_line( position = position_dodge(.2), size = 1, alpha = .15) + 
  geom_point( position = position_dodge(.2)) +
  scale_x_discrete(labels = c('1','2','3','4','5','6','7')) +
  scale_y_continuous(breaks = seq(.2, 1.0, .1), limits = c(.2,1)) +
  theme_classic() +
  theme(axis.line = element_line(color = 'black', size = 3),
        axis.ticks = element_line(colour = "black", size = 2),
        axis.ticks.length=unit(-.25, "cm"),
        text = element_text(size = 20)) + 
  geom_line(data = meanVals, aes(x= lag, y = mean_lag), linewidth = 4) + 
  geom_point(data = meanVals, aes(x= lag, y = mean_lag), size = 10) +
  geom_errorbar(data = meanVals, aes(x=lag, y = mean_lag, 
                                     ymin = mean_lag - ste_lag,
                                     ymax = mean_lag + ste_lag),
                width = .5, linewidth = 2) + 
  ylab('T2 detection p(T2|T1)')
  # ylim(c(.4, 1))
  # rect(1, 5, 3, 7, col="white")
  dev.off()
  
  #apply linear mixed effects model
  anova_result <- lmer(value ~ lag + (1 | subID), data = plotDat)
  anova(anova_result)
  
  
  ####      Figure 2a #######################################################
# T1 detection accuracy 
  plotDat <- df[,!grepl( 'T2.T1',colnames(df)) & grepl('T1',colnames(df))] 
  plotDat$subID <- seq(1,length(plotDat[,1]), 1)
  plotDat <- plotDat %>%
    pivot_longer(cols = starts_with('lag'),
                 names_to = "lag", 
                 values_to = "value")
  
  meanVals <- plotDat %>% 
    group_by(lag) %>% 
    summarise(mean_lag = mean(value), 
              ste_lag = sd(value)/sqrt(13) * qt(.975, 13))
  meanVals$subID = 99
  
  
  xVals = unique(plotDat$lag)
  png("pilot_T1.png",         # File name
      width=500, height=500)
  plotDat %>% ggplot(aes(x=lag, y=value, group = subID)) +
    geom_line( position = position_dodge(.2), size = 1, alpha = .15) + 
    geom_point( position = position_dodge(.2)) +
    scale_x_discrete(labels = c('1','2','3','4','5','6','7')) +
    scale_y_continuous(breaks = seq(.2, 1.0, .1), limits = c(.2,1)) +
    theme_classic() +
    theme(axis.line = element_line(color = 'black', size = 3),
          axis.ticks = element_line(colour = "black", size = 2),
          axis.ticks.length=unit(-.25, "cm"),
          text = element_text(size = 20)) + 
    geom_line(data = meanVals, aes(x= lag, y = mean_lag), linewidth = 4) + 
    geom_point(data = meanVals, aes(x= lag, y = mean_lag), size = 10) +
    geom_errorbar(data = meanVals, aes(x=lag, y = mean_lag, 
                                       ymin = mean_lag - ste_lag,
                                       ymax = mean_lag + ste_lag),
                  width = .5, linewidth = 2) + 
    ylab('T2 detection p(T2|T1)') 
  # rect(1, 5, 3, 7, col="white")
  dev.off()
  
  #linear mixed effects model 
  anova_result <- lmer(value ~ lag + (1 | subID), data = plotDat)
  anova(anova_result)
  

#I want the ending stim presentation time
df$stimTim = unlist(lapply(allDat, function(x) x$curStimTim[240]))
print(paste0('the mean ending stim presentation time was: ', round(mean(df$stimTim),2)))
#hard code to average over different T1 measures at different lags
df$allT1 = rowMeans(df[,seq(6,26,3)])


####    main task design Table                   #### 
mainDesign = data.frame('Lag' = c(1,5), 
                        'Both Targets' = c('84 (84)', '84 (50)'),
                        'T1 Only' = c('28 (0)', '28 (0)'),
                        'T2 Only' = c('28', '28'),
                        'No Target' = c('10', '10'))

mainDesign %>% 
  kbl(align = 'c') %>% 
  kable_classic(full_width = F, 
                font_size = 20) %>%
  column_spec(1, border_right = T) %>%
  row_spec(1, align = 'c')%>%
  footnote(general = "Number of trials presented during the AB task. Parenthetical numbers indicate the subset of T1 items subsequently presented during the memory test.",
           general_title = "Table 1: ",
           footnote_as_chunk = T, title_format = c("italic", "underline")
  )  
                        
                        

#### Experiment 2 plotting and final analysis #################################
#NOTE: summary variables were extracted using python 
#see notebook: ABmemNotebook.ipynb 
#also see helper python functions: AB_mem_helper_funcs.py
#bring in the data for the new version of the experiment: 
df = read.csv(paste0(repoPath, '/Exp2Summary.csv'))
#implement criterion 1:
df = df[df$ABT1overall>.5, ]
#implement criterion 2:
df = df[df$TargALL_UVSD_d >.1, ]
#implement criterion 3:
df = df[df$Lag5_T1P_T2P_T2>.5, ]

sum(df$DistALL_UVSD_x>18.3)
sum(df$Dist_UVSD_x>31.4)

ABmainTable = data.frame('Targets' = c('Both', 'Both', 'T1', 'T1', 'T2', 'T2', 'Neither', 'Neither'), 
                         'Lag'     = rep(c(1,5), 4), 
                         'T1'      = rep(".", 8), 
                         'T2'      = rep(".", 8), 
                         'T2_T1'   = rep(".", 8))
lags = c('Lag1', 'Lag5')
presAbs = c('P', 'A')

ti = 1
  for(t1p in 1:2){
    for(t2p in 1:2){
      for(L in 1:2){
      T1Vals = df[[paste(lags[L],'_T1',presAbs[t1p],'_T2',presAbs[t2p],'_T1', sep="")]]
      T2Vals = df[[paste(lags[L],'_T1',presAbs[t1p],'_T2',presAbs[t2p],'_T2', sep="")]]
      T2_T1Vals = df[[paste(lags[L],'_T1',presAbs[t1p],'_T2',presAbs[t2p],'_T2_T1', sep="")]]
      ABmainTable$T1[ti] = paste(as.character(round(mean(T1Vals), 2)), ' (',
                                 as.character(round(sd(T1Vals), 2)), ')', sep="")
      ABmainTable$T2[ti] = paste(as.character(round(mean(T2Vals), 2)), ' (',
                                 as.character(round(sd(T2Vals), 2)), ')', sep="")
      ABmainTable$T2_T1[ti] = paste(as.character(round(mean(T2_T1Vals), 2)), ' (',
                                 as.character(round(sd(T2_T1Vals), 2)), ')', sep="")
      ti = ti+1
    }
  }
  }
ABmainTable %>% 
  kbl(align = 'c') %>% 
  kable_classic(full_width = F, 
                font_size = 20) %>%
  column_spec(2, border_right = T) %>%
  row_spec(1, align = 'c')%>%
  footnote(general = "Accuracy and (standard deviation) for different trial types for T1, T2, and T2|T1 ",
           general_title = "Table 2: ",
           footnote_as_chunk = T, title_format = c("italic", "underline")
  )

# 
# T2T1.plotDat = df[,c("X", "Lag1_T1P_T2P_T2_T1", "Lag5_T1P_T2P_T2_T1")] %>%
#   pivot_longer(cols = c("Lag1_T1P_T2P_T2_T1", "Lag5_T1P_T2P_T2_T1"),
#                names_to = "lag", 
#                values_to = "value")
# 
# 
# T2T1.plotDat %>% ggplot(aes(x=lag, y=value)) +
#   geom_boxplot() +
#   geom_line(aes(group = X), position = position_dodge(.2)) + 
#   geom_point(aes(group = X), position = position_dodge(.2))

#### plot T1 detection accuracy for Exp 2 ####################################

plotDat <- df[,c("Lag1_T1P_T2P_T1", "Lag5_T1P_T2P_T1")] 
plotDat$subID <- seq(1,length(plotDat[,1]), 1)
plotDat <- plotDat %>%
  pivot_longer(cols = starts_with('lag'),
               names_to = "lag", 
               values_to = "value")

meanVals <- plotDat %>% 
  group_by(lag) %>% 
  summarise(mean_lag = mean(value), 
            ste_lag = sd(value)/sqrt(46)* qt(.975, 46))
meanVals$subID = 99


xVals = unique(plotDat$lag)
png("main_T1.pdf",         # File name
    width=500, height=500)
plotDat %>% ggplot(aes(x=lag, y=value, group = subID)) +
  geom_line( position = position_dodge(.2), size = 1, alpha = .15) + 
  geom_point( position = position_dodge(.2)) +
  scale_x_discrete(labels = c('1','5')) +
  scale_y_continuous(breaks = seq(.4, 1.0, .1), limits = c(.4,1)) +
  theme_classic() +
  theme(axis.line = element_line(color = 'black', size = 3),
        axis.ticks = element_line(colour = "black", size = 2),
        axis.ticks.length=unit(-.25, "cm"),
        text = element_text(size = 20)) + 
  geom_line(data = meanVals, aes(x= lag, y = mean_lag), linewidth = 4) + 
  geom_point(data = meanVals, aes(x= lag, y = mean_lag), size = 10) +
  geom_errorbar(data = meanVals, aes(x=lag, y = mean_lag, 
                                     ymin = mean_lag - ste_lag,
                                     ymax = mean_lag + ste_lag),
                width = .5, linewidth = 2) + 
  ylab('T1 detection p(T1)')
# rect(1, 5, 3, 7, col="white")
dev.off()


t.test(df$Lag1_T1P_T2P_T1 - df$Lag5_T1P_T2P_T1)

#### plot T2|T1 accuracy to show main AB effect ####################################

plotDat <- df[,c("Lag1_T1P_T2P_T2_T1", "Lag5_T1P_T2P_T2_T1")] 
plotDat$subID <- seq(1,length(plotDat[,1]), 1)
plotDat <- plotDat %>%
  pivot_longer(cols = starts_with('lag'),
               names_to = "lag", 
               values_to = "value")

meanVals <- plotDat %>% 
  group_by(lag) %>% 
  summarise(mean_lag = mean(value), 
            ste_lag = sd(value)/sqrt(46)* qt(.975, 46))
meanVals$subID = 99


xVals = unique(plotDat$lag)
png("main_T2.pdf",         # File name
    width=500, height=500)
plotDat %>% ggplot(aes(x=lag, y=value, group = subID)) +
  geom_line( position = position_dodge(.2), size = 1, alpha = .15) + 
  geom_point( position = position_dodge(.2)) +
  scale_x_discrete(labels = c('1','5')) +
  scale_y_continuous(breaks = seq(.2, 1.0, .1), limits = c(.2,1)) +
  theme_classic() +
  theme(axis.line = element_line(color = 'black', size = 3),
        axis.ticks = element_line(colour = "black", size = 2),
        axis.ticks.length=unit(-.25, "cm"),
        text = element_text(size = 20)) + 
  geom_line(data = meanVals, aes(x= lag, y = mean_lag), linewidth = 4) + 
  geom_point(data = meanVals, aes(x= lag, y = mean_lag), size = 10) +
  geom_errorbar(data = meanVals, aes(x=lag, y = mean_lag, 
                                     ymin = mean_lag - ste_lag,
                                     ymax = mean_lag + ste_lag),
                width = .5, linewidth = 2) + 
  ylab('T2 detection p(T2|T1)')
# rect(1, 5, 3, 7, col="white")
dev.off()

t.test(df$Lag1_T1P_T2P_T2_T1 - df$Lag5_T1P_T2P_T2_T1)



#### memory for targets versus distractors #######################################################
plotDat <- df[,c("DistALL_UVSD_d", "TargALL_UVSD_d")] 
plotDat$subID <- seq(1,length(plotDat[,1]), 1)
plotDat <- plotDat %>%
  pivot_longer(cols = ends_with('_d'),
               names_to = "type", 
               values_to = "value")

meanVals <- plotDat %>% 
  group_by(type) %>% 
  summarise(mean_lag = mean(value), 
            ste_lag = sd(value)/sqrt(46) * qt(.975, 46))
meanVals$subID = 99


png("main_distTargMem.pdf",         # File name
    width=500, height=500)
plotDat %>% ggplot(aes(x=type, y=value, group = subID)) +
  geom_line( position = position_dodge(.2), size = 1, alpha = .15) + 
  geom_point( position = position_dodge(.2)) +
  scale_x_discrete(labels = c('AB distractors','AB targets')) +
  scale_y_continuous(breaks = seq(0, 3.0, .5), limits = c(0,3)) +
  theme_classic() +
  theme(axis.line = element_line(color = 'black', size = 3),
        axis.ticks = element_line(colour = "black", size = 2),
        axis.ticks.length=unit(-.25, "cm"),
        text = element_text(size = 20)) + 
  # geom_line(data = meanVals, aes(x= type, y = mean_lag), linewidth = 4) + 
  geom_point(data = meanVals, aes(x= type, y = mean_lag), size = 10) +
  geom_errorbar(data = meanVals, aes(x=type, y = mean_lag, 
                                     ymin = mean_lag - ste_lag,
                                     ymax = mean_lag + ste_lag),
                width = .5, linewidth = 2) + 
  ylab("memory performance (d')")
# rect(1, 5, 3, 7, col="white")
dev.off()

#evaluating output from 3 distribution model
t.test(df$DistALL_UVSD_d - df$TargALL_UVSD_d)

######






#### Compare memory for blink, non-blink, and lag 5 targets #######################################################
plotDat <- df[,c("L1blink_UVSD_d", "L1cor_UVSD_d", "L5_UVSD_d")] 
plotDat$subID <- seq(1,length(plotDat[,1]), 1)
plotDat <- plotDat %>%
  pivot_longer(cols = ends_with('_d'),
               names_to = "type", 
               values_to = "value")

meanVals <- plotDat %>% 
  group_by(type) %>% 
  summarise(mean_lag = mean(value), 
            ste_lag = sd(value)/sqrt(46) * qt(.975, 46))
meanVals$subID = 99


png("main_blinkCorMem.pdf",         # File name
    width=500, height=500)
plotDat %>% ggplot(aes(x=type, y=value, group = subID)) +
  geom_line( position = position_dodge(.2), size = 1, alpha = .15) + 
  geom_point( position = position_dodge(.2)) +
  scale_x_discrete(labels = c('AB blink','AB no blink', 'lag 5')) +
  scale_y_continuous(breaks = seq(0, 3.0, .5), limits = c(0,3)) +
  theme_classic() +
  theme(axis.line = element_line(color = 'black', size = 3),
        axis.ticks = element_line(colour = "black", size = 2),
        axis.ticks.length=unit(-.25, "cm"),
        text = element_text(size = 20)) + 
  # geom_line(data = meanVals, aes(x= type, y = mean_lag), linewidth = 4) + 
  geom_point(data = meanVals, aes(x= type, y = mean_lag), size = 10) +
  geom_errorbar(data = meanVals, aes(x=type, y = mean_lag, 
                                     ymin = mean_lag - ste_lag,
                                     ymax = mean_lag + ste_lag),
                width = .5, linewidth = 2) + 
  ylab("memory performance (d')")
# rect(1, 5, 3, 7, col="white")
dev.off()

#apply linear mixed effects model: 
memory.lm = lmer(value~type+(1|subID), data = plotDat)
anova(memory.lm)

#pairwise tests
t.test(df$L1blink_UVSD_d - df$L1cor_UVSD_d)

mean(df$L1blink_UVSD_d - df$L1cor_UVSD_d) / sd(df$L1blink_UVSD_d - df$L1cor_UVSD_d)

t.test(df$L1cor_UVSD_d - df$L5_UVSD_d)
t.test(df$L1blink_UVSD_d - df$L5_UVSD_d)

#does memory predict the AB magnitude? looking at it in several ways, 
#and none indicate an effect
(df$L1cor_UVSD_d + df$L1blink_UVSD_d) / 2 -> df$L1_d_mean
df %>% ggplot(aes(x = L1_d_mean, y = blinkMag)) + geom_point()

(df$L1cor_UVSD_d + df$L1blink_UVSD_d) / 2 -> df$L1_d_mean
df$lagMem = df$L1_d_mean - df$L5_UVSD_d
df %>% ggplot(aes(x = lagMem, y = blinkMag)) + geom_point()


cor(df$blinkMag, df$TargALL_UVSD_d)
df$blinkMemDiff = df$L1blink_UVSD_d - df$L1cor_UVSD_d
df %>% ggplot(aes(x = blinkMemDiff, y= blinkMag)) +
  geom_point()

#### overall table of summary values for memory performance ####

memoryTable = data.frame('targetType' = c('target', 'T1|blink', 'T1|T2', 'lag5', 'dist', 'novel'),
                         'accuracy' = rep('', 6), 
                         'd'          = rep('', 6), 
                         'SDratio'   = rep('', 6),
                         'con1'      = rep('', 6),
                         'con2'      = rep('', 6),
                         'con3'      = rep('', 6),
                         'con4'      = rep('', 6),
                         'con5'      = rep('', 6),
                         'con6'      = rep('', 6),
                         'trialCount'= rep('', 6))
typesPer = c('targMem', 'PP1blink_mem', 'PP1cor_mem', 'PP5cor_mem', 'distMem', 'novel_mem')
typesUVSD = c('TargALL_UVSD_', 'L1blink_UVSD_', 'L1cor_UVSD_', 'L5_UVSD_', 'Dist_UVSD_', 'nov')
confShort = c('TargALL', 'L1blink', 'L1cor', 'L5', 'DistALL', 'NovALL')
for(ti in 1:6) { 
  memoryTable$accuracy[ti] = paste(as.character(round(mean(df[[typesPer[ti]]]), 2)), ' (',
                                   as.character(round(sd(df[[typesPer[ti]]]), 2)), ')', sep="")
  confRatings = df[,apply(as.matrix(names(df)), 1, function(x) grepl(confShort[ti], x))]
  confRatings = confRatings[,1:6]
  counts = rowSums(confRatings)
  for(ci in 1:6){
    memoryTable[ti,4+ci] = paste(as.character(round(mean(confRatings[,ci]/counts),2 )), ' (',
          as.character(round(sd(confRatings[,ci]/counts),2 )), ')', sep="")
  }
  memoryTable[ti,11] = paste(as.character(round(mean(counts))), ' (',
                             as.character(round(sd(counts))), ')', sep="")
  if(ti<6) {
  memoryTable$d[ti] = paste(as.character(round(mean(df[[paste(typesUVSD[ti], 'd',sep = '')]]),2)), ' (',
                           as.character(round(sd(df[[paste(typesUVSD[ti], 'd',sep = '')]]),2)), ')', sep="")
  memoryTable$SDratio[ti] = paste(as.character(round(mean(df[[paste(typesUVSD[ti], 's',sep = '')]]),2)), ' (',
                            as.character(round(sd(df[[paste(typesUVSD[ti], 's',sep = '')]]),2)), ')', sep="")
  }
}



memoryTable %>% 
  kbl(align = 'c') %>% 
  kable_classic(full_width = F, 
                font_size = 20) %>%
  column_spec(1, border_right = T) %>%
  row_spec(1, align = 'c')%>%
  footnote(general = "Memory accuracy, UVSD model fit, and confidence rating values separated by AB task trial type. 
                      target = memory for items that served as targets during AB task
                      T1|blink = memory for T1 items that were followed by T2 misses at lag 1
                      T1|T2 = memory for T1 items that were followed by T2 hits at lag 1
                      dist = memory for items that served as distractors during AB task
                      novel = correct rejection rate for novel items presented only during the memory test
                      mean (standard deviation)",
           general_title = "Table 3: ",
           footnote_as_chunk = T, title_format = c("italic", "underline")
  )


#### Stimulus characteristic question ###########################
#Was there something different about the T1 images on blink trials versus non-blink trials? 
#Specifically, if low familiarity images require more time to process, then 
#they might be more likely to result in blinks. 

#There is a significant difference in familairity of T1 images from blink trials
#versus non-blink trials, but the mean difference is small
t.test(df$blinkFam - df$nonBlinkFam)

mean(df$blinkFam - df$nonBlinkFam) / sd(df$blinkFam - df$nonBlinkFam)

mean(df$blinkFam)
mean(df$nonBlinkFam)
df$blinkMemEffect = df$L1blink_UVSD_d - df$L1cor_UVSD_d
df$blinkFamEffect =  df$nonBlinkFam - df$blinkFam 
df %>% ggplot(aes(x = blinkFamEffect, y = blinkMemEffect)) + 
  geom_point()
cor.test(df$blinkFamEffect, df$blinkMemEffect)

#individuals more affected by familiarity of T1 images were not more likely to 
#have bigger blink memory effects



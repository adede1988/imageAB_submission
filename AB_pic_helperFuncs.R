#helper functions for analysis of cognitive challenge task
#Adam Dede; Feb, 2021
#adam.osman.dede@gmail.com
#This task has three components with distinct trial types in each:
#1) picture task
#2) letter task
#3) affective rating task

#assign colors to the two groups 
groupCols <<- c("black", "#C66222") #useful page for finding colors: https://htmlcolorcodes.com/

twoDistBarPlot <- function(plotvals, SDvals) {
  par(mar = c(5,4,4,5))
  opar <- par(lwd = 20)
  barplot(plotvals, beside = TRUE, ylim=c(1, 1.8), col = c('white', 'white'), 
          cex.axis = 1.0, border = 'black', 
          axes = FALSE, cex.names = 1.25, xpd = FALSE) -> twoDistPlot
  
  axis(1,lwd.ticks = 0, labels=FALSE, lwd = 15)
  axis(1, at=c(-5,1), labels = F, pos = 1, lwd = 15, tck = 0)
  axis(2, tck =.03, pos = 0, at=seq(1,1.8,.1), labels= seq(1,1.8,.1), lwd = 15)
  # text(x = 0, y = c(0, .5, 1, 1.5), par('usr')[1], labels=c(0, .5, 1, 1.5), srt = 0, pos = 2, xpd = T)
  #error bars
  arrows(twoDistPlot, plotvals+SDvals, 
         twoDistPlot, plotvals, angle = 90, code = 1, length = .10, lwd = 10)
  arrows(twoDistPlot, plotvals-SDvals, 
         twoDistPlot, plotvals, angle = 90, code = 1, length = .10, lwd = 10)
}



readCCdat = function(filePath) {
  #set the program to open the specified folder
  setwd(filePath)
  #list out character vectors of all the files in the specified folder
  datFiles = list.files()  
  #import all the data into a list of dataframes
  allData <- lapply(datFiles, read.csv)
  #initialize a summary dataframe to be carried around through all analyses 
  n = length(allData)
  df <- data.frame('year' = replicate(n,0),
                   'month' = replicate(n,0),
                   'day' = replicate(n,0),
                   'hour' = replicate(n,0),
                   'min' = replicate(n,0),
                   stringsAsFactors=FALSE)
  df[,] = t(matrix(unlist(lapply(allData, function(x) dateSplit(x$date[1]))), nrow = 5))
  return(list(allData, df))
} 


plotBasicAB <- function(allDat,df, statOffset = 3){
  #statOffset input allows for flexibility in plotting: 
  #1 = T1
  #2 = T2
  #3 = T2|T1
  
  #For each trial type there are a few things I want: 
  goalStats = c('T1', 'T2', 'T2.T1')
  
  #all this is for augmenting the summary data frame: 
  #only need to do this once, so check if it's already been done
  if(sum(grepl("lag", names(df))==T)==0){
    abDat = lapply(allDat, function(x) subset(x, phase=="main" & skip1==0 & skip2 ==0))
    #get names of all different trial types
    abTypes = sort(unique(abDat[[1]]$cndTyp))
    #combine the trial types and goal stats
    colsOut = unlist(lapply(abTypes, function(x) unlist(lapply(goalStats, function(y) paste(x,y, sep='.')))))
    #make a data.frame with each row representing a participant and columns with colsOut as names
    abSummary = unlist(lapply(abDat, function(x) c(unlist(lapply(abTypes, function(y) getTrialTypeSum(x, y, goalStats))))))
    abSummary = setNames(data.frame(t(matrix(abSummary, ncol = length(abDat), nrow = length(colsOut)))), colsOut)
    #add the new rows to the overall dataframe 
    df = cbind(df, abSummary)
  }
  #now, I want to plot summary statistic over lags 1-7
  #subset down to the relevant data
  lagDat = df[, grepl("lag" , names(df))]
  lagDat = lagDat[, seq(1+statOffset-1,21,3)] #hard coded method for accessing T1, T2, or T2|T1 data
  #get means
  lagDatMeans = colMeans(lagDat)
  print(lagDatMeans)
  #get confidence intervals
  lagDatCI = apply(lagDat, 2, conInt)
  #layout is 1X1 (single plot) margins are 4=bottom then counter clockwise
  # par(mfrow=c(1,1), mar=c(4,6,2,4))
  # par(xpd=TRUE) #allow drawing outside plot bounds
  #initialize the plot
  plot(c(0), axes = F, ylab = goalStats[statOffset], xlab = "lag", cex.lab = 1, ylim=c(.5, 1), xlim=c(.8,7.05))
  #lines for group 1
  # for(ii in 1:nrow(df)){
  #   myLines(lagDat[ii,], replicate(7,0), 2)
  #   
  # }
  myLines(lagDatMeans,lagDatCI,1)
  
  #axes
  axis(1, at=c(1:7), labels = c(1:7), pos = .5, lwd=15,  cex.axis = 1, tck = .03)
  axis(1, at=c(-1,1), labels = F, pos = .5, lwd = 15, tck = 0)
  axis(2, at=seq(.5,1,.1), labels= seq(.5,1,.1), pos = .5, lwd = 15, cex.axis = 1, tck = .03)
  
  #add a legend
  # legend(5, .5, legend=c("low/low", "high/low"),
  #        col=c(groupCols[1], groupCols[2]), cex=1.5, pch = 15)
  
  return(df)
}

#calculates a confidence interval size based on a t-value
#default confidence interval is .95
conInt <- function(x, confidence = .95) {
  x = x[!is.na(x)] #NA remove
  return ( (sd(x)/sqrt(length(x))))#*qt(1-(1-confidence)/2, length(x)-1)) #standard error * t critical
}

#plots lines 
#assumes integer separation in x between y values
myLines <- function(y,CI, group){
  if(group==1){
    x = seq(1, length(y))+(group-1)*.05 #x-coordinates
    lines(x, y, lwd = 10, col = groupCols[group]) #make line segments
    points(x, y, pch = 15, col = groupCols[group], cex = 10) #put scatter points on 
    arrows(x, y, y1=y+CI, angle = 90, lwd = 10, col = groupCols[group]) #up error bars
    arrows(x, y, y1=y-CI, angle = 90, lwd = 10, col = groupCols[group]) #down error bars
  } else {
    x = seq(1, length(y))#x-coordinates
    
    lines(x, y, lwd = 10, col = rgb(.9,.9,.9))
    # points(x, y, pch = 15, col = rgb(0,0,0, alpha = 0.5),
    #         cex = 3) #put scatter points on
  }
}

#pulls out the desired summary stats in the desired order for a particular letter task trial type
#NOTE: hard coding of the desired stats is necessary since they have different calculation needs
#But the given order is maintained. 
getTrialTypeSum<-function(df, trialTyp, stats){
  df = subset(df, df$cndTyp == trialTyp) #subset to only the current trial type
  df = subset(df, !is.na(df$cor1)) #get rid of accidental early press trials
  df$cor1 = as.numeric(df$cor1) #convert to numeric in case it was character vector
  res = array(rep(NA,length(stats))) #1d vector for storage of results
  for (ss in 1:length(stats)){ #maintain input order by going through in order
    if (stats[ss]=="T1"){
      res[ss] = mean(df$cor1)
    }
    if (stats[ss]=="T2"){
      res[ss] = mean(df$cor2)
    }
    if (stats[ss]=="T2.T1"){
      temp = subset(df, df$cor1==1)
      res[ss] = mean(temp$cor2)
    }
    
  }
  return(res)
}


#split date stamps out into year, month, day, hour, minute
dateSplit <- function(inDate){
  s1 = strsplit(inDate, '_')
  year = as.numeric(s1[[1]][1])
  month = s1[[1]][2]
  day = as.numeric(s1[[1]][3])
  hour = as.numeric(substring(s1[[1]][4], 1, 2))
  min = as.numeric(substring(s1[[1]][4], 3, 4))
  return (c(year,month,day,hour,min))
}

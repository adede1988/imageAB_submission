#helper functions for dual/single AB task 

readAB_blk <- function(curBlk, bossNorms, bossNorms2){
  #remove extra columns
  curBlk <- curBlk[,colSums(is.na(curBlk)) != length(curBlk[,1])]
  curBlk <- curBlk[,colSums(apply(curBlk, c(1,2), function(x) identical("", x)))
                   != length(curBlk[,1]) ]
  outDat = data.frame()
  #Trial Types:
  T1 = c(0, 1)
  T2 = c(0, 1)
  L = c(1,5)
  TTs = expand.grid(T1, T2, L)
  pres = c('T1', 'T2', 'L')
  names(TTs) <- pres
  curBlk$T1cor = 0
  curBlk$T2cor = 0
  curBlk$T1pnts = 0 
  curBlk$T2pnts = 0
  curBlk$T1Fam = 0
  curBlk$T1Fam = sapply(curBlk$T1_Identity, function(X){findT1Fam(X, bossNorms, bossNorms2)})
  #make a correct / incorrect column for T1 and T2 depending on trial type
  for(t in 1:dim(TTs)[1]){
    idx = which(curBlk$T1_PA==TTs$T1[t] & 
                  curBlk$T2_PA==TTs$T2[t] & 
                  curBlk$Lag==TTs$L[t])
    T1Rsp = curBlk$key_resp_2.keys[idx]
    if( sum(grepl("key_resp_3.keys", names(curBlk)))==1){
      T2Rsp = curBlk$key_resp_3.keys[idx]
      doT2 = T
    } else{
      doT2 = F
    }
    if(TTs$T1[t] == 1){
      curBlk$T1cor[idx] = as.numeric(T1Rsp == 'right')
      curBlk$T1pnts[idx] = curBlk$T1cor[idx] + as.numeric(T1Rsp == 'left') * -1
    } else{
      curBlk$T1cor[idx] = as.numeric(T1Rsp == 'left')
      curBlk$T1pnts[idx] = curBlk$T1cor[idx] + as.numeric(T1Rsp == 'right') * -3
    }
    if( doT2){
      if(TTs$T2[t] == 1){
        curBlk$T2cor[idx] = as.numeric(T2Rsp == 'right')
        curBlk$T2pnts[idx] = curBlk$T2cor[idx] + as.numeric(T2Rsp == 'left') * -1
      } else{
        curBlk$T2cor[idx] = as.numeric(T2Rsp == 'left')
        curBlk$T2pnts[idx] = curBlk$T2cor[idx] + as.numeric(T2Rsp == 'right') * -3
      }
      
      cutidx = curBlk[idx,] #block data cut to the current trial type
      outDat[t,'T2cor']=mean(curBlk$T2cor[idx])
      outDat[t,'T2T1cor'] = mean(cutidx$T2cor[cutidx$T1cor==1])   
      outDat[t, 'T2'] = TTs$T2[t] 
      
    } else {
      outDat[t,'T2cor']= NA
      outDat[t,'T2T1cor'] = NA  
      outDat[t, 'T2'] = TTs$T2[t] 
    }
    outDat[t,'T1cor']=mean(curBlk$T1cor[idx])
    outDat[t, 'T1'] = TTs$T1[t] 
    outDat[t, 'L']  = TTs$L[t] 
    outDat[t, 'trialCount'] = length(idx)
    
  }
  
  
  
  outDat$blinkT1Fam = mean(curBlk$T1Fam[curBlk$T1_PA==1 & curBlk$T2_PA==1 &
                                          curBlk$T1cor == 1 & curBlk$T2cor ==0], 
                           na.rm = T)
  outDat$nonBlinkT1Fam = mean(curBlk$T1Fam[curBlk$T1_PA==1 & curBlk$T2_PA==1 &
                                          curBlk$T1cor == 1 & curBlk$T2cor ==1], 
                           na.rm = T)
  
  outDat$pnts = sum(curBlk$T1pnts) + sum(curBlk$T2pnts)
  if(sum(grepl('PROLIFIC_PID', names(curBlk)))>0){
    outDat$subID = curBlk$PROLIFIC_PID[1]
  }else{
    outDat$subID = curBlk$subjectID[1]
  }
  outDat$date = curBlk$date[1]
  
  return(outDat)
}

readMEM_blk <- function(curBlk, ABdat){
  #remove extra columns
  curBlk <- curBlk[,colSums(is.na(curBlk)) != length(curBlk[,1])]
  curBlk <- curBlk[,colSums(apply(curBlk, c(1,2), function(x) identical("", x)))
                   != length(curBlk[,1]) ]
  #remove extra columns
  ABdat <- ABdat[,colSums(is.na(ABdat)) != length(ABdat[,1])]
  ABdat <- ABdat[,colSums(apply(ABdat, c(1,2), function(x) identical("", x)))
                 != length(ABdat[,1]) ]
  #detect if a single target or double target block: 
  if( sum(grepl("key_resp_3.keys", names(ABdat)))==1){
    doT2 = T
  } else {
    doT2 = F
  }
  outDat = data.frame('count' = rep(0,24), #trial count
                      'conf'  = rep(c(1,2,3,4,5,6),4), #confidence rating
                      'trialType' = c(rep('blink', 6), rep('noBlink',6),
                                      rep('lag5', 6), rep('novel', 6)))
  curBlk$pnts = 0
  #loop the mem trials and ID them in the AB task to determine trial type:
  #blink, no blink, lag 5, novel
  for(tt in 1:dim(curBlk)[1]){
    conf = curBlk$key_resp_4.keys[tt]
    link = which(ABdat$trialID == curBlk$ABtrialID[tt])
    if(length(link)>1){ #only possible in practice!
      link = link[length(link)]
    }
    if(length(link)>0) { #this is an item that was presented previously 
      if(curBlk$priorLag[tt] == 5) { #lag 5
        ABtrial = ABdat[link,]
        T1rsp = ABtrial$key_resp_2.keys
        # if('right' == T1rsp){
        outDat$count[12+conf] = outDat$count[12+conf] + 1
        # }
      } else { #LAG 1: blink or no blink trial
        ABtrial = ABdat[link,]
        T1rsp = ABtrial$key_resp_2.keys
        if(doT2){
          T2rsp = ABtrial$key_resp_3.keys
          if('right' == T1rsp & 'right' == T2rsp){ #no Blink
            outDat$count[6+conf] = outDat$count[6+conf] + 1
          } else if('right' == T1rsp & 'left' == T2rsp){ # blink
            outDat$count[conf] = outDat$count[conf] + 1
          }
        }else{ #all "Blink" since there's no T2, so they never see T2
          # if('right' == T1rsp){
            outDat$count[conf] = outDat$count[conf] + 1
          # }
        }
        
      }
      
    } else { #this item is novel
      outDat$count[18+conf] = outDat$count[18+conf] + 1
    }
    
    if(curBlk$CondType[tt] == "old"){
      if(conf>3){
        curBlk$pnts[tt] = conf - 3
      }else{
        curBlk$pnts[tt] = conf - 4
      }
    } else { 
      if(conf<4){
        curBlk$pnts[tt] = abs(conf - 4) 
      } else {
        if(conf==6){
          curBlk$pnts[tt] = -9
        }else if(conf == 5){
          curBlk$pnts[tt] = -4
        } else{
          curBlk$pnts[tt] = -1
        }
        
      }
    }
    
    
  }
  outDat$dualTarget = doT2
  outDat$pnts = sum(curBlk$pnts)
  if(sum(grepl('PROLIFIC_PID', names(curBlk)))>0){
    outDat$subID = curBlk$PROLIFIC_PID[1]
  }else{
    outDat$subID = curBlk$subjectID[1]
  }
  outDat$date = curBlk$date[1]
  outDat$trialCount = dim(curBlk)[1]
  return(outDat)
  
}



uvsdPred <- function(params){
  #params is a data frame with rows equal to number of trial types
  #cols: d', s, c1, c2, c3, c4, c5
  #hard code to keep criteria the same across trial types
  #trying to get predicted probabilities of each confidence rating for 
  #each trial type
  crits = c(-999, unname(unlist(params[1,3:dim(params)[2]])), 999)
  outPreds = pnorm(matrix(rep(crits,dim(params)[1]), nrow = dim(params)[1], byrow = T), params[,1], params[,2])
  outPreds = apply(outPreds, 1, diff)
  
  return(outPreds)
  
}

uvsdLike <- function(params, conf){
  #params is c1,c2,c3,c4,c5, followed by nX2 d/s params
  #n is number of trial type conditions
  
  #reorganize counts to be 6 (confidence) X (conditions)
  numTypes = length(unique(conf$trialType))
  types = unique(conf$trialType)
  numConf = length(unique(conf$conf))
  
  #set ordering of criteria: 
  if(any(diff(params)[1:(numConf-2)]<0)){
    return(Inf)
  }
  # #guard rail against negative memory strength
  if(any(params[seq(numConf,length(params),2)]<0 ) & #0
     any(params[seq(numConf,length(params),2)]>4 )){ #4
    return(Inf)
  }
  # #guard rail to limit range of S parameter
  if(any(params[seq(numConf+1,length(params),2)]>2 | #2
         params[seq(numConf+1,length(params),2)]<.5)){ #.5
    return(Inf)
  }
  
  params[(length(params)+1):(length(params)+2)] =c(0,1) #hard code novel items (always last!)
  #create a data frame of parameters with a separate row for each type of trial
  newParams = data.frame(d = params[seq(numConf,length(params),2)], 
                         s = params[seq(numConf+1,length(params),2)])
  for(cc in 1:(numConf-1)){
    newParams[,paste0('c',cc)] = params[cc]
  }
  
  
  
  
  predVals = uvsdPred(newParams)
  trialCounts = matrix(nrow = numConf, ncol = numTypes)
  for(tt in 1:length(types)){
    trialCounts[,tt] = conf$count[conf$trialType==types[tt]]
  }
  trialCounts[trialCounts==0] = 1/sum(trialCounts)
  if(any(is.nan(log(predVals)))){
    print(log(predVals))
    print(predVals)
    print(params)
    return(Inf)
  }
  LL = sum(log(predVals) * trialCounts)
  return(-LL)
  
}


fitUVSD <- function(conf){
  
  typeCount <- conf %>% group_by(trialType) %>% summarize(sum(count)) 
  
  numTypes = length(unique(conf$trialType))
  numConf = length(unique(conf$conf))
  types = unique(conf$trialType)
  for(tt in 1:length(types)){
    if(typeCount[typeCount$trialType == types[tt],2]==0){
      conf = conf %>% filter(trialType != types[tt])
      
    }
  }
  numTypes = length(unique(conf$trialType))
  numConf = length(unique(conf$conf))
  types = unique(conf$trialType)
  params = c(seq(-1,1, length.out = numConf-1), rep(c(1,1), numTypes-1)) 
  
  tmp = function(params){
    uvsdLike(params, conf)
  }
  
  fit <- optim(params, tmp)
  fit$par
  params = fit$par
  params[(length(params)+1):(length(params)+2)] =c(0,1) #hard code novel items
  newParams = data.frame(d = params[seq(numConf,length(params),2)], 
                         s = params[seq(numConf+1,length(params),2)])
  for(cc in 1:(numConf-1)){
    newParams[,paste0('c',cc)] = params[cc]
  }
  
  
  
  predVals = uvsdPred(newParams)
  nullPred = array(data = 1/numConf, dim = dim(predVals))
 
  trialCounts = matrix(nrow = numConf, ncol = numTypes)
  for(tt in 1:length(types)){
    trialCounts[,tt] = conf$count[conf$trialType==types[tt]]
  }
  trialCounts[trialCounts==0] = 1/sum(trialCounts)
  expected = predVals * matrix(rep(colSums(trialCounts),numConf), nrow = numConf, ncol = numTypes, byrow = T)
  #correct for huge deviations caused by very small values:
  expected[expected<1/sum(trialCounts)] = 1/sum(trialCounts)
  expectedNull = nullPred * matrix(rep(colSums(trialCounts),numConf), nrow = numConf, ncol = numTypes, byrow = T)
  #correct for huge deviations caused by very small values:
  expectedNull[expectedNull<1/sum(trialCounts)] = 1/sum(trialCounts)
  chiVal = sum((trialCounts - expected)^2 / expected)
  chiNull = sum((trialCounts - expectedNull)^2 / expectedNull)
  pval = 1 - pchisq(chiVal, (dim(expected)[1]-1)*(dim(expected)[2]-1))
  
  outStats = data.frame()
  outStats[1,'chi'] = chiVal
  outStats[1, 'chiP'] = pval
  outStats[1, 'varExp'] = 1 - (chiVal / chiNull)
  
  for(tt in 1:length(types)){
    outStats[1,paste0(types[tt], '_d')] = params[numConf+(tt-1)*2] 
    outStats[1,paste0(types[tt], '_s')] = params[(numConf+(tt-1)*2)+1] 
  }
  return(outStats)
}

simpConf <- function(conf, confSet){
  numTypes = length(unique(conf$trialType))
  numConf = length(unique(conf$conf))
  types = unique(conf$trialType)
  start = T
  if(confSet == 4){
    compressConf = list(1,c(2,3),c(4,5),6)
  } else if (confSet == 3){
    compressConf = list(1, c(2,3,4,5),6)
  } else if (confSet == 6){
    compressConf = list(1, 2,3,4,5,6)
  }
  for(tt in 1:numTypes){
    for(cc in 1:length(compressConf)){
      cur = conf[conf$trialType == types[tt] & conf$conf %in% compressConf[[cc]],]
      cur$count[1] = sum(cur$count)
      cur = cur[1,]
      cur$conf = cc
      if(start){
        outConf = cur
        start = F
      } else{
        outConf = rbind(outConf,cur)
      }
    }
  }
  return(outConf)
}

findT1Fam <- function(name, bossNorms, bossNorms2){
  idx = which(bossNorms['FILENAME'] == name)

if(any(is.nan(bossNorms$Fam_Mean[idx]))){
  idx = which(bossNorms2['FILENAME']==name)
  fam = bossNorms2$Fam_Mean[idx]
}else{
  fam = bossNorms$Fam_Mean[idx]
}
if(length(fam) == 1){
  fam = fam[1]
}else{
  fam = NaN
}
return(fam)
}



readABmemFile <- function(sub, numConf, bossNorms, bossNorms2){
  rawDat = read.csv(sub)
  blkNames = unique(rawDat$TrialType)
  blkNames = blkNames[-1]
  #skip analyzing practice blocks
  blkNames = blkNames[!grepl('prac', blkNames)]
  outDat = list()
  subjectSummary = data.frame()
  for(ii in 1:length(blkNames)){
    curBlk = rawDat %>% filter(TrialType == blkNames[ii])
    if(grepl("AB", blkNames[ii])){
      cbRes <- readAB_blk(curBlk, bossNorms, bossNorms2) 
      cbRes$blkName = blkNames[ii]
      #false alarm rates for T1 and T2 (1 - accuracy for target absent trials)
      T1FA =1 - sum(cbRes[cbRes$T1==0,'T1cor'] * 
                   cbRes[cbRes$T1==0,'trialCount']) / 
        sum(cbRes[cbRes$T1==0,'trialCount'])
      T2FA =1 - sum(cbRes[cbRes$T2==0,'T2cor'] * 
                   cbRes[cbRes$T2==0,'trialCount']) / 
        sum(cbRes[cbRes$T2==0,'trialCount'])
      
      #grab data for present present trials:
      ppRes = cbRes[cbRes$T1==1 & cbRes$T2==1,]
      subjectSummary[1,paste0(blkNames[ii], '_lag1_T2T1')] = ppRes$T2T1cor[ppRes$L==1]
      subjectSummary[1,paste0(blkNames[ii], '_lag1_T1')] = ppRes$T1cor[ppRes$L==1]
      subjectSummary[1,paste0(blkNames[ii], '_lag5_T2T1')] = ppRes$T2T1cor[ppRes$L==5]
      subjectSummary[1,paste0(blkNames[ii], '_lag5_T1')] = ppRes$T1cor[ppRes$L==5]
      subjectSummary[1,paste0(blkNames[ii], '_T1FA')] = T1FA
      subjectSummary[1,paste0(blkNames[ii], '_T2FA')] = T2FA
      subjectSummary[1,paste0(blkNames[ii], '_T1blinkFam')] = cbRes$blinkT1Fam[1]
      subjectSummary[1,paste0(blkNames[ii], '_T1nonblinkFam')] = cbRes$nonBlinkT1Fam[1]
      
      outDat[[ii]] = cbRes
    } else { 
      ABdat <- paste(sub("_.*", "", blkNames[ii]), "AB", sep = "_")
      ABdat = rawDat%>% filter(TrialType == ABdat)
      conf <- readMEM_blk(curBlk, ABdat)
      
      #blockwise accuracy (more stable than fitting SDT model to small dataset)
      #blink: 
      tmpConf = conf %>% filter(trialType == 'blink')
      subjectSummary[1,paste0(blkNames[ii], '_blink_acc')] = 
                                  sum(tmpConf$count[4:6]) / sum(tmpConf$count)
      #noBlink: 
      tmpConf = conf %>% filter(trialType == 'noBlink')
      subjectSummary[1,paste0(blkNames[ii], '_noBlink_acc')] = 
        sum(tmpConf$count[4:6]) / sum(tmpConf$count)
      #lag 5: 
      tmpConf = conf %>% filter(trialType == 'lag5')
      subjectSummary[1,paste0(blkNames[ii], '_lag5_acc')] = 
        sum(tmpConf$count[4:6]) / sum(tmpConf$count)
      #lag 1: 
      tmpConf = conf %>% filter(trialType == 'blink' | trialType == 'noBlink') %>%
        group_by(conf) %>% summarize(count = sum(count))
      subjectSummary[1,paste0(blkNames[ii], '_lag1_acc')] = 
        sum(tmpConf$count[4:6]) / sum(tmpConf$count)
      
      
      conf = simpConf(conf, numConf)
      
      
      outDat[[ii]] <- conf
      outDat[[ii]]$blkName = blkNames[ii]
      
      
    }
    subjectSummary[1,paste0(blkNames[ii], '_pnts')] = outDat[[ii]]$pnts[1]
  }
  subjectSummary[1, 'ID'] = outDat[[1]]$subID[1]
  subjectSummary[1, 'date'] = outDat[[1]]$date[1]
  subjectSummary[1, 'age'] = rawDat$age[1]
  subjectSummary[1, 'sex'] = rawDat$sex[1]
  subjectSummary[1, 'nativeLang'] = rawDat$native.language[1]
  #what final info do I want? 
  return(list(subjectSummary, outDat))
  
  
  
}














# 
# 
# 
# uvsdPredSimp <- function(params){
#   #params is a data frame with rows equal to number of trial types
#   #cols: d', s, c1, c2, c3,
#   #hard code to keep criteria the same across trial types
#   #trying to get predicted probabilities of each confidence rating for
#   #each trial type
#   
#   outPreds = matrix(nrow = 4, ncol = dim(params)[1])
#   
#   #con = 4:
#   
#   outPreds[4,] = 1-pnorm(params$c3, mean = params$d, sd = params$s)
#   outPreds[3,] = 1-pnorm(params$c2, mean = params$d, sd = params$s)-colSums(outPreds, na.rm = T)
#   outPreds[2,] = 1-pnorm(params$c1, mean = params$d, sd = params$s)-colSums(outPreds, na.rm = T)
#   outPreds[1,] = 1-colSums(outPreds, na.rm = T)
#   
#   return(outPreds)
#   
# }
# 
# uvsdLikeSimp <- function(params, conf){
#   #params is c1,c2,c3, followed by nX2 d/s params
#   #n is number of trial type conditions
#   
#   #reorganize counts to be 6 (confidence) X (conditions)
#   numTypes = length(unique(conf$trialType))
#   types = unique(conf$trialType)
#   
#   #add in hard coded sd ratio values:
#   tmp = params[1:3]
#   for(ii in 1:(numTypes-1)){
#     tmp[3 + (ii-1)*2 + 1] = params[ii+3]
#     tmp[3 + (ii-1)*2 + 2] = 1
#   }
#   params = tmp
#   #set ordering of criteria:
#   if(any(diff(params)[1:2]<0)){
#     
#     return(Inf)
#   }
#   
#   params[(length(params)+1):(length(params)+2)] =c(0,1) #hard code novel items
#   newParams = data.frame(d = params[seq(4,length(params),2)],
#                          s = params[seq(5,length(params),2)],
#                          c1 = rep(params[1], numTypes),
#                          c2 = rep(params[2], numTypes),
#                          c3 = rep(params[3], numTypes))
#   
#   
#   
#   predVals = uvsdPredSimp(newParams)
#   trialCounts = matrix(nrow = 4, ncol = numTypes)
#   for(tt in 1:length(types)){
#     trialCounts[,tt] = conf$count[conf$trialType==types[tt]]
#   }
#   trialCounts[trialCounts==0] = 1/sum(trialCounts)
#   if(any(is.nan(log(predVals)))){
#     print(log(predVals))
#     print(predVals)
#     print(params)
#     return(Inf)
#   }
#   LL = sum(log(predVals) * trialCounts)
#   return(-LL)
#   
# }
# 
# 
# fitUVSDSimp <- function(conf){
#   numTypes = length(unique(conf$trialType))
#   types = unique(conf$trialType)
#   params = c(0, .3, .6, rep(1, numTypes-1))
#   
#   tmp = function(params){
#     uvsdLikeSimp(params, conf)
#   }
#   
#   fit <- optim(params, tmp)
#   params = fit$par
#   
#   
#   tmp = params[1:3]
#   for(ii in 1:(numTypes-1)){
#     tmp[3 + (ii-1)*2 + 1] = params[ii+3]
#     tmp[3 + (ii-1)*2 + 2] = 1
#   }
#   params = tmp
#   params[(length(params)+1):(length(params)+2)] =c(0,1) #hard code novel items
#   newParams = data.frame(d = params[seq(4,length(params),2)],
#                          s = params[seq(5,length(params),2)],
#                          c1 = rep(params[1], numTypes),
#                          c2 = rep(params[2], numTypes),
#                          c3 = rep(params[3], numTypes))
#   
#   
#   
#   predVals = uvsdPredSimp(newParams)
#   trialCounts = matrix(nrow = 4, ncol = numTypes)
#   for(tt in 1:length(types)){
#     trialCounts[,tt] = conf[conf$trialType==types[tt],1]
#   }
#   trialCounts[trialCounts==0] = 1/sum(trialCounts)
#   expected = predVals * matrix(rep(colSums(trialCounts),4), nrow = 4, ncol = numTypes, byrow = T)
#   # chiStats = chisq.test(c(trialCounts), c(expected))
#   
#   chiVal = sum((trialCounts - expected)^2 / expected)
#   chiP = 1 - pchisq(chiVal, (4*numTypes) - 1 )
#   
#   outStats = data.frame()
#   outStats[1,'chi'] = chiVal
#   outStats[1, 'chiP'] = chiP
#   print(p)
#   for(tt in 1:length(types)){
#     outStats[1,paste0(types[tt], '_d')] = params[4+(tt-1)*2]
#   }
#   return(outStats)
# }


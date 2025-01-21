# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 09:13:33 2022

@author: Adam Dede
"""

import numpy as np
import pandas as pd
import math
from scipy.stats import norm
from scipy.optimize import minimize
import copy
from functools import partial

#path = r"C:\Users\dtf8829\Documents\ABmem\data\PPT15_ABstudy_2022-06-30_09h08.37.319.csv"


def findT1(trialDat):
    pics = trialDat.loc[trialDat.index.str.contains('Img')]
    return(pics.iloc[int(trialDat['T1Loc'])])

def refMem2AB(trialDat, ABDat): 
    trialLink = pd.Series(np.zeros(3), index = ['ABT1cor', 'ABT2cor', 'lag'])
    try: 
        ii = np.where((ABDat['T1imageID'].values == trialDat['Picture']) & (ABDat['T1_PA'].values == 'P') )[0][0]
        trialLink['ABT1cor'] = ABDat.iloc[ii,:]['T1cor']
        trialLink['ABT2cor'] = ABDat.iloc[ii,:]['T2cor']
        trialLink['lag'] = ABDat.iloc[ii,:]['Lag']
        return(trialLink)
    except:
        trialLink['ABT1cor'] = np.nan
        trialLink['ABT2cor'] = np.nan
        trialLink['lag'] = np.nan
        return(trialLink)


def UVSDlike(memDat, c, d, s, trialTypes, returnVal):
    #main function for dealing with the confidence data and generating 
    #predicted values based on model parameters
    
    if any(np.diff(c)<0): #guardrail to ensure monotonic increasing criteria
        return(99999999999)
    #how many of each trial were there
    typeCounts = np.zeros(len(trialTypes))
    #what were the empirical confidence ratings
    conRatEmp = np.zeros([len(trialTypes), 6])
    #what is this model's prediction of confidence ratings
    conRatPre = np.zeros([len(trialTypes), 6])
    for ii,key in enumerate(trialTypes):
        
        #go through all the possible subtypes contributing to this trial type: 
        for val in trialTypes[key]: 
            typeCounts[ii] += sum(memDat['ABrefType'] == val)
            #get counts of empirical confidence ratings: 
            conRatEmp[ii,:] += [sum(memDat.loc[memDat['ABrefType']==val, 'MemResponse.keys']==resp) for resp in range(1,7)]
        if ii==len(d): #this is the novel distribution! 
            conRatPre[ii,:] = np.diff(np.append(np.append(0, norm.cdf(c, loc=0, scale=1)),1))
        else:  
            conRatPre[ii,:] = np.diff(np.append(np.append(0, norm.cdf(c, loc=d[ii], scale=s[ii])),1)) 
    
    returnConRatEmp = copy.copy(conRatEmp)
    #bump everything to avoid zeros
    conRatEmp += 1 
    conRatPre[conRatPre==0] = 1/sum(sum(conRatEmp)) 
    typeCounts[typeCounts==0] = 1
    
    #get expected values for calculating chiSquared
    conRatExp = conRatPre * typeCounts.reshape(-1,1) 
    
    #prepare variables for calculating the summary stats
    likelihood = 0
    conRatEmp = conRatEmp.flatten()
    conRatPre = conRatPre.flatten()
    conRatExp = conRatExp.flatten()
    
    #look at each response type and add it into the likelihood
    for val, reps in zip(conRatPre, conRatEmp):
        #add the log of its probability for each time it was observed
        likelihood += np.sum(np.repeat(np.log(val), reps))
    
    chiSqVal = sum( (conRatEmp - conRatExp)**2 / conRatExp )
    if returnVal==1: 
        return(-likelihood)
    elif returnVal==2: 
        return(chiSqVal)
    else: 
        return(returnConRatEmp)

def UVSDwrap(memDat, c, d, s, trialTypes): 
    #wrapper function to condense the parameters into a single set 
    #and run the minimization 
    curFunc = lambda x: UVSDlike(memDat, 
                                 x[0:len(c)], 
                                 x[len(c):len(c)+len(d)], 
                                 x[len(c)+len(d):], 
                                 trialTypes, 1)
    curFunc2 = lambda x: UVSDlike(memDat, 
                                  x[0:len(c)], 
                                  x[len(c):len(c)+len(d)], 
                                  x[len(c)+len(d):], 
                                  trialTypes, 2)
    curFunc3 = lambda x: UVSDlike(memDat, 
                                  x[0:len(c)], 
                                  x[len(c):len(c)+len(d)], 
                                  x[len(c)+len(d):], 
                                  trialTypes, 3)
    
    params = np.concatenate([c,d,s])
    paramsOut = minimize(curFunc, params)['x']
    chiSqVal = curFunc2(paramsOut)
    confVals = curFunc3(paramsOut)
    return([paramsOut, chiSqVal, confVals])
    
    

def fitUVSD(memDat, out): 
    #fitting the UVSD model with 2 d values (one for targets, one for distracters)
    c = np.linspace(-1,2,5) #initialize criteria
    d = [2, .25] #initialize d 
    s = [1.5, 1] #initialize standard deviation values
    #make a dictionary of trial types, last is always novel 
    trialTypes = dict(zip(('TargALL_UVSD', 'DistALL_UVSD', 'NovALL_UVSD'), 
                          (['PP1cor', 'PP1blink', 'skip', 'PP1miss', 'PP1T2', 'PP5cor'], ['distracter'], ['novel'])))
    (params, chiSqVal, confRat) = UVSDwrap(memDat, c, d, s, trialTypes)
    for ii in range(len(trialTypes)-1): 
        out[list(trialTypes.keys())[ii] + '_d'] = params[5+ii]
        out[list(trialTypes.keys())[ii] + '_s'] = params[5+len(d)+ii]
        
    for ii, key in enumerate(trialTypes): 
        for conRating in list(range(6)): 
            out[key[0:-4] + str(conRating+1)] = confRat[ii,conRating]
    
    
    #fit it again, but this time with 4 d values (split targets into blink/cor lag1 and lag5 overall)
    c = np.linspace(-1,2,5) #initialize criteria
    d = [2, 2, 2, 1] #initialize d 
    s = [1.2,1.2,1.2,1] #initialize s
    #make a dictionary of trial types
    trialTypes = dict(zip(('L1cor_UVSD', 'L1blink_UVSD', 'L5_UVSD', 'Dist_UVSD', 'Nov_UVSD'), 
                          (['PP1cor'], ['PP1blink'], ['skip', 'PP5cor'], ['distracter'], ['novel'])))
    (params, chiSqVal, confRat) = UVSDwrap(memDat, c, d, s, trialTypes)
    for ii in range(len(trialTypes)-1): 
        out[list(trialTypes.keys())[ii] + '_d'] = params[5+ii]
        out[list(trialTypes.keys())[ii] + '_s'] = params[5+len(d)+ii]
        
    for ii, key in enumerate(trialTypes): 
        for conRating in list(range(6)): 
            out[key[0:-4] + str(conRating+1)] = confRat[ii,conRating]
            
    return(out)

def findT1Fam(name, bossNorms, bossNorms2):
    idx = np.where(bossNorms['FILENAME'] == name)
    
    if np.isnan(bossNorms['Fam_Mean'].iloc[idx]).any():
        idx = np.where(bossNorms2['FILENAME'].str.contains(name))[0]
        fam = bossNorms2['Fam_Mean'].iloc[idx]
    else:
        fam = bossNorms['Fam_Mean'].iloc[idx]
    
    if len(fam) == 1:
        fam = fam.iloc[0]
    else:
        fam = np.nan
    return(fam)
        
        
def readDataFile(path): 
    df = pd.read_csv(path)
    bossNorms = pd.read_csv(r"C:\Users\dtf8829\Documents\ABmem\BOSS_norms.csv")
    bossNorms2 = pd.read_csv(r"C:\Users\dtf8829\Documents\ABmem\BOSS_norms2.csv")
    #get the AB data
    ABDat = df[df['TrialType']=='AB']
    
    #need to make a column that indicates what the T1 image was
    ABDat['T1imageID'] = ABDat['T1_Identity']
    f = partial(findT1Fam, bossNorms=bossNorms, bossNorms2=bossNorms2)
    ABDat['T1fam'] = ABDat['T1_Identity'].apply(f)
    
    #set T1cor to 0 (wrong) for all 
    ABDat['T1cor'] = 0
    #set to 1 for trials on which T1 was present and the participant pressed 'left'
    ABDat.loc[(ABDat['key_resp_2.keys'] == 'left') & (ABDat['T1_PA'] == 'P'), 'T1cor'] = 1
    #set to 1 for trials on which T1 was absent and the participant pressed 'right'
    ABDat.loc[(ABDat['key_resp_2.keys'] == 'right') & (ABDat['T1_PA'] == 'A'), 'T1cor'] = 1
    
    
    #set T2cor to 0 (wrong) for all 
    ABDat.loc[:,'T2cor'] = 0
    #set to 1 for trials on which T1 was present and the participant pressed 'left'
    ABDat.loc[(ABDat['key_resp_3.keys'] == 'left') & (ABDat['T2_PA'] == 'P'), 'T2cor'] = 1
    #set to 1 for trials on which T1 was absent and the participant pressed 'right'
    ABDat.loc[(ABDat['key_resp_3.keys'] == 'right') & (ABDat['T2_PA'] == 'A'), 'T2cor'] = 1
    
   
    
    out = pd.Series()
    #Overall AB for T1
    out['ABT1overall'] = np.mean(ABDat['T1cor'])
    out['ABT2overall'] = np.mean(ABDat['T2cor'])
    
    #cherck for participants who got the instructions for responding backwards
    if out['ABT1overall'] < .5:
        T1Temp = np.zeros(len(ABDat['T1cor']))
        T1Temp[ABDat['T1cor'] == 0] = 1
        T1Temp[ABDat['T1cor'] == 1] = 0
        ABDat['T1cor'] = T1Temp
        
        T2Temp = np.zeros(len(ABDat['T1cor']))
        T2Temp[ABDat['T2cor'] == 0] = 1
        T2Temp[ABDat['T2cor'] == 1] = 0
        ABDat['T2cor'] = T2Temp
        
        #Overall AB for T1
        out['ABT1overall'] = np.mean(ABDat['T1cor'])
        out['ABT2overall'] = np.mean(ABDat['T2cor'])
    
    
    ABDat['PerformanceType'] = 'lag5'
    ABDat.loc[(ABDat['T2_PA'] == 'P') & (ABDat['T1_PA'] == 'P') & 
              (ABDat['T1cor'] == 1) & (ABDat['Lag']==1) & 
              (ABDat['T2cor']==1), 'PerformanceType'] = 'nonBlink'
    ABDat.loc[(ABDat['T2_PA'] == 'P') & (ABDat['T1_PA'] == 'P') & 
              (ABDat['T1cor'] == 1) & (ABDat['Lag']==1) & 
              (ABDat['T2cor']==0), 'PerformanceType'] = 'Blink'
    
    out['blinkFam'] = np.nanmean(ABDat.loc[ABDat['PerformanceType'] == 'Blink', 'T1fam'])
    out['nonBlinkFam'] = np.nanmean(ABDat.loc[ABDat['PerformanceType'] == 'nonBlink', 'T1fam'])
    #get the % correct for T1, T2, and T2|T1 for all trial types
    
    for cnd in ABDat['CondType'].unique():
        out[cnd +'_T1'] = np.mean(ABDat.loc[ABDat['CondType'] == cnd, 'T1cor'])
        out[cnd +'_T2'] = np.mean(ABDat.loc[ABDat['CondType'] == cnd, 'T2cor'])
        out[cnd +'_T2_T1'] = np.mean(ABDat.loc[ (ABDat['CondType'] == cnd) & (ABDat['T1cor'] == 1), 'T2cor'])
    

    
    #get the memory data
    memDat = df[df['TrialType'] == 'Mem']
 
    #create linkage between memory and attentional blink data
    f = lambda x: refMem2AB(x, ABDat)
    x = memDat.apply(f, axis = 1)
    memDat = pd.concat([memDat,x], axis = 1)
    memDat.dropna(axis = 1, thresh = 30, inplace = True)
    
   
    
    #key trial types: 
    memDat['ABrefType'] = 'skip'
    #previous lag1 P P both correct 'PP1cor'
    memDat.loc[(memDat['T1PA']=='P') & (memDat['ABT1cor']==1) & (memDat['ABT2cor']==1) & (memDat['Lag']==1), 'ABrefType'] = 'PP1cor'
    #previous lag1 P P both miss 'PP1miss'
    memDat.loc[(memDat['T1PA']=='P') & (memDat['ABT1cor']==0) & (memDat['ABT2cor']==0) & (memDat['Lag']==1), 'ABrefType'] = 'PP1miss'
    #previous lag1 P P blink 'PP1blink'
    memDat.loc[(memDat['T1PA']=='P') & (memDat['ABT1cor']==1) & (memDat['ABT2cor']==0) & (memDat['Lag']==1), 'ABrefType'] = 'PP1blink'
    #previous lag1 P P T2 correct 'PP1T2'
    memDat.loc[(memDat['T1PA']=='P') & (memDat['ABT1cor']==0) & (memDat['ABT2cor']==1) & (memDat['Lag']==1), 'ABrefType'] = 'PP1T2'
    #previous lag5 P P both correct 'PP5cor'
    memDat.loc[(memDat['T1PA']=='P') & (memDat['ABT1cor']==1) & (memDat['ABT2cor']==1) & (memDat['Lag']==5), 'ABrefType'] = 'PP5cor'
    #distracter
    memDat.loc[memDat['ImageType']=='Distractor', 'ABrefType'] = 'distracter'
    #novel
    memDat.loc[memDat['ImageType']=='Novel', 'ABrefType'] = 'novel'
    
    #assessing memory
    memDat['oldNew'] = 0
    memDat.loc[memDat['MemResponse.keys']>3, 'oldNew'] = 1
    memDat['memCor'] = 0
    memDat.loc[(memDat['ABrefType'] != 'novel') & (memDat['oldNew'] == 1), 'memCor'] = 1
    memDat.loc[(memDat['ABrefType'] == 'novel') & (memDat['oldNew'] == 0), 'memCor'] = 1
    
#    #test to see if the subject switched round the memory responses 
#    if np.mean(memDat.loc[memDat['ABrefType'] == 'PP5cor', 'memCor']) < .5:
#        tempCon = np.zeros(len(memDat['MemResponse.keys']))
#        conOrig = memDat['MemResponse.keys']
#        tempCon[conOrig==1] = 6
#        tempCon[conOrig==2] = 5
#        tempCon[conOrig==3] = 4
#        tempCon[conOrig==4] = 3
#        tempCon[conOrig==5] = 2
#        tempCon[conOrig==6] = 1
#        memDat['MemResponse.keys'] = tempCon
#        memDat['oldNew'] = 0
#        memDat.loc[memDat['MemResponse.keys']>3, 'oldNew'] = 1
#        memDat['memCor'] = 0
#        memDat.loc[(memDat['ABrefType'] != 'novel') & (memDat['oldNew'] == 1), 'memCor'] = 1
#        memDat.loc[(memDat['ABrefType'] == 'novel') & (memDat['oldNew'] == 0), 'memCor'] = 1
        
    
    out['targMem'] = np.mean(memDat.loc[(memDat['ABrefType'] != 'distracter') & (memDat['ABrefType'] != 'skip'), 'memCor'])
    out['distMem'] = np.mean(memDat.loc[(memDat['ABrefType'] == 'distracter') | (memDat['ABrefType'] == 'novel'), 'memCor'])
    out['T1Overall'] = np.mean(memDat.loc[(memDat['ABT1cor'] == 1), 'memCor'])
    out['T2Overall'] = np.mean(memDat.loc[(memDat['ABT2cor'] == 1), 'memCor'])
    
    for cnd in memDat['ABrefType'].unique():
        if (cnd != 'skip'): 
            out[cnd + '_mem'] = np.mean(memDat.loc[(memDat['ABrefType'] == cnd), 'memCor'])
    
    #fit an unequal variance signal detection model to the data
    #do this both for targets v. distractors v. novel items
    # and for: pp1cor v. pp1blink v. lag5 v. distractor v. novel items
    out = fitUVSD(memDat, out)
    
    out.name = df['date'][0]
    
    return(out)
    
    

    
#test = readDataFile(path)
    
    
    
    
    
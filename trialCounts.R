
#code for trial counting and confidence data table making: 

#assumes that a 6-point confidence scale is being kept rather than compressed

###########            trial counts and images used: ################## 
rawDat = read.csv(file_list[10])
ABdat = rawDat[grepl('_AB', rawDat$TrialType), ]
allPics = pivot_longer(ABdat, cols = 
                         names(test)[grepl('Img', names(test)) & 
                                       !grepl('Total', names(test))],
                       names_to = 'img', values_to = 'file') %>% 
  select(c('img', 'file'))

sum(grepl('tapler', unique(allPics$file)))

uniImages = unique(allPics$file)
uniCounts = sapply(uniImages, function(x) sum(grepl(x, allPics$file)))
uniImages = data.frame('imgID' = uniImages, 'count' = uniCounts)

#

memDat = rawDat[grepl('_MEM', rawDat$TrialType) & ! grepl('prac', rawDat$TrialType), ]
memDat = select(memDat, names(memDat)[!is.na(unlist(memDat[4,]))])
sum(memDat$CondType == 'new')
sum(memDat$priorLag == 1)

### get info for Table 3 #############################################

memoryTable = data.frame('targetType' = c('dual', 'single', 'dualLag1', 'dualLag5', 'singleLag1', 'singleLag5', 'T1|blink', 'T1|T2', 'novel'),
                         'accuracy' = rep('', 9), 
                         'd'          = rep('', 9), 
                         'SDratio'   = rep('', 9),
                         'con1'      = rep('', 9),
                         'con2'      = rep('', 9),
                         'con3'      = rep('', 9),
                         'con4'      = rep('', 9),
                         'con5'      = rep('', 9),
                         'con6'      = rep('', 9),
                         'trialCount'= rep('', 9))
typesPer = c('dual_', 'single_', 'dualLag1_', 'dualLag5_', 'singleLag1_', 'singleLag5_', 'blink_', 'noBlink_', 'novel_')
typesUVSD = c('dual_', 'single_', 'dualLag1_', 'dualLag5_', 'singleLag1_', 'singleLag5_', 'blink_', 'noBlink_', 'novel_')
varNames = c('acc', 'd', 's', 'count')
outNames = c('accuracy', 'd', 'SDratio', 'trialCount')
confShort = c('oldD', 'oldS', 'lag1D', 'lag5D', 'lag1S', 'lag5S', 'blink', 'noBlink', 'new')
for(ti in 1:9) { 
  for(vi in 1:4) {
    curName = paste0(typesPer[ti], varNames[vi])
    memoryTable[ti, outNames[vi]] = paste(as.character(round(mean(subSums[,curName]), 2)), ' (',
                                          as.character(round(sd(subSums[,curName]), 2)), ')', sep="")
  }
  for(ci in 1:6) {
    curName = paste0(confShort[ti], ci)
    curOut = paste0('con', ci)
    memoryTable[ti, curOut] = paste(as.character(round(mean(subSums[,curName]), 2)), ' (',
                                    as.character(round(sd(subSums[,curName]), 2)), ')', sep="")
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


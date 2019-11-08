###simulated worlds############################################################################
#wrapper functions

getSims <- function(mysimNu){
  
  #get all params
  getParams(parameterDF,mysimNu)
  
  if(samplingType=="Occurrence"){
    Occ <- getOccuHistory(occProb = OccProb, Seffects, Teffects, Ieffects,
                          Tcovariate,Scovariate)
    Obs <- getRepeatVisits(Occ=Occ, NVisits=NVisits, DetProb = DetProb, SiteDetEffects,YearDetEffects,IntDetEffects, Scovariate,Tcovariate,Noise)
    return(Obs)
    
  } else if(samplingType=="Abundance"){
    Abund <- getAbundHistory(lambda = lambda, Seffects, Teffects, Ieffects,
                             Tcovariate,Scovariate)
    Obs <- getRepeatAbundVisits(Abund=Abund, NVisits=NVisits, DetProb = DetProb, SiteDetEffects,YearDetEffects,IntDetEffects, Scovariate,Tcovariate,Noise)
    return(Obs)
  }
  
}


getModels <- function(Obs,type=1){
  
  modelData <- getSpartaFormat(Obs,focalSpecies,subsetPositive)
  bugs.data <- getBugsData(modelData)
  modelSummary <- runModel(bugs.data,modelData,focalSpecies,myModel)
  modelSummary$Type <- type

  return(modelSummary)
}


####simulation summaries###############################################################################

getObsSummaries <- function(df){
  require(plyr)
  
  #summaries per simulation
  perSim <- ldply(df,summarise,
                  
                  #true occurence
                  Occ_year1 = mean(Occurence[Year==min(Year)]),
                  Occ_last = mean(Occurence[Year==max(Year)]),
                
                  #true occurence of focal species
                  Occ_focal_year1 = mean(Occurence[Species==focalSpecies & Year==min(Year)]),
                  Occ_focal_last = mean(Occurence[Species==focalSpecies & Year==max(Year)]),
                  
                  #observation of all species
                  Obs_year1 = mean(Visit1[Year==min(Year)]),
                  Obs_last = mean(Visit1[Year==max(Year)]),
                  
                  #observation of focal species
                  Obs_focal_year1 = mean(Visit1[Species==focalSpecies & Year==min(Year)]),
                  Obs_focal_last = mean(Visit1[Species==focalSpecies & Year==max(Year)]),
  
                  #how many have all zeros
                  Occ_Zeros = all(Occurence==0),
                  Occ_Ones = all(Occurence==1))
  
  #summary overall
  colMeans(perSim)
  
}


combineObs <- function(x){
  
  temp <- ldply(x, function(x){
        
  focalDF <- subset(x,Species==focalSpecies)
  
  #trend estimate
  focalDFSummary <- ddply(focalDF,.(Year),summarise,nuSites=sum(Occurence))
  true.trend <- lm(nuSites/NSites ~ Year,data=focalDFSummary)$coef[2]

  true.first <- mean(x$Occurence[x$Year==min(x$Year) & x$Species==focalSpecies])
  true.last <- mean(x$Occurence[x$Year==max(x$Year) & x$Species==focalSpecies])
  true.change <- (true.last-true.first)/true.first 
  true.odds <- log((true.last/(1-true.last))/(true.first/(1-true.first)))
  
  return(data.frame(true.first,true.last,true.change,true.odds,true.trend))
  })
  
  return(temp)
  
  }


gettrueSummaries <- function(x){
  
  require(reshape2)
  temp <- melt(x)
  
  outTrueSummary <- ddply(temp,.(variable),summarise,
                         lowerQ = quantile(value,0.25),
                         medianQ = quantile(value,0.5),
                         upperQ = quantile(value,0.75))
  
  names(outTrueSummary)[1] <- "Param"
  
  return(outTrueSummary)
  
}

###rawData analysis#############################################################################################

#https://stats.idre.ucla.edu/r/faq/how-can-i-estimate-the-standard-error-of-transformed-regression-parameters-in-r-using-the-delta-method/
rawdataAnalysis <- function(Obs){
  
  modelData <- getSpartaFormat(Obs,focalSpecies=focalSpecies,subsetPositive = TRUE)
  
  #glmer - mixed effect model
  require(lme4)
  #transform into probabilities
  require(boot)
  e <- 0.00000001

    
    lmer1 <- glmer(Obs ~ (factor(Year)-1)+(1|Site),data=subset(modelData,!is.na(Obs)),
                 family=binomial)
  
    #extract output
    
    #first and last
    lmerPsi_First <- summary(lmer1)$coefficients[1,1]
    lmerPsi_First_se <- summary(lmer1)$coefficients[1,2]
    len <- nrow(summary(lmer1)$coefficients)
    lmerPsi_Last <- summary(lmer1)$coefficients[len,1]
    lmerPsi_Last_se <- summary(lmer1)$coefficients[(nrow(summary(lmer1)$coefficients)),2]
    
    #change parameters
    lmerPsi_Odds <- lmerPsi_Last - lmerPsi_First
    lmerPsi_Odds_se <- sqrt(lmerPsi_Last_se^2 + lmerPsi_First_se^2)
    
    #estimate se using delta method
    require(msm)
    lmerPsi_First_se <- deltamethod(~ exp(x1)/(1+exp(x1)), lmerPsi_First, vcov(lmer1)[1,1])
    lmerPsi_Last_se <- deltamethod(~ exp(x1)/(1+exp(x1)), lmerPsi_Last, vcov(lmer1)[len,len])
    
    #change first and last into probability
    lmerPsi_First <- inv.logit(lmerPsi_First)
    lmerPsi_Last <- inv.logit(lmerPsi_Last)
    lmerPsi_Change <- (lmerPsi_Last-lmerPsi_First)/lmerPsi_First
    
    temp <- data.frame(raw.first=lmerPsi_First,
                       raw.last=lmerPsi_Last,
                       raw.change=lmerPsi_Change,
                       raw.odds=lmerPsi_Odds)
                       #raw.first.se=lmerPsi_First_se,
                       #raw.last.se=lmerPsi_Last_se,
                       #raw.odds.se=lmerPsi_Odds_se)
    
    #simple trend model
    
    #trend in proportion of sites occupied??
    
    lmer1 <- glmer(Obs ~ Year+(1|Site),data=subset(modelData,!is.na(Obs)),
                   family=binomial)

    lmerPsi_trend <- summary(lmer1)$coefficients[2,1]
    lmerPsi_trend_se <- summary(lmer1)$coefficients[2,2]
    lmerPsi_trend_z <- summary(lmer1)$coefficients[2,3]
    lmerPsi_trend_sig <- summary(lmer1)$coefficients[2,4]
    
    temp2 <- data.frame(raw.trend=lmerPsi_trend,
                       raw.Ztrend=lmerPsi_trend_z,
                       raw.Sigtrend=lmerPsi_trend_sig)
    

  return(cbind(temp,temp2))
}

#also add a raw poisson analysis################

getrawSummaries <- function(outRaw){
  
  require(reshape2)
  outRawMelted <- melt(outRaw)
  outRawSummary <- ddply(outRawMelted,.(variable),summarise,
                         lowerQ = quantile(value,0.25),
                         medianQ = quantile(value,0.5),
                         upperQ = quantile(value,0.75))
  
  names(outRawSummary)[1] <- "Param"
  
  return(outRawSummary)
  
}

###combining models###############################################################################

combineOut <- function(out){
  
  outCombined <- do.call(rbind,out)
  outCombined$Param[which(outCombined$Param=="first.psi")] <- "model.first"
  outCombined$Param[which(outCombined$Param=="last.psi")] <- "model.last"
  outCombined$Param[which(outCombined$Param=="mean.psi.change")] <- "model.change"
  outCombined$Param[which(outCombined$Param=="odds.psi.change")] <- "model.odds"
  outCombined$Param[which(outCombined$Param=="regres.psi")] <- "model.trend"
  
  #add simulation number
  j = 0
  outCombined$simRep <- NA
  for(i in 1:length(outCombined$Param)){
    j = ifelse(outCombined$Param[i]=="mean.psi",j+1,j)
    outCombined$simRep[i] <- j
  }
  
  return(outCombined)
  
}

###compare results##########################################################################

checkResults <- function(out,Obs){
  
  allModels <- combineOut(out)
  allTruth <- combineObs(Obs)
  
  temp <- subset(allModels,Param=="model.last")
  temp <- cbind(temp,allTruth)
  g1 <- ggplot(temp)+
          geom_point(aes(x=mean,y=true.last))+
          geom_abline(slope=1,intercept=0)+
          xlim(0,1)+ylim(0,1)
  
  temp <- subset(allModels,Param=="model.change")
  temp <- cbind(temp,allTruth)
  g2<- ggplot(temp)+
          geom_point(aes(x=mean,y=true.change))+
          geom_abline(slope=1,intercept=0)
  
  temp <- subset(allModels,Param=="model.odds")
  temp <- cbind(temp,allTruth)
  g3 <- ggplot(temp)+
          geom_point(aes(x=mean,y=true.odds))+
          geom_abline(slope=1,intercept=0)
  
  print(g1)
  print(g2)
  print(g3)
  
}


getModelSummaries <- function(out){
  
  outCombined <- combineOut(out)
  
  outSummaries <- subset(outCombined,Param %in%c ("model.change","model.odds",
                                                  "model.last","model.first",
                                                  "model.trend"))
  outSummaries <- ddply(outSummaries,.(Param),summarise,
                        lowerQ=quantile(mean,0.25),
                        medianQ=quantile(mean,0.5),
                        upperQ=quantile(mean,0.75))
  
  return(outSummaries)
}



plotComparison <- function(modelSummaries,rawSummaries,trueSummaries){
  
  temp <- rbind(modelSummaries,rawSummaries,trueSummaries)

  temp$model <- sapply(temp$Param,function(x)strsplit(x,"\\.")[[1]][1])
  temp$type <- sapply(temp$Param,function(x)strsplit(x,"\\.")[[1]][2])
  
  print(ggplot(temp)+
          geom_crossbar(aes(x=model,y=medianQ,ymin=lowerQ,ymax=upperQ),width=0.4)+
          facet_wrap(~type,scales="free")+
          ylab("change in occupancy")+
          theme_bw())
  
  return(temp)
}


# plotSE <- function(out,outRaw){
#   
#   #for occupancy detection model
#   outCombined <- combineOut(out)
#   outSummaries <- subset(outCombined,Param %in%c ("model.odds","model.last"))
#   outSummaries <- ddply(outSummaries,.(Param),summarise,
#                         lowerQ=quantile(sd,0.25),
#                         medianQ=quantile(sd,0.5),
#                         upperQ=quantile(sd,0.75))
#   
#   #from raw data
#   Param <- c("basic.odds","basic.last")
#   lowerQ <- c(quantile(outRaw$lmerPsi_Odds_se,0.25),
#               quantile(outRaw$lmerPsi_Last_se,0.25))
#   
#   medianQ <- c(quantile(outRaw$lmerPsi_Odds_se,0.5),
#                quantile(outRaw$lmerPsi_Last_se,0.5))
#   
#   upperQ <- c(quantile(outRaw$lmerPsi_Odds_se,0.75),
#               quantile(outRaw$lmerPsi_Last_se,0.75))
#   
#   temp <- data.frame(Param,lowerQ,medianQ,upperQ)
#   
#   allSE <- rbind(outSummaries,temp)
#   allSE$Model <- sapply(allSE$Param,function(x)strsplit(x,"\\.")[[1]][1])
#   allSE$Param <- sapply(allSE$Param,function(x)strsplit(x,"\\.")[[1]][2])
#   
#   print(ggplot(allSE)+
#           geom_crossbar(aes(x=Model,y=medianQ,ymin=lowerQ,ymax=upperQ),width=0.4)+
#           facet_wrap(~Param,scales="free")+
#           ylab("Parameter")+
#           theme_bw())
#   
# }

###difference plots##################################################################

plotDifference <- function(Obs,out,outRaw){
  
  require(reshape2)
  outObs <- combineObs(Obs)
  outCombined <- combineOut(out)
  outModels <- acast(outCombined,simRep~Param,value.var="mean")
  allData <- cbind(outModels,outObs,outRaw)
  
  #first occupancy
  firstData <- allData[,grepl("first",names(allData))]
  firstData$model_vs_true <- firstData$model.first - firstData$true.first
  firstData$raw_vs_true <- firstData$raw.first - firstData$true.first
  
  #last occupancy
  lastData <- allData[,grepl("last",names(allData))]
  lastData$model_vs_true <- lastData$model.last - lastData$true.last
  lastData$raw_vs_true <- lastData$raw.last - lastData$true.last
  
  #change occupancy
  changeData <- allData[,grepl("change",names(allData))]
  changeData$model_vs_true <- changeData$model.change - changeData$true.change
  changeData$raw_vs_true <- changeData$raw.change - changeData$true.change
  
  #odd occupancy
  oddsData <- allData[,grepl("odds",names(allData))]
  oddsData$model_vs_true <- oddsData$model.odds - oddsData$true.odds
  oddsData$raw_vs_true <- oddsData$raw.odds - oddsData$true.odds

  #melt data frame
  require(reshape2)
  firstDataM <- melt(firstData[,c("model_vs_true","raw_vs_true")])
  firstDataM$Param <- "first"
  lastDataM <- melt(lastData[,c("model_vs_true","raw_vs_true")])
  lastDataM$Param <- "last"
  changeDataM <- melt(changeData[,c("model_vs_true","raw_vs_true")])
  changeDataM$Param <- "change"
  oddsDataM <- melt(oddsData[,c("model_vs_true","raw_vs_true")])
  oddsDataM$Param <- "odds"
  
  
  #combine all
  df <- rbind(firstDataM,lastDataM,changeDataM,oddsDataM)
  names(df)[1:2] <- c("Comparison","Difference")
  df$Param <- factor(df$Param, levels=c("first","last","change","odds"))
  require(ggplot2)
  g1 <- ggplot(df)+
    geom_boxplot(aes(x=Comparison,y=Difference),outlier.shape = NA)+
    facet_wrap(~Param,scales="free",ncol=1)+
    geom_hline(yintercept=0,color="red",linetype="dashed")+
    coord_flip()+
    theme_bw()
  print(g1)
}

###results saving###########################################################################

saveResults <- function(out,type="fixed_interceptonly"){
  
  allModels <- combineOut(out)
  print(summary(allModels$Rhat))
  allModels <- subset(allModels,Rhat<1.1)
  
  g1 <- ggplot(subset(allModels,Param=="model.last"))+
          geom_point(aes(x=mean,y=truePsiLast))+
          geom_abline(slope=1,intercept=0)
          #xlim(0,1)+ylim(0,1)
  
  g2 <- ggplot(subset(allModels,Param=="model.change"))+
          geom_point(aes(x=mean,y=truePsiChange))+
          geom_abline(slope=1,intercept=0)
          #xlim(0,1)+ylim(0,1)

  g3 <- ggplot(subset(allModels,Param=="models.odds"))+
          geom_point(aes(y=mean,x=truePsiOdds))+
          geom_hline(yintercept=0)

  require(cowplot)
  plot_grid(g1,g2,g3)
  ggsave(paste0("plots/",type,"_",systemTime,".png"))
  
}
#random###########################################################################

#https://www.andrewheiss.com/blog/2016/04/25/convert-logistic-regression-standard-errors-to-odds-ratios-with-r/

get.or.se <- function(model) {
  or = exp(summary(model)$coefficients[,1])
  var.diag = diag(vcov(model))
  or.se = sqrt(or^2 * var.diag) 
  return(or.se)
}


####power analysis########################################################

#number of times a significant trend was detected

powerAnalysis <- function(out,outRaw){

#raw analysis
  powerRaw <- nrow(subset(outRawT,raw.Sigtrend<0.05))
 
#jags model
  outCombined <- combineOut(out)
  outCombined$Sig <- 0
  outCombined$Sig[outCombined$X2.5.<0 & outCombined$X97.5.<0] <- 1
  outCombined$Sig[outCombined$X2.5.>0 & outCombined$X97.5.>0] <- 1
  powerModel <- sum(outCombined$Sig[outCombined$Param=="model.trend"])
  N=length(outCombined$Sig[outCombined$Param=="model.trend"])

  #combine
  temp <- data.frame(model=c("raw","model"),power=c(powerRaw/N,powerModel/N))
  
  #plot
  require(ggplot2)
  print(ggplot(temp,aes(x=model,y=power))+
    geom_bar(stat="identity"))
  
  return(temp)
  
}



###########################################################################









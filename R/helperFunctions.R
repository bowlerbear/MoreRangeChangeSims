library(plyr)
###simulated worlds############################################################################
#wrapper functions

getSims <- function(){
  
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

getModels <- function(Obs,bias="standard",myModel){
  
  modelData <- getSpartaFormat(Obs,focalSpecies,subsetPositive)
  bugs.data <- getBugsData(modelData)
  modelSummary <- runModel(bugs.data,modelData,focalSpecies,myModel)
  modelSummary$Bias <- bias

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

###raw data analysis#############################################################################################

#https://stats.idre.ucla.edu/r/faq/how-can-i-estimate-the-standard-error-of-transformed-regression-parameters-in-r-using-the-delta-method/
rawdataAnalysis <- function(Obs){
  
  modelData <- getSpartaFormat(Obs,focalSpecies=focalSpecies,subsetPositive = TRUE)
  
  #glmer - mixed effect model
  require(lme4)
  #transform into probabilities
  require(boot)
  e <- 0.00000001

  #simplify data into presence/absence for each year
    siteData <- ddply(modelData,.(Year,Site),summarise,
                      P=length(unique(Occurence[!is.na(Obs) & Obs>0])))
  
    #fixed effect year model    
    lmer1 <- glm(P ~ (factor(Year)-1),data=siteData,family=binomial)
  
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
    lmer1 <- glm(P ~ Year,data=siteData,family=binomial)

    #predict occurence in final year
    last <- predict(lmer1,newdata=data.frame(Year=NYears),type="response")
    #predict occurence in first year
    first <- predict(lmer1,newdata=data.frame(Year=1),type="response")
    lmerPsi_trend = (last-first)/NYears
    
    #power stats
    lmerPsi_trend_z <- summary(lmer1)$coefficients[2,3]
    lmerPsi_trend_sig <- summary(lmer1)$coefficients[2,4]
    
    
    
    temp2 <- data.frame(raw.trend=lmerPsi_trend,
                       raw.Ztrend=lmerPsi_trend_z,
                       raw.Sigtrend=lmerPsi_trend_sig)
    
    
    #total sites poisson trend model
    siteData <- ddply(modelData,.(Year),summarise,
                      nSites=length(unique(Site[!is.na(Obs) & Obs>0])))
    siteData$totalSites <- length(unique(modelData$Site))

    glm1 <- glm(cbind(nSites,totalSites-nSites) ~ Year,data=siteData,family=binomial)
    #predict occurence in final year
    last <- predict(glm1,newdata=data.frame(Year=NYears),type="response")
    #predict occurence in first year
    first <- predict(glm1,newdata=data.frame(Year=1),type="response")
    glmPsi_trend = (last-first)/NYears
    
    glmPsi_trend_z <- summary(glm1)$coefficients[2,3]
    glmPsi_trend_sig <- summary(glm1)$coefficients[2,4]
    
    temp3 <- data.frame(rawP.trend=glmPsi_trend,
                        rawP.Ztrend=glmPsi_trend_z,
                        rawP.Sigtrend=glmPsi_trend_sig)

  return(cbind(temp,temp2,temp3))
}

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

getModelSummaries <- function(out){
  
  outCombined <- combineOut(out)
  
  outSummaries <- subset(outCombined,Param %in%c ("model.change","model.odds",
                                                  "model.last","model.first",
                                                  "model.trend"))
  
  outSummaries <- ddply(outSummaries,.(Param,Bias),summarise,
                        lowerQ=quantile(mean,0.25),
                        medianQ=quantile(mean,0.5),
                        upperQ=quantile(mean,0.75),
                        lowerQ_sd=quantile(sd,0.25),
                        medianQ_sd=quantile(sd,0.5),
                        upperQ_sd=quantile(sd,0.75))
  
  return(outSummaries)
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

  #trend in occupancy
  trendData <- allData[,grepl("trend",names(allData))]
  trendData$model_vs_true <- trendData$model.trend - trendData$true.trend
  trendData$raw_vs_true <- trendData$raw.trend - trendData$true.trend
  
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
  trendDataM <- melt(trendData[,c("model_vs_true","raw_vs_true")])
  trendDataM$Param <- "trend"
  
  #combine all
  df <- rbind(firstDataM,lastDataM,changeDataM,oddsDataM,trendDataM)
  names(df)[1:2] <- c("Comparison","Difference")
  df$Param <- factor(df$Param, levels=c("first","last","change","odds","trend"))
  
  require(ggplot2)
  g1 <- ggplot(df)+
    geom_boxplot(aes(x=Comparison,y=Difference),outlier.shape = NA)+
    facet_wrap(~Param,scales="free",ncol=1)+
    geom_hline(yintercept=0,color="red",linetype="dashed")+
    coord_flip()+
    theme_bw()
  print(g1)
}



getDifferences <- function(Obs,out,outRaw){
  
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
  
  #trend in occupancy
  trendData <- allData[,grepl("trend",names(allData))]
  trendData$model_vs_true <- trendData$model.trend - trendData$true.trend
  trendData$raw_vs_true <- trendData$raw.trend - trendData$true.trend
  
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
  trendDataM <- melt(trendData[,c("model_vs_true","raw_vs_true")])
  trendDataM$Param <- "trend"
  
  #combine all
  df <- rbind(firstDataM,lastDataM,changeDataM,oddsDataM,trendDataM)
  names(df)[1:2] <- c("Comparison","Difference")
  df$Param <- factor(df$Param, levels=c("first","last","change","odds","trend"))
  
  return(df)
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



###end########################################################################









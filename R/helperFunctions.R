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
  
  #get actual change
  e <- 0.00000001
  modelSummary$truePsiFirst <- mean(Obs$Occurence[Obs$Year==min(Obs$Year) & Obs$Species==focalSpecies])
  modelSummary$truePsiLast <- mean(Obs$Occurence[Obs$Year==max(Obs$Year) & Obs$Species==focalSpecies])
  modelSummary$truePsiChange <- modelSummary$truePsiLast/modelSummary$truePsiFirst 
  modelSummary$truePsiOdds <- log((modelSummary$truePsiLast)/(1-modelSummary$truePsiLast))/((modelSummary$truePsiFirst)/(1-modelSummary$truePsiFirst))
  
  return(modelSummary)
}

####simulation summaries###############################################################################

getObsSummaries <- function(df){
  require(plyr)
  
  
  #summaries per simulation
  perSim <- ldply(df,summarise,
                  #year1
                  Occ_year1 = mean(Occurence[Year==min(Year)]),
                  Occ_focal_year1 = mean(Occurence[Species==focalSpecies & Year==min(Year)]),
                  Obs_year1 = mean(Visit1[Year==min(Year)]),
                  Obs_focal_year1 = mean(Visit1[Species==focalSpecies & Year==min(Year)]),
                  
                  #last year
                  Occ_last = mean(Occurence[Year==max(Year)]),
                  Occ_focal_last = mean(Occurence[Species==focalSpecies & Year==max(Year)]),
                  Obs_last = mean(Visit1[Year==max(Year)]),
                  Obs_focal_last = mean(Visit1[Species==focalSpecies & Year==max(Year)]))
  
  #summary overall
  colMeans(perSim)
  
}

###rawData analysis#############################################################################################

#https://stats.idre.ucla.edu/r/faq/how-can-i-estimate-the-standard-error-of-transformed-regression-parameters-in-r-using-the-delta-method/
rawdataAnalysis <- function(Obs){
  
  modelData <- getSpartaFormat(Obs,focalSpecies=focalSpecies,subsetPositive = TRUE)
  
  #glmer - mixed effect model
  e <- 0.00000001
  require(lme4)
  #transform into probabilities
  require(boot)
  
  lmer1 <- glmer(Obs ~ (factor(Year)-1)+(1|Site),data=subset(modelData,!is.na(Obs)),
                 family=binomial)
  
  #first and last
  lmerPsiFirst <- summary(lmer1)$coefficients[1,1]
  lmerPsiFirst_se <- summary(lmer1)$coefficients[1,2]
  len <- nrow(summary(lmer1)$coefficients)
  lmerPsiLast <- summary(lmer1)$coefficients[len,1]
  lmerPsiLast_se <- summary(lmer1)$coefficients[(nrow(summary(lmer1)$coefficients)),2]
  
  #change parameters
  lmerPsiOdds <- lmerPsiLast - lmerPsiFirst
  lmerPsiOdds_se <- sqrt(lmerPsiLast_se^2 + lmerPsiFirst_se^2)
  
  #estimate se using delta method
  require(msm)
  lmerPsiFirst_se <- deltamethod(~ exp(x1)/(1+exp(x1)), lmerPsiFirst, vcov(lmer1)[1,1])
  lmerPsiLast_se <- deltamethod(~ exp(x1)/(1+exp(x1)), lmerPsiLast, vcov(lmer1)[len,len])
  
  #change first and last into probability
  lmerPsiFirst <- inv.logit(lmerPsiFirst)
  lmerPsiLast <- inv.logit(lmerPsiLast)
  lmerPsiChange <- lmerPsiLast/lmerPsiFirst
  
  temp <- data.frame(lmerPsiFirst,lmerPsiLast,lmerPsiChange,lmerPsiOdds,
                     lmerPsiFirst_se,lmerPsiLast_se,lmerPsiOdds_se)
  
  return(temp)
}

getrawSummaries <- function(outRaw){
  Param <- c("basic.mean","basic.odds","basic.last","basic.first")
  lowerQ <- c(quantile(outRaw$lmerPsiChange,0.25),
              quantile(outRaw$lmerPsiOdds,0.25),
              quantile(outRaw$lmerPsiLast,0.25),
              quantile(outRaw$lmerPsiFirst,0.25))
  
  medianQ <- c(quantile(outRaw$lmerPsiChange,0.5),
               quantile(outRaw$lmerPsiOdds,0.5),
               quantile(outRaw$lmerPsiLast,0.5),
               quantile(outRaw$lmerPsiFirst,0.5))
  
  upperQ <- c(quantile(outRaw$lmerPsiChange,0.75),
              quantile(outRaw$lmerPsiOdds,0.75),
              quantile(outRaw$lmerPsiLast,0.75),
              quantile(outRaw$lmerPsiFirst,0.75))
  
  data.frame(Param,lowerQ,medianQ,upperQ)
  
}

###comparing results###############################################################################

combineOut <- function(out){
  
  outCombined <- do.call(rbind,out)
  outCombined$Param[which(outCombined$Param=="first.psi")] <- "model.first"
  outCombined$Param[which(outCombined$Param=="last.psi")] <- "model.last"
  outCombined$Param[which(outCombined$Param=="mean.psi.change")] <- "model.change"
  outCombined$Param[which(outCombined$Param=="odds.psi.change")] <- "model.odds"
  return(outCombined)
  
}

checkResults <- function(out){
  
  allModels <- combineOut(out)
  #allModels <- subset(allModels,Rhat<1.1)
  
  print(ggplot(subset(allModels,Param=="model.last"))+
          geom_point(aes(x=mean,y=truePsiLast))+
          geom_abline(slope=1,intercept=0)+
          xlim(0,1)+ylim(0,1))
  
  print(ggplot(subset(allModels,Param=="model.change"))+
          geom_point(aes(x=mean,y=truePsiChange))+
          geom_abline(slope=1,intercept=0))
  
  print(ggplot(subset(allModels,Param=="model.odds"))+
          geom_point(aes(x=mean,y=truePsiOdds))+
          geom_abline(slope=1,intercept=0))
  
}


getModelSummaries <- function(out){
  
  outCombined <- combineOut(out)
  
  outSummaries <- subset(outCombined,Param %in%c ("model.change","model.odds",
                                                  "model.last","model.first"))
  outSummaries <- ddply(outSummaries,.(Param),summarise,
                        lowerQ=quantile(mean,0.25),
                        medianQ=quantile(mean,0.5),
                        upperQ=quantile(mean,0.75))
  
  #get raw data
  Param <- c("raw.mean","raw.odds","raw.last","raw.first")
  lowerQ <- c(quantile(outCombined$truePsiChange[outCombined$Param=="model.change"],0.25),
              quantile(outCombined$truePsiOdds[outCombined$Param=="model.odds"],0.25),
              quantile(outCombined$truePsiLast[outCombined$Param=="model.last"],0.25),
              quantile(outCombined$truePsiFirst[outCombined$Param=="model.first"],0.25))
  
  medianQ <- c(quantile(outCombined$truePsiChange[outCombined$Param=="model.change"],0.5),
               quantile(outCombined$truePsiOdds[outCombined$Param=="model.odds"],0.5),
               quantile(outCombined$truePsiLast[outCombined$Param=="model.last"],0.5),
               quantile(outCombined$truePsiFirst[outCombined$Param=="model.first"],0.5))
  
  upperQ <- c(quantile(outCombined$truePsiChange[outCombined$Param=="model.change"],0.75),
              quantile(outCombined$truePsiOdds[outCombined$Param=="model.odds"],0.75),
              quantile(outCombined$truePsiLast[outCombined$Param=="model.last"],0.75),
              quantile(outCombined$truePsiFirst[outCombined$Param=="model.first"],0.75))
  
  rawDF <- data.frame(Param,lowerQ,medianQ,upperQ)
  
  return(rbind(outSummaries,rawDF))
}



plotComparison <- function(modelSummaries,rawSummaries){
  
  temp <- rbind(modelSummaries,rawSummaries)
  
  temp$type <- as.character(sapply(temp$Param,function(x){
    
    if(grepl("odds",x)){ 
      "odds"
    }else if (grepl("last",x)){
      "last occ"
    }else if (grepl("first",x)){
      "first occ"
    } else {
      "change"
    }}))
  
  temp$model <- c("OD","OD","OD","OD","truth","truth","truth","truth",
                  "naive","naive","naive","naive")
  
  temp$model <- factor(temp$model, levels=c("truth","naive","OD"))
  temp$type <- factor(temp$type, levels=c("first occ","last occ","change","odds"))
  
  print(ggplot(temp)+
          geom_crossbar(aes(x=model,y=medianQ,ymin=lowerQ,ymax=upperQ),width=0.4)+
          facet_wrap(~type,scales="free")+
          ylab("change in occupancy")+
          theme_bw())
  
  return(temp)
}


plotSE <- function(out,outRaw){
  
  #for occupancy detection model
  outCombined <- combineOut(out)
  outSummaries <- subset(outCombined,Param %in%c ("model.odds","model.last"))
  outSummaries <- ddply(outSummaries,.(Param),summarise,
                        lowerQ=quantile(sd,0.25),
                        medianQ=quantile(sd,0.5),
                        upperQ=quantile(sd,0.75))
  
  #from raw data
  Param <- c("basic.odds","basic.last")
  lowerQ <- c(quantile(outRaw$lmerPsiOdds_se,0.25),
              quantile(outRaw$lmerPsiLast_se,0.25))
  
  medianQ <- c(quantile(outRaw$lmerPsiOdds_se,0.5),
               quantile(outRaw$lmerPsiLast_se,0.5))
  
  upperQ <- c(quantile(outRaw$lmerPsiOdds_se,0.75),
              quantile(outRaw$lmerPsiLast_se,0.75))
  
  temp <- data.frame(Param,lowerQ,medianQ,upperQ)
  
  allSE <- rbind(outSummaries,temp)
  allSE$Model <- sapply(allSE$Param,function(x)strsplit(x,"\\.")[[1]][1])
  allSE$Param <- sapply(allSE$Param,function(x)strsplit(x,"\\.")[[1]][2])
  
  print(ggplot(allSE)+
          geom_crossbar(aes(x=Model,y=medianQ,ymin=lowerQ,ymax=upperQ),width=0.4)+
          facet_wrap(~Param,scales="free")+
          ylab("Parameter")+
          theme_bw())
  
}

###difference plots##################################################################

plotDifference <- function(out,outRaw){
  outCombined <- combineOut(out)
  
  #first occupancy
  firstData <- subset(outCombined,Param=="model.first")
  firstData$basic.first <- outRaw$lmerPsiFirst
  names(firstData)[which(names(firstData)=="mean")] <- "model.first"
  firstData$model_vs_true <- firstData$model.first - firstData$truePsiFirst
  firstData$basic_vs_true <- firstData$basic.first - firstData$truePsiFirst
  
  #last occupancy
  lastData <- subset(outCombined,Param=="model.last")
  lastData$basic.last <- outRaw$lmerPsiLast
  names(lastData)[which(names(lastData)=="mean")] <- "model.last"
  lastData$model_vs_true <- lastData$model.last - lastData$truePsiLast
  lastData$basic_vs_true <- lastData$basic.last - lastData$truePsiLast
  
  #change occupancy
  changeData <- subset(outCombined,Param=="model.change")
  changeData$basic.change <- outRaw$lmerPsiChange
  names(changeData)[which(names(changeData)=="mean")] <- "model.change"
  changeData$model_vs_true <- changeData$model.change - changeData$truePsiChange
  changeData$basic_vs_true <- changeData$basic.change - changeData$truePsiChange
  
  #odd occupancy
  oddsData <- subset(outCombined,Param=="model.odds")
  oddsData$basic.odds <- outRaw$lmerPsiOdds
  names(oddsData)[which(names(oddsData)=="mean")] <- "model.odds"
  oddsData$model_vs_true <- oddsData$model.odds - oddsData$truePsiOdds
  oddsData$basic_vs_true <- oddsData$basic.odds - oddsData$truePsiOdds

  #melt data frame
  require(reshape2)
  firstDataM <- melt(firstData[,c("model_vs_true","basic_vs_true")])
  firstDataM$Param <- "first"
  lastDataM <- melt(lastData[,c("model_vs_true","basic_vs_true")])
  lastDataM$Param <- "last"
  changeDataM <- melt(changeData[,c("model_vs_true","basic_vs_true")])
  changeDataM$Param <- "change"
  oddsDataM <- melt(oddsData[,c("model_vs_true","basic_vs_true")])
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












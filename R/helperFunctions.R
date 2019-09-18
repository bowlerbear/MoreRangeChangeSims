checkResults <- function(out){
  
  allModels <- do.call(rbind,out)
  print(summary(allModels$Rhat))
  allModels <- subset(allModels,Rhat<1.05)
  
  print(ggplot(subset(allModels,Param=="first.psi"))+
    geom_point(aes(x=mean,y=truePsiFirst))+
    geom_abline(slope=1,intercept=0)+
    xlim(0,1)+ylim(0,1))
  
  print(ggplot(subset(allModels,Param=="mean.psi.change"))+
    geom_point(aes(x=mean,y=truePsiChange))+
    geom_abline(slope=1,intercept=0)+
    xlim(0,1)+ylim(0,1))
  
  print(ggplot(subset(allModels,Param=="odds.psi.change"))+
    geom_point(aes(x=mean,y=truePsiOdds))+
    geom_abline(slope=1,intercept=0)+
    xlim(0,1)+ylim(0,1))

  
  #difference plots
  allModels$firstDiff <- allModels$truePsiFirst - allModels$mean
  print(ggplot(subset(allModels,Param=="first.psi"))+
          geom_point(aes(y=firstDiff,x=truePsiFirst))+
          geom_hline(yintercept=0))
  
  allModels$Change <- allModels$truePsiChange - allModels$mean
  print(ggplot(subset(allModels,Param=="mean.psi.change"))+
          geom_point(aes(y=Change,x=truePsiChange))+
          geom_hline(yintercept=0))
  
}

saveResults <- function(out,type="fixed_interceptonly"){
  
  allModels <- do.call(rbind,out)
  print(summary(allModels$Rhat))
  allModels <- subset(allModels,Rhat<1.05)
  
  g1 <- ggplot(subset(allModels,Param=="first.psi"))+
          geom_point(aes(x=mean,y=truePsiFirst))+
          geom_abline(slope=1,intercept=0)
          #xlim(0,1)+ylim(0,1)
  
  g2 <- ggplot(subset(allModels,Param=="mean.psi.change"))+
          geom_point(aes(x=mean,y=truePsiChange))+
          geom_abline(slope=1,intercept=0)
          #xlim(0,1)+ylim(0,1)
  
  allModels$firstDiff <- allModels$truePsiFirst - allModels$mean
  g3 <- ggplot(subset(allModels,Param=="first.psi"))+
          geom_point(aes(y=firstDiff,x=truePsiFirst))+
          geom_hline(yintercept=0)
  
  allModels$Change <- allModels$truePsiChange - allModels$mean
  g4 <- ggplot(subset(allModels,Param=="mean.psi.change"))+
          geom_point(aes(y=Change,x=truePsiChange))+
          geom_hline(yintercept=0)
  
  
  require(cowplot)
  plot_grid(g1,g2,g3,g4,nrow=2)
  ggsave(paste0("plots/",type,".png"))
  
}


getObsSummaries <- function(df){
require(plyr)
  
    
  #summaries per simulation
  perSim <- ldply(df,summarise,
                  #year1
                  nuOcc_year1 = mean(Occurence[Year==min(Year)]),
                  meanOcc_focal_year1 = mean(Occurence[Species==focalSpecies & Year==min(Year)]),
                  nuObs_year1 = mean(Visit1[Year==min(Year)]),
                  nuObs_focal_year1 = mean(Visit1[Species==focalSpecies & Year==min(Year)]),
                  #last year
                  nuOcc_last = mean(Occurence[Year==max(Year)]),
                  meanOcc_focal_last = mean(Occurence[Species==focalSpecies & Year==max(Year)]),
                  nuObs_last = mean(Visit1[Year==max(Year)]),
                  nuObs_focal_last = mean(Visit1[Species==focalSpecies & Year==max(Year)]))
  
                  #summary overall
                  colMeans(perSim)
  
}

#getSims <- function(x){
#  Occ <- getOccuHistory(occProb = OccProb, spatialbeta=Seffects, interactionbeta=Ieffects, spatialcovariates=Scovariate,temporalbeta=Teffects, temporalcovariates=Tcovariate)
#  Obs <- getRepeatVisits(Occ=Occ, NVisits=NVisits, DetProb = DetProb, SiteDetEffects,YearDetEffects,IntDetEffects, Scovariate,Tcovariate)
#  return(Obs)
#}

getSims <- function(x){
  Abund <- getAbundHistory(lambda = lambda, spatialbeta=Seffects, interactionbeta=Ieffects, spatialcovariates=Scovariate,temporalbeta=Teffects, temporalcovariates=Tcovariate)
  Obs <- getRepeatAbundVisits(Abund=Abund, NVisits=NVisits, DetProb = DetProb, SiteDetEffects,YearDetEffects,IntDetEffects, Scovariate,Tcovariate)
  return(Obs)
  }

getModels <- function(Obs,type="Standard"){
  modelData <- getSpartaFormat(Obs,focalSpecies=focalSpecies,subsetPositive = TRUE)
  bugs.data <- getBugsData(modelData)
  modelSummary <- runModel(bugs.data,modelData,focalSpecies,myModel)
  modelSummary$Type <- type
  
  #get actual change
  e <- 0.00000001
  modelSummary$truePsiFirst <- mean(Obs$Occurence[Obs$Year==min(Obs$Year) & Obs$Species==focalSpecies])
  modelSummary$truePsiLast <- mean(Obs$Occurence[Obs$Year==max(Obs$Year) & Obs$Species==focalSpecies])
  modelSummary$truePsiChange <- modelSummary$truePsiLast -modelSummary$truePsiFirst 
  modelSummary$truePsiOdds <- ((modelSummary$truePsiLast+e)/(1-modelSummary$truePsiLast+e))/((modelSummary$truePsiFirst+e)/(1-modelSummary$truePsiFirst+e))
  
  return(modelSummary)
}

getModelSummaries <- function(out){
  
  #get model results
  outCombined <- do.call(rbind,out)
  outCombined$Param[which(outCombined$Param=="psi.fs[10]")] <- "last.psi"
  outSummaries <- subset(outCombined,Param %in%c ("psi.mean","psi.odds","psi.last"))
  outSummaries <- ddply(outSummaries,.(Param),summarise,
                        lowerQ=quantile(mean,0.25),
                        medianQ=quantile(mean,0.5),
                        upperQ=quantile(mean,0.75))
  
  #get raw data
  Param <- c("raw.mean","raw.odds","raw.last")
  lowerQ <- c(quantile(outCombined$truePsiChange[outCombined$Param=="mean.psi.change"],0.25),
              quantile(outCombined$truePsiOdds[outCombined$Param=="odds.psi.change"],0.25),
              quantile(outCombined$truePsiLast[outCombined$Param=="last.psi"],0.25))
  
  medianQ <- c(quantile(outCombined$truePsiChange[outCombined$Param=="mean.psi.change"],0.5),
              quantile(outCombined$truePsiOdds[outCombined$Param=="odds.psi.change"],0.5),
              quantile(outCombined$truePsiLast[outCombined$Param=="last.psi"],0.5))
  
  upperQ <- c(quantile(outCombined$truePsiChange[outCombined$Param=="mean.psi.change"],0.75),
              quantile(outCombined$truePsiOdds[outCombined$Param=="odds.psi.change"],0.75),
              quantile(outCombined$truePsiLast[outCombined$Param=="last.psi"],0.75))
  
  rawDF <- data.frame(Param,lowerQ,medianQ,upperQ)
  
  return(rbind(outSummaries,rawDF))
}

rawdataAnalysis <- function(Obs){
  modelData <- getSpartaFormat(Obs,focalSpecies=focalSpecies,subsetPositive = TRUE)
  
  #glmer - mixed effect model
  e <- 0.00000001
  require(lme4)
  #transform into probabilities
  require(boot)
  
  lmer1 <- glmer(Obs ~ (factor(Year)-1)+(1|Site),data=subset(modelData,!is.na(Obs)),family=binomial)
  lmerPsiFirst <- inv.logit(summary(lmer1)$coefficients[1,1])
  lmerPsiLast <- inv.logit(summary(lmer1)$coefficients[(nrow(summary(lmer1)$coefficients)),1])
  lmerPsiChange <- lmerPsiLast -lmerPsiFirst 
  lmerPsiOdds <- ((lmerPsiLast+e)/(1-lmerPsiLast+e))/((lmerPsiFirst+e)/(1-lmerPsiFirst+e))
  temp <- data.frame(lmerPsiFirst,lmerPsiLast,lmerPsiChange,lmerPsiOdds)
  return(temp)
}

getrawSummaries <- function(outRaw){
  Param <- c("basic.mean","basic.odds","basic.last")
  lowerQ <- c(quantile(outRaw$lmerPsiChange,0.25),
              quantile(outRaw$lmerPsiOdds,0.25),
              quantile(outRaw$lmerPsiLast,0.25))
  
  medianQ <- c(quantile(outRaw$lmerPsiChange,0.5),
               quantile(outRaw$lmerPsiOdds,0.5),
               quantile(outRaw$lmerPsiLast,0.5))
  
  upperQ <- c(quantile(outRaw$lmerPsiChange,0.75),
              quantile(outRaw$lmerPsiOdds,0.75),
              quantile(outRaw$lmerPsiLast,0.75))
  
  data.frame(Param,lowerQ,medianQ,upperQ)
  
}


plotComparison <- function(modelSummaries,rawSummaries){
  
  temp <- rbind(modelSummaries,rawSummaries)
  
  temp$type <- as.character(sapply(temp$Param,function(x){
    
    if(grepl("odds",x)){ 
      "odds"
    }else if (grepl("last",x)){
        "last Occ"
    }else {
        "change"
      }}))
      
  temp$model <- c("OD","OD","OD","truth","truth","truth","naive","naive","naive")
  
  temp$model <- factor(temp$model, levels=c("truth","naive","OD"))
  temp$type <- factor(temp$type, levels=c("last Occ","change","odds"))
  
  print(ggplot(subset(temp,type!="odds"))+
          geom_crossbar(aes(x=model,y=medianQ,ymin=lowerQ,ymax=upperQ),width=0.4)+
          facet_wrap(~type,scales="free")+
          ylab("change in occupancy")+
          theme_bw())
  
  return(temp)
}





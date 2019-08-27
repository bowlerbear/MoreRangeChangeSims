getObsSummaries <- function(df){
require(plyr)
  
    
  #summaries per simulation
  perSim <- ldply(df,summarise,
                  #year1
                  nuOcc_year1 = mean(Occurence[Year==1]),
                  meanOcc_focal_year1 = mean(Occurence[Species==focalSpecies & Year==1]),
                  nuObs_year1 = mean(Visit1[Year==1]),
                  nuObs_focal_year1 = mean(Visit1[Species==focalSpecies & Year==1]),
                  #last year
                  nuOcc_last = mean(Occurence[Year==max(Year)]),
                  meanOcc_focal_last = mean(Occurence[Species==focalSpecies & Year==max(Year)]),
                  nuObs_last = mean(Visit1[Year==max(Year)]),
                  nuObs_focal_last = mean(Visit1[Species==focalSpecies & Year==max(Year)]))
  
                  #summary overall
                  colMeans(perSim)
  
}

getSims <- function(x){
  Occ <- getOccuHistory(occProb = OccProb, spatialbeta=Seffects, interactionbeta=Ieffects, spatialcovariates=Scovariate,temporalbeta=Teffects, temporalcovariates=Tcovariate)
  Obs <- getRepeatVisits(Occ=Occ, NVisits=NVisits, DetProb = DetProb, SiteDetEffects,YearDetEffects,IntDetEffects, Scovariate,Tcovariate)
  return(Obs)
}

getModels <- function(Obs,type="Standard"){
  modelData <- getSpartaFormat(Obs,focalSpecies=focalSpecies,subsetPositive = TRUE)
  bugs.data <- getBugsData(modelData)
  modelSummary <- runModel(bugs.data,modelData,focalSpecies)
  modelSummary$Type <- type
  
  #get actual change
  modelSummary$truePsiFirst <- mean(Obs$Occurence[Obs$Year==1 & Obs$Species==focalSpecies])
  modelSummary$truePsiLast <- mean(Obs$Occurence[Obs$Year==max(Obs$Year) & Obs$Species==focalSpecies])
  modelSummary$truePsiChange <- modelSummary$truePsiLast -modelSummary$truePsiFirst 
  modelSummary$truePsiOdds <- (modelSummary$truePsiLast/(1-modelSummary$truePsiLast))/(modelSummary$truePsiFirst/(1-modelSummary$truePsiFirst))
    
  return(modelSummary)
}

getModelSummaries <- function(out){
  
  outCombined <- do.call(rbind,out)
  outSummaries <- ddply(outCombined,.(Type,Param),
                        summarise,
                        meanValue = mean(mean),
                        medValue = median(mean))
  
  return(outSummaries)
}

plotSummaries <- function(out){
  
  outCombined <- do.call(rbind,out)
  require(ggplot2)
  ggplot(subset(outCombined,Param %in% c("prop.p","prop.muZ","trend.muZ")))+
    geom_violin(aes(x=Type,y=mean),draw_quantiles=c(0.25,0.5,0.75))+
    facet_wrap(~Param)+
    theme_bw()+ylab("Estimate")
}

pulsedVisits <- function(Obs, nStandardVisits,samplingbias=FALSE){
  
  #decide on atlas year
  atlasYear <- sample(unique(Obs$Year),1)
  
  #assume fewer visits outside these atlas years
  visits <- sample(which(grep("Visit",names(Obs))),nStandardVisits)
  
  Obs[!Obs$Year%in%atlasYear,visits]<-NA
  
  return(Obs)
}
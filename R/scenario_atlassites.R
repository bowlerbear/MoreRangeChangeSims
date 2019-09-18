reduceSites <- function(Obs, nStandardSites,samplingbias=FALSE){
  
  #decide on atlas year
  #atlasYear <- sample(unique(Obs$Year),1)
  atlasYear <- max(Obs$Year)-1
  
  #assume a random sample of sites are not sampling in all years except one
  if(samplingBias==FALSE){
    sitesSurveys <- sample(unique(Obs$Site),nStandardSites)#core survey sites
  }else if (samplingBias==TRUE){ 
    #assume a bias sample (low quality) are not sampled in all years except one
  sitesSurveys <- 1:NSites[order(Scovariate)][1:nStandardSites]
  }
  
  #apply bias
  Obs[((!Obs$Site %in% sitesSurveys) & (!Obs$Year %in% atlasYear)),
      which(grep("Visit",names(Obs)))] <- NA
  
  return(Obs)
}
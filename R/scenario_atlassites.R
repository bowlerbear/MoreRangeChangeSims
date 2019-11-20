reduceSites <- function(Obs, nStandardSites){
  
  #decide on atlas year
  #atlasYear <- sample(unique(Obs$Year),1)
  atlasYear <- max(Obs$Year)-1
  
  #assume a random sample of sites are not sampling in all years except one
  sitesSurveys <- sample(unique(Obs$Site),nStandardSites)#core survey sites
  #low for sampling bias or not???
  
  #apply bias
  Obs[((!Obs$Site %in% sitesSurveys) & (!Obs$Year %in% atlasYear)),
      grep("Visit",names(Obs))] <- NA
  
  return(Obs)
  
}
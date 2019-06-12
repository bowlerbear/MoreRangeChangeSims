#group together sites into pairs
groupSite <- function(NSites=NSites, Obs=Obs){
  siteTranslate <- data.frame(Site=1:NSites,newSites=rep(1:(NSites/2),each=2))
  
  #rename sites in the Obs data frame
  Obs$originalSite <- Obs$Site
  Obs$Site <- siteTranslate$newSites[match(Obs$Site,siteTranslate$Site)]
  
  return(Obs)
}

reduceVisits <- function(Obs=Obs_grouped){
  require(plyr)
  Obs_reduced<- ddply(Obs,.(Year,Site),function(x){
    siteLevels <- unique(x$originalSite)
    visitLevels <- grep('Visit',names(x))
    for(v in visitLevels){
      x[x$originalSite==sample(siteLevels,1),v] <- NA
    }
    #remove original Site column
    x <- x[,-grep("originalSite",names(x))]
    return(x)
  })
}
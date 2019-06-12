#apply it
applyReducedVisits <- function(Obs,prVisit=0.5,siteEffect=1.5,Scovariate){
  Obs_reduced <- Obs
  Obs_reduced$prVisited <- NA
  for(i in 1:nrow(Obs)){
      Obs_reduced$prVisted[i] <- inv.logit(logit(prVisit) + siteEffect*Scovariate[Obs_reduced$Site[i]])
  }
  
  Obs_reduced$SiteYear <- paste(Obs_reduced$Site,Obs_reduced$Year,sep="_")
  for(j in unique(Obs_reduced$SiteYear)){
    for(i in grep("Visit",names(Obs_reduced))){
      Obs_reduced[Obs_reduced$SiteYear==j,i][rbinom(1,1,min(Obs_reduced$prVisted[Obs_reduced$SiteYear==j]))==0] <- NA
    }
  }
  
  Obs_reduced <- Obs_reduced[,-which(names(Obs_reduced=="prVisted"))]
  Obs_reduced <- Obs_reduced[,-which(names(Obs_reduced=="SiteYear"))]
  
  return(Obs_reduced)
}



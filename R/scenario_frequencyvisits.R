#fewer visits to lower quality sites
applyReducedVisits <- function(Obs,type="static",prVisit=0.73,siteEffect=0,timeEffect=0){
  
  Obs_reduced <- Obs
  Obs_reduced$prVisited <- 1
  
  #add quality covariate
  Obs_reduced$Scovariate <- Scovariate[match(Obs_reduced$Site,1:length(Scovariate))]
  
  #site with highest quality is always visited
  
  #linear reduction
  if(type=="static"){
  for(i in 1:nrow(Obs_reduced)){
      Obs_reduced$prVisited[i] <- inv.logit(logit(prVisit) + 
                                              siteEffect*Scovariate[Obs_reduced$Site[i]])
  }
  }else if(type=="dynamic"){
    for(i in 1:nrow(Obs_reduced)){
      Obs_reduced$prVisited[i] <- inv.logit(logit(prVisit) + 
                                              siteEffect*Scovariate[Obs_reduced$Site[i]]+
                                              timeEffect*Obs_reduced$Year[i]*Scovariate[Obs_reduced$Site[i]])
    
    }
  }
  
  Obs_reduced$SiteYear <- paste(Obs_reduced$Site,Obs_reduced$Year,sep="_")
  
  for(j in unique(Obs_reduced$SiteYear)){
    for(i in grep("Visit",names(Obs_reduced))){
      Obs_reduced[Obs_reduced$SiteYear==j,i][rbinom(1,1,min(Obs_reduced$prVisited[Obs_reduced$SiteYear==j]))==0] <- NA
    }
  }

  Obs_reduced <- Obs_reduced[,-which(names(Obs_reduced)=="prVisited")]
  Obs_reduced <- Obs_reduced[,-which(names(Obs_reduced)=="Scovariate")]
  Obs_reduced <- Obs_reduced[,-which(names(Obs_reduced)=="SiteYear")]
  
  return(Obs_reduced)
  
}


#probability of a visit at a site is 50% greater when focal species recorded there in previous year
#biasEffect = first prob is prob of visit when species not seen, second prob is prob of visit when species seen
applyBiasedVisits <- function(Obs,biasEffect=c(0.33,0.66)){
  
  #for each year, identify whether the focal species was seen in previous year 
  presentData <- subset(Obs,(Visit1==1|Visit2==1|Visit3==1|Visit4==1|Visit5==1) &
                          Species==focalSpecies)
  #get first year seen at each site
  presentData <- ddply(presentData,.(Site),summarise,firstYear=min(Year))
  
  #add on to sheet whether species has been previously seen at each site
  Obs_reduced <- Obs
  Obs_reduced$firstYear <- presentData$firstYear[match(Obs_reduced$Site,presentData$Site)]
  
  Obs_reduced$SiteYear <- paste(Obs_reduced$Site,Obs_reduced$Year,sep="_")
  for(j in unique(Obs_reduced$SiteYear)){
    for(i in grep("Visit",names(Obs_reduced))){
      if(Obs_reduced$Year[Obs_reduced$SiteYear==j][1] <
         Obs_reduced$firstYear[Obs_reduced$SiteYear==j][1]){
        Obs_reduced[Obs_reduced$SiteYear==j,i][rbinom(1,1,biasEffect[1])==0] <- NA
      } else{
        Obs_reduced[Obs_reduced$SiteYear==j,i][rbinom(1,1,biasEffect[2])==0] <- NA
      }
    }
  }
  
  Obs_reduced <- Obs_reduced[,-which(names(Obs_reduced)=="firstYear")]
  Obs_reduced <- Obs_reduced[,-which(names(Obs_reduced)=="SiteYear")]
  
  return(Obs_reduced)
  
}


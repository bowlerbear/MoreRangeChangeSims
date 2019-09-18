pulsedVisits <- function(Obs, nStandardVisits,samplingbias=FALSE){
  
  #only in the atlas year do we have the maximum number of visits, 
  #otherwise, we just have the standard number of visits
  
  #decide on atlas year
  #atlasYear <- sample(unique(Obs$Year),1)
  atlasYear <- max(Obs$Year)-1
  
  #assume fewer visits outside these atlas years
  visits <- sample(grep("Visit",names(Obs)),nStandardVisits)
  
  Obs[!Obs$Year%in%atlasYear,visits] <- NA
  
  return(Obs)
}

pulsedSpreadVisits <- function(Obs,SCovariate,lowProb){
  
  #low quality sites are not visited in non-atlas years
  Obs$Scovariate <- Scovariate[match(Obs$Site,1:length(Scovariate))]
  
  #identify lowest 50% quality sites
  Obs$Low <- ifelse(Obs$Scovariate<quantile(Obs$Scovariate,0.50),1,0)
  
  #decide on atlas year
  #atlasYear <- sample(unique(Obs$Year),1)
  atlasYear <- max(Obs$Year)-1
  
  #if its not an atlas year, and its a low quality site, 
  #it has a lower probability of being visits
  #matrix of visit probabilities
  myRows <- nrow(subset(Obs,!Year%in%atlasYear & Low==1))
  myCols <- length(grep("Visit",names(Obs)))
  visited <- matrix(lowProb,nrow=myRows,ncol=myCols)
  visited <- apply(visited,1:2,function(x)rbinom(1,1,x))
  #set to NA when no visit
  visited[visited!=1] <- NA
  Obs[!Obs$Year%in%atlasYear & Obs$Low ==1,grep("Visit",names(Obs))] <- 
    Obs[!Obs$Year%in%atlasYear & Obs$Low ==1,grep("Visit",names(Obs))]*visited
  
  return(Obs)
}


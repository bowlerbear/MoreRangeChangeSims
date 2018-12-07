# This carries out some simple simulations and fitting


source("R/SimOcc.R")

ProbOcc <- seq(from=0.2, to=0.9, length=100)
ProbObs <- seq(from=0.9, to=0.2, length=100)
FirstSims <- SimOcc(PrOcc=ProbOcc, PrObs=ProbOcc, NVisits = 10)


# A workflow

#  1. Use N species, and S sites
NSpecies <- 4
NSites <- 8
ACovariate <- rnorm(NSites)


# 2. Simulate occupancy

getOccuHistory <- function(occProb,alpha,beta){
  
  alpha <- rnorm(NSpecies, 0, OccProb)#random variation among species in their intercept
  beta <- rnorm(NSpecies, 0.1,0.4) #random variation among species in their slopes

  #linear predictor to relate to predicted expected occurence
  lgtPrOcc <- alpha + ACovariate%o%beta
  PrOcc <- 1/(1 + exp(-lgtPrOcc))

  #perform binomial sampling
  Occ <- data.frame(apply(PrOcc, 2, function(pr) rbinom(length(pr), 1, pr)))
  
  #rearrange data frame
  rownames(Occ) <- paste0("Site", seq_along(PrOcc[,1]))
  colnames(Occ) <- paste0("Species", seq_along(PrOcc[1,]))
  
  #return to the data frame
  return(Occ)

}


# From the occupancy, simulate observations for different visits

#  (this is still being written)

getRepeatVisits <- function(nVisits=5, PrObs = DetProb){

    Occ$NVisits <- rep(nVisits, nrow(Occ))
    
    #simulate the visits
    Obs <- sapply(rownames(Occ), function(site, Occ, PrObs = DetProb) {
    nvisits <- Occ[site,"NVisits"]
    occ <- unlist(Occ[site, names(Occ)!="NVisits"])
    obs <- replicate(nvisits, rbinom(length(occ), 1, PrObs))
    obs <- t(obs*occ)
    colnames(obs) <- names(occ)
    obs
    }, Occ=Occ, simplify=FALSE)

    #unlist and reorganise the data frame    
    Obs <- plyr::ldply(Obs)
    names(Obs) <- gsub(".id", "Site", names(Obs))
    
    #Add Visit Number
    Obs$Visit <- rep(1:nVisits,NSpecies)
    
    #melt data frame
    library(reshape)
    Obs <- melt(Obs,id=c("Site","Visit"))
    names(Obs)[which(names(Obs)=="variable")] <- "Species"
    names(Obs)[which(names(Obs)=="value")] <- "obs"
    
    #Remove numbers from site and species names
    Obs$Site <- as.numeric(gsub("Site","",Obs$Site))
    Obs$Species <- as.numeric(gsub("Species","",Obs$Species))
    
    #return data frame
    return(Obs)
    
}

plotHistory <- function(obs){
  require(ggplot)
  ggplot(obs)+
    geom_point(aes(x=Visit,y=Species,color=factor(obs)),size=2)+
    scale_colour_manual(values=c("black","white"))+
    theme(legend.position="none")+
    facet_wrap(~Site)
}



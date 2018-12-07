# A workflow

#   This is the workflow we want to use. Below is a simple example
# Create parameters
# Simulate actual occupancy
# Simulate data
# Degrade data
# (Calculate diagnostics)
# Format data for sparta
# Fit model in sparta
# Summarise simulations




# This is ugly: we should make the whole thing an R package
sapply(dir("R/"), function(file) source(paste0("R/", file)))

###################
# Create parameters
#  (this needs cleaning up!)
#  1. Use N species, and S sites
ProbOcc <- seq(from=0.2, to=0.9, length=100)
ProbObs <- seq(from=0.9, to=0.2, length=100)
FirstSims <- SimOcc(PrOcc=ProbOcc, PrObs=ProbOcc, NVisits = 10)

NSpecies <- 4
NSites <- 8
ACovariate <- rnorm(NSites)
alpha <- rnorm(NSpecies, 0, 0.1)#random variation among species in their intercept
beta <- rnorm(NSpecies, 0.1,0.4) #random variation among species in their slopes

# Simulate actual occupancy
sapply(dir("R/"), function(file) source(paste0("R/", file)))
ActualOcc <- getOccuHistory(alpha=alpha, beta=beta, covariates=ACovariate)




# 2. Simulate occupancy


# Simulate data
# Degrade data
# (Calculate diagnostics)
# Format data for sparta
# Fit model in sparta
# Summarise simulations



#  (this is still being written)
Occ$NVisits <- rep(5, nrow(Occ))

# repeat this section for a number of years
nYear <- 5

All.Obs <- NULL

for (i in 1:nYear){
  
  Obs <- sapply(rownames(Occ), function(site, Occ, PrObs = 0.5) {
    nvisits <- Occ[site,"NVisits"]
    occ <- unlist(Occ[site, names(Occ)!="NVisits"])
    obs <- replicate(nvisits, rbinom(length(occ), 1, PrObs))
    obs <- t(obs*occ)
    colnames(obs) <- names(occ)
    obs
  }, Occ=Occ, simplify=FALSE)
  
  Obs <- plyr::ldply(Obs)
  
  # add a year column
  Obs$year <- i
  
  # rearrange
  Obs <- Obs[, c(1, ncol(Obs), 2:(ncol(Obs)-1))]
  names(Obs) <- gsub(".id", "Site", names(Obs))
  
  All.Obs <- rbind(All.Obs, Obs)
}

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


# small edits to match sims output with sparta input data - CO 06/12/2018

library(sparta)

# to run the model from visit level data, use the occDetFunc function in sparta
# function requires 2 data elements (usually created by the formatOccData function)

# data element 1: spp_vis - a dataframe with visit in first column
# and taxa for remaining columns, TRUE or FALSE depending om observations

# first, take the obs table and make the "site" column a visit (combination of
# site, date and year), this needs to be unique.  Currently combine row name and site name.

spp_vis <- All.Obs

# visit: name combined to represent "site_date_year"
spp_vis$Site <- paste(spp_vis$Site, "_", rownames(spp_vis), "_", spp_vis$year, sep = "")
colnames(spp_vis)[1] <- "visit"

# change 1s and 0s to TRUE and FALSE
spp_vis[, 3:ncol(spp_vis)] <- spp_vis[, 3:ncol(spp_vis)] == 1


# extract this now before removing the year column
occDetData <- spp_vis[, c("visit", "year")]

# remove year column
spp_vis <- spp_vis[, c(1, 3:ncol(spp_vis))]

# data element 2: occDetData - dataframe giving the site, list length and TP (time period)
occDetData <- cbind(occDetData, apply(X = spp_vis[, 3:ncol(spp_vis)], MARGIN = 1, FUN = sum))

colnames(occDetData) <- c("visit", "TP", "L")
  
occDetData$site <- sub("_.*", "", occDetData$visit)

# currently there are visits with no observations, don't have these in real datasets
# remove these
occDetData <- occDetData[!occDetData$L == 0, ]

# filter off these visits with 0 LL from the spp_vis
spp_vis <- spp_vis[spp_vis$visit %in% occDetData$visit, ]

# check that these fit into the sparta function

out <- occDetFunc(taxa_name = "Species1", 
                  occDetdata = occDetData, 
                  spp_vis = spp_vis, 
                  n_iterations = 10, 
                  nyr = 1, 
                  burnin = 2, 
                  thinning = 1, 
                  n_chains = 3, 
                  modeltype = c("ranwalk", "halfcauchy", "catlistlength"))

# this currently works 06/12/2018


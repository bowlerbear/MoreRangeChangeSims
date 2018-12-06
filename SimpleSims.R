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

alpha <- rnorm(NSpecies, 0,1)
beta <- rnorm(NSpecies, 0.1,0.4)

lgtPrOcc <- alpha + ACovariate%o%beta
PrOcc <- 1/(1 + exp(-lgtPrOcc))

Occ <- data.frame(apply(PrOcc, 2, function(pr) rbinom(length(pr), 1, pr)))
rownames(Occ) <- paste0("Site", seq_along(PrOcc[,1]))
colnames(Occ) <- paste0("Species", seq_along(PrOcc[1,]))

Occ[,1] <- Occ[,1]*0
# From the occupancy, simulate observations for different visits

#  (this is still being written)
Occ$NVisits <- rep(5, nrow(Occ))
Obs <- sapply(rownames(Occ), function(site, Occ, PrObs = 0.5) {
  nvisits <- Occ[site,"NVisits"]
  occ <- unlist(Occ[site, names(Occ)!="NVisits"])
  obs <- replicate(nvisits, rbinom(length(occ), 1, PrObs))
  obs <- t(obs*occ)
  colnames(obs) <- names(occ)
  obs
}, Occ=Occ, simplify=FALSE)

Obs <- plyr::ldply(Obs)
names(Obs) <- gsub(".id", "Site", names(Obs))



# small edits to match sims outoput with sparta input data - CO

library(sparta)

# to run the model from visit level data, use the occDetFunc function in sparta


# data element 1: spp_vis - a dataframe with visit in first column
# and taxa for remaining columns, TRUE or FALSE depending om observations

# first, take the obs table and make the "site" column a visit (combination of
# site and date), this needs to be unique.  Currently combine row name and site name.

spp_vis <- Obs

spp_vis$Site <- paste(spp_vis$Site, "_", rownames(spp_vis), sep = "")
colnames(spp_vis)[1] <- "visit"

# change 1s and 0s to TRUE and FALSE

spp_vis[, 2:ncol(spp_vis)] <- spp_vis[, 2:ncol(spp_vis)] == 1


# data element 2: occDetData - dataframe giving the site, list length and year

# calculate list length of each visit
occDetData <- spp_vis$visit

occDetData <- cbind(occDetData, apply(X = spp_vis[, 2:ncol(spp_vis)], MARGIN = 1, FUN = sum))

# need to add in the year

 colnames(occDetData) <- c("visit", "year", "LL")
  


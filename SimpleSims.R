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


# From the occupancy, simulate observations for different visits

#  (this is still being written)
Occ$NVisits <- rep(5, nrow(Occ))
Obs <- apply(Occ, 1, function(Occ, PrObs = 0.5) {
  Obs <- rep(Occ["NVisits"])
  
  
})





#
# Function to simulate simple occupancy
#  Inputs:
#    PrOcc - Probability of occupancy. If a scalar, all sites have same probability.
#            If a vector, simulates as many sites as length of PrOcc, with different
#            probabilities
#    PrObs - Probability of observation. If a scalar, all sites have same probability.
#            If a vector, simulates as many sites as length of PrOcc, with different
#            probabilities
#  NVisits - number of visits. Defaults to 4
#   NPlots - number of plots.  Ignored if length(PrOcc) > 1. Defaults to 100
# Outputs:
#  Data fraome with three columns:
#    Occupancy: occupancy (Bernoulli)
#    Observations: number of presences
#    Visits: Number of visits

SimOcc <- function(PrOcc, PrObs, NVisits = 4, NPlots = 100) {
  if(length(PrOcc)!=length(PrObs)) stop("PrOcc and PrObs should be same length")
  if(length(PrOcc)==1) PrOcc <- rep(PrOcc, NPlots)
  if(length(PrObs)==1) PrObs <- rep(PrObs, NPlots)
  if(length(NVisits)==1) NVisits <- rep(NVisits, NPlots)
  
  Occ <- rbinom(length(PrOcc), 1, PrOcc)
  Obs <- rbinom(length(PrOcc), NVisits, PrObs)*Occ
  data.frame(Occupancy=Occ, Observations=Obs, Visits=NVisits)
}



# This function simulated "true" occupancy
#   Arguments:
#   alpha - intercept (on logit scale) for occupancy
#   beta - slopes for covariate effects. Vector or data frame/matrix. 
#          Can be a scale, when the same effect is used over all species and covariates
#          If a vector, one per species
#          If a data frame, species in columns, covaraites in rows

#   covariates - vector or data frame/marix of covariate values
#  Returns a data frame of "true" occupancies with sites in rows and species in columns
getOccuHistory <- function(occProb=0, beta=0, covariates=NULL) {
  if(length(occProb)==1) warning("same intercept used for all species")
  if(!is.null(covariates)) {
    if(length(beta)==1) {
      warning("same slope used for all covariates and species")
      } else {
      if(is.vector(covariates)) { # if only one covariate
        if(length(beta)!=length(occProb)) stop("Need as many betas as alphas (= number of species)")
      } else {
        if(ncol(beta)!=length(occProb)) stop("Need as many beta columns as species")
        #if(nrow(beta)!=ncol(covariates)) stop("Need as many beta rows as covariates")
      }
    }
  }

  #linear predictor to relate to predicted expected occurence
  lgtPrOcc <- occProb + covariates%o%beta
  PrOcc <- 1/(1 + exp(-lgtPrOcc))
  
  #perform binomial sampling
  Occ <- data.frame(apply(PrOcc, 2, function(pr) rbinom(length(pr), 1, pr)))
  
  #rearrange data frame
  rownames(Occ) <- paste0("Site", seq_along(PrOcc[,1]))
  colnames(Occ) <- paste0("Species", seq_along(PrOcc[1,]))
  
  #return to the data frame
  return(Occ)
}

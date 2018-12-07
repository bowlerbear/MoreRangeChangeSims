# This function simulated "true" occupancy
#   Arguments:
#   alpha - intercept (on logit scale) for occupancy
#   beta - slopes for covariate effects
#   covariates - vector or data frame/marix of covariate values
#  Returns a data frame of "true" occupancies with sites in rows and species in columns
getOccuHistory <- function(alpha=0, beta=0, covariates=NULL) {
  if(length(alpha)==1) warning("same intercept used for all species")
  
    if(!is.null(covariates)) {
    if(!is.vector(covariates)){
      if(length(beta)==1) warning("same slope used for all covariates")
      if(length(beta)!=ncol(covariates)) stop("Need as many betas as covariates")
    } else {
      if(length(beta)>1) stop("Only one covariate, so must have only one beta")
    }
  }
  #linear predictor to relate to predicted expected occurence
  lgtPrOcc <- alpha + covariates%o%beta
  PrOcc <- 1/(1 + exp(-lgtPrOcc))
  
  #perform binomial sampling
  Occ <- data.frame(apply(PrOcc, 2, function(pr) rbinom(length(pr), 1, pr)))
  
  #rearrange data frame
  rownames(Occ) <- paste0("Site", seq_along(PrOcc[,1]))
  colnames(Occ) <- paste0("Species", seq_along(PrOcc[1,]))
  
  #return to the data frame
  return(Occ)
}

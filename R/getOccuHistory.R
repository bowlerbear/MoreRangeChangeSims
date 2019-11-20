
####single-time slice simulation##################################################
#
# This function simulated "true" occupancy
#   Arguments:
#   alpha - intercept (on logit scale) for occupancy
#   beta - slopes for covariate effects. Vector or data frame/matrix. 
#          Can be a scale, when the same effect is used over all species and covariates
#          If a vector, one per species
#          If a data frame, species in columns, covariates in rows

#   covariates - vector or data frame/marix of covariate values
#  Returns a data frame of "true" occupancies with sites in rows and species in columns
# getOccuHistory <- function(occProb=0, beta=0, covariates=NULL) {
#   if(length(occProb)==1) warning("same intercept used for all species")
#   if(!is.null(covariates)) {
#     if(length(beta)==1) {
#       warning("same slope used for all covariates and species")
#       } else {
#       if(is.vector(covariates)) { # if only one covariate
#         if(length(beta)!=length(occProb)) stop("Need as many betas as alphas (= number of species)")
#       } else {
#         if(ncol(beta)!=length(occProb)) stop("Need as many beta columns as species")
#         #if(nrow(beta)!=ncol(covariates)) stop("Need as many beta rows as covariates")
#       }
#     }
#   }
# 
#   #linear predictor to relate to predicted expected occurence
#   lgtPrOcc <- occProb + covariates%o%beta
#   PrOcc <- 1/(1 + exp(-lgtPrOcc))
#   
#   #perform binomial sampling
#   Occ <- data.frame(apply(PrOcc, 2, function(pr) rbinom(length(pr), 1, pr)))
#   
#   #rearrange data frame
#   rownames(Occ) <- paste0("Site", seq_along(PrOcc[,1]))
#   colnames(Occ) <- paste0("Species", seq_along(PrOcc[1,]))
#   
#   #return to the data frame
#   return(Occ)
# }


#####occupancy time series

getOccuHistory <- function(occProb, Seffects, Teffects, Ieffects,
                           Tcovariate,Scovariate) {
  
  # #warnings messages:
  # if(length(occProb)==1) warning("same intercept used for all species")
  # 
  # if(!is.null(spatialcovariates)) {
  #   if(length(spatialbeta)==1) {
  #     warning("same spatial slope used for all species")
  #   } 
  #   else {
  #     if(is.vector(spatialcovariates)) { # if only one covariate
  #       if(length(spatialbeta)!=length(occProb)) stop("Need as many betas as alphas (= number of species)")
  #     } else { #multiple covariates
  #       if(ncol(spatialbeta)!=length(occProb)) stop("Need as many beta columns as species")
  #     }
  #   }
  # }
  # 
  # if(!is.null(temporalcovariates)) {
  #   if(length(temporalbeta)==1) {
  #     warning("same temporal slope used for all species")
  #   } 
  #   else {
  #     if(is.vector(temporalcovariates)) { # if only one covariate
  #       if(length(temporalbeta)!=length(occProb)) stop("Need as many betas as alphas (= number of species)")
  #     } else {
  #       if(ncol(temporalbeta)!=length(occProb)) stop("Need as many beta columns as species")
  #     }
  #   }
  # }
  
  #linear predictor to relate to predicted expected occurence
  lgtPrOcc <- array(data=NA,dim=c(NSpecies,NSites,NYears))
  for(s in 1:NSpecies){
    for(j in 1:NSites){
      for(t in 1:NYears){
        
        #on logit scale
        lgtPrOcc[s,j,t] <- logit(OccProb[s]) + 
                            Scovariate[j]* Seffects[s] + 
                              Tcovariate[t] * Teffects[s] +
                              Scovariate[j] * Tcovariate[t] * Ieffects[s]   
        
      }
    }
  }
  
  #convert into probabilities
  library(boot)
  PrOcc <-inv.logit(lgtPrOcc)
  #hist(PrOcc)
  
  #perform binomial sampling
  Occ <- apply(PrOcc, c(1,2), function(pr) rbinom(length(pr), 1, pr))
  
  #rearrange data frame
  dimnames(Occ)[[1]] <- paste0("Year", 1:NYears)
  dimnames(Occ)[[2]] <- paste0("Species", 1:NSpecies)
  dimnames(Occ)[[3]] <- paste0("Site", 1:NSites)
  
  #return to the array
  return(Occ)
}

#####################################################################################

getAbundHistory <- function(lambda, Seffects, Teffects, Ieffects,
                            Tcovariate,Scovariate) {
  
  #linear predictor to relate to predicted expected occurence
  lgAbund <- array(data=NA,dim=c(NSpecies,NSites,NYears))
  for(s in 1:NSpecies){
    for(j in 1:NSites){
      for(t in 1:NYears){
        #on logit scale
        lgAbund[s,j,t] <- log(lambda[s]) + 
          Scovariate[j]*Seffects[s] + 
          Tcovariate[t]*Teffects[s] +
          Scovariate[j]*Tcovariate[t]*Ieffects[s]    
        
      }
    }
  }
  
  #convert into probabilities
  PrAbund <- exp(lgAbund)
  #hist(PrAbund)
  
  #perform binomial sampling
  Abund <- apply(PrAbund, c(1,2), function(pr) sapply(pr,function(x)rpois(1,x)))
  
  #rearrange data frame
  dimnames(Abund)[[1]] <- paste0("Year", 1:NYears)
  dimnames(Abund)[[2]] <- paste0("Species", 1:NSpecies)
  dimnames(Abund)[[3]] <- paste0("Site", 1:NSites)
  
  #return to the array
  return(Abund)
}
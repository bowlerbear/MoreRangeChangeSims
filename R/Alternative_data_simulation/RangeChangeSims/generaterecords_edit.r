generate_records <- function(nYrs=10, true_data, decline=0, which.decline=1, sites_to_use=c(1,0), 
                             pSVS=0.07, mv=20, vrs=F, stoch=T, nSites = 1000, drop = FALSE) {
  
  
  #get proportion of single visit sites
  if(length(pSVS)==1) pSVS <- c(pSVS[1], 0)
  
  zeroDecline <- function(i){
    records <- recording_cycle(pSVS=pSVS[1]+i*pSVS[2], true_data=true_data,
                               max_vis=mv, VisRichSites=vrs, stochastic=stoch)
    records$Year <- i
    return(records)
  }
  
  if (sum(decline)==0){
    records <- lapply(1:nYrs, zeroDecline)  #can start loop at zero
    
  } 
  
  if(drop == FALSE){
    # This part of the loop is triggered under two circumstances
    # 1: Testing power
    # 2: when 1 or more nonfocal (nuisance) species are declining (scenario F)
    # To test power, we simulate a declining population
    # For the IUCN criterion A, a decline of 30% over 10 years leads to a Vulnerable listing.
    
    # When nonfocal species are declining, which.decline specifies which ones
    if(length(which.decline) == 1) if(which.decline < 1) { # double if avoids a warning message
      #which.decline is a proportion of species, to be selected at random
      #convert the proportion into a number
      nsp_declining <- round(which.decline * (ncol(true_data)-1))
      # select that number at random from the nonfocal species
      which.decline <- sample(2:ncol(true_data), size=nsp_declining, replace=F)
      
      # is the focal also declining? (i.e. are we testing power under scenario F?)
      if(decline[1] == 0) # NO
        decline <- rep(decline[2], nsp_declining)
      else {#Yes
        decline <- c(decline[1], rep(decline[2], nsp_declining))
        which.decline <- c(1, which.decline)
      }
    }        
    if(length(which.decline) != length(decline)) stop('incompatible lengths')
    
    # CO edited: 09/11/2016: persitence rate changed to 1-nyrs to try to ensure that the decline is consistent
    # across the time period rather- yr 1 is the original occurrence level
    
    # each year, every record has a finite probability of persisting to the next year
    ann_persit_rate <- (1-decline)^(1/(nYrs-1))
    # for making populations go extinct, we need to iterate across a function nYr times
    # each time we pass in the extant sites and get back the ones that survived to the next year
    #annual_survival <- function(occ, rate) occ[rbinom(n=length(occ),1,rate)==1]
    #occ <- annual_survival(occ, ann_persit_rate) # I can't get this to work inside sapply, but does work inside afor loop
    records <- list()
    
    # new 25/7/16: count the true number of sites!
    true_occ <- matrix(data=NA, nrow=length(which.decline), ncol=nYrs)
    
    # CO added: 09/11/2016: first year want no decline
    
    #for(i == 1) {  #can start loop at zero
    #  for(j in 1:length(which.decline)){
    #  occ <- which(true_data[,which.decline[j]]==1) #site numbers where focal is present
    #  extinctions <- occ[rbinom(n=length(occ),1,ann_persit_rate[j])==0] # the site numbers at which extinctions will happen this year
    #  true_data[extinctions,which.decline[j]] <- 0 # set them to extinct
    # }
    
    records[[1]] <- recording_cycle(pSVS=pSVS[1]+1*pSVS[2], true_data=true_data, max_vis=mv, VisRichSites=vrs, stochastic=stoch) #if starting loop at zero, use records[[i+1]]
    records[[1]]$Year <- 1
    true_occ[,1] <- sum(true_data[,'focal'])/nSites# new 25/7/16: proportion of occ sites
    #}
    #records <- do.call(rbind, records)
    #}
    
    for(i in 2:nYrs) {  #can start loop at zero
      for(j in 1:length(which.decline)){
        occ <- which(true_data[,which.decline[j]]==1) #site numbers where focal is present
        extinctions <- occ[rbinom(n=length(occ),1,ann_persit_rate[j])==0] # the site numbers at which extinctions will happen this year
        true_data[extinctions,which.decline[j]] <- 0 # set them to extinct
      }
      
      records[[i]] <- recording_cycle(pSVS=pSVS[1]+i*pSVS[2], true_data=true_data, max_vis=mv, VisRichSites=vrs, stochastic=stoch) #if starting loop at zero, use records[[i+1]]
      records[[i]]$Year <- i
      true_occ[,i] <- sum(true_data[,'focal'])/nSites# new 25/7/16: proportion of occ sites
    }
    #records <- do.call(rbind, records)
  }# end of drop is false loop
  #records <- lapply(1:nYrs, FUN=recording_cycle, nVisits=nVisits, true_data=true_data) #simple version with constant number of visits
  # replace melt with rbind 08/05/2014 TA&GP
  # records <- melt(records, id.vars=1:3) #simply appends the two list elements into a simple 
  # names(records)[4] <- 'Year'
  
  
  if(drop == TRUE){
    records <- list()
    true_occ <- matrix(data=NA, nrow=length('focal'), ncol=nYrs)
    
    
    for(i in 1:4) {  #can start loop at zero
      decline = 0
      ann_persit_rate <- (1-decline)^(1/nYrs)
      occ <- which(true_data[,'focal']==1) #site numbers where focal is present
      extinctions <- occ[rbinom(n=length(occ),1,ann_persit_rate)==0] # the site numbers at which extinctions will happen this year
      true_data[extinctions,'focal'] <- 0 # set them to extinct
      
      records[[i]] <- recording_cycle(pSVS=pSVS[1]+i*pSVS[2], true_data=true_data, max_vis=mv, VisRichSites=vrs, stochastic=stoch) #if starting loop at zero, use records[[i+1]]
      records[[i]]$Year <- i
      true_occ[,i] <- sum(true_data[,'focal'])/nSites# new 25/7/16: proportion of occ sites
      
    }
    
    
    
    for(i in 5:10) {  #can start loop at zero
      decline = 0.6
      ann_persit_rate <- (1-decline)^(1/1)
      
      # dropping year
      occ <- which(true_data[,'focal']==1) #site numbers where focal is present
      extinctions <- occ[rbinom(n=length(occ),1,ann_persit_rate)==0] # the site numbers at which extinctions will happen this year
      true_data5 <- true_data
      true_data5[extinctions,'focal'] <- 0
      
      records[[i]] <- recording_cycle(pSVS=pSVS[1]+5*pSVS[2], true_data=true_data5, max_vis=mv, VisRichSites=vrs, stochastic=stoch) #if starting loop at zero, use records[[i+1]]
      records[[i]]$Year <- i
      true_occ[,i] <- sum(true_data5[,'focal'])/nSites# new 25/7/16: proportion of occ sites
      
      
    }
    
    
  }
  
  
  records <- do.call(rbind, records)
  attr(records, "true_occ") <- true_occ # new 25/7/16
  return(records)
  
}
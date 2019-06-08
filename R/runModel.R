runModel <- function(bugs.data,modelData,focalSpecies){
  
  #fit model
  require(rjags)
  require(R2WinBUGS)
  require(jagsUI)
  
  #define model params
  ni <- 10000   ;   nb <- 2000   ;   nt <- 5   ;   nc <- 3
  
  #specify parameters to monitor
  params <- c("psi.fs","mean.psi","mean.psi.change")
  
  require(reshape2)
  zst <- acast(subset(modelData,Species==focalSpecies), Site~Year, value.var="Obs",fun=max)
  zst [is.infinite(zst)] <- 0
  inits <- function(){list(z = zst)}
  
  #run model
  out <- jags(bugs.data, inits=inits, params, "R/basic_odm.txt", n.thin=nt,
              n.chains=nc, n.burnin=nb,n.iter=ni,parallel=T)
  
  outSummary <- data.frame(out$summary)
  outSummary$Param <-row.names(out$summary)
  
  #get summary parameters
  return(outSummary)
  
}
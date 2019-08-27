runModel <- function(bugs.data,modelData,focalSpecies){
  
  #fit model
  require(rjags)
  require(R2WinBUGS)
  require(jagsUI)
  
  #define model params
  ni <- 10000   ;   nb <- 2000   ;   nt <- 5   ;   nc <- 3
  
  #specify parameters to monitor
  params <- c("prop.muZ","prop.p","first.psi","mean.psi","mean.p",
              "mean.psi.change","odds.psi.change",
              "psi.fs")
  
  require(reshape2)
  zst <- acast(subset(modelData,Species==focalSpecies), Site~Year, value.var="Obs",fun=max,na.rm=T)
  zst [is.infinite(zst)] <- 0
  inits <- function(){list(z = zst)}
  
  #run model
  out <- jags(bugs.data, inits=inits, params, "R/basic_odm.txt", n.thin=nt,
              n.chains=nc, n.burnin=nb,n.iter=ni,parallel=T,verbose=FALSE)
  
  outSummary <- data.frame(out$summary)
  outSummary$Param <-row.names(out$summary)
  
  #get summary parameters
  return(outSummary)
  
}
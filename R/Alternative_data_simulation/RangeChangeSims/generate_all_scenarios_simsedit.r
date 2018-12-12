generate_all_scenarios <- function(nSites=1000, nSpecies=25, nYrs=10, pSVS=0.07, mv=10, vrs=F, stoch=T, Scenarios='BCDF', combos=F,
                                   save_data=F, id='', pFocal=list(Occ=0.5, DetP=0.5),p_short=list(init=0.6,final=0.9), pDetMod=0.2, decline=0, drop = FALSE){
  # This is a wrapper for creating data
  # 25 October: added 'decline' for testing power
  #  3 December: removed B1 & replaced E
  # 25 January 2013: I added B3: new sites are biased toward the focal species
  #  4 February: added F: decline in a nonfocal species
  #  7 February: added the option of user-controlled subset of scenarios
  # 13 February: I added the optional argument pFocal to allow exploration of changing it from outside
  # 14 February: I dropped D1 and F1
  # 8 April: removed B3n
  
  # the actual data   
  true_data <- create_data(nSites=nSites, nSpecies=nSpecies, pFocal=pFocal)
  
  if(grepl('A', Scenarios)){ ### CO added this 04/11/2016, only want even recording scenario
  #cntrlsp <- median_occ(true_data) #the identity of a control species
  records_A = generate_records(nYrs=nYrs, pSVS=pSVS, true_data=true_data, decline=decline, mv=mv, vrs=vrs, stoch=stoch, drop=drop)
  records <- list(A_EvenRcrdng=records_A)
  }
  
  if(grepl('B', Scenarios)){ # run the B Scenarios
    #records$B1_IncrnSite = subset_recs_by_X(records_A, X=records_A$Site, init=0.5), # first year has 50% sites, final year has all
    records$B2_IncrnVisit = subset_recs_by_X(records_A, X=records_A$Visit, init=0.5) # first year has 50% visits, final year has all
    records$B3f_IncrnVBiasFc = subset_recs_Biased(records_A, init=0.5, target=true_data[,1])
    #records$B3n_IncrnVBiasNf = subset_recs_Biased(records_A, init=0.5, target=true_data[,cntrlsp]) # bias wrt a species of median detectability      
  }
  if(grepl('C', Scenarios)){ # run the C Scenarios
    records$C1_pShortLEven = shorten_lists(records_A, p_short=mean(unlist(p_short))) #a fixed proportion will be short, split evenly between incidental and others 
    records$C2_pShortLIncr = shorten_lists(records_A, p_short=p_short) #proportion of short lists increases over time
  }
  if(grepl('D', Scenarios)){ # run the D Scenarios
    #records$D1_SelectvEven = bias_focal(records_A, disadvantage=pDetMod) #Remove records of the focal species at random (bias *against* detection) 
    records$D2_SelectvIncr = bias_focal(records_A, disadvantage=c(pDetMod,0)) #Bias against detection decreases to zero (i.e. increasing apparency)
  }  
  if(grepl('E', Scenarios)){ # generally not used - after 25/4 it is unuseable since we've dropped nVisits as a parameter
    #records$E1n_IncidentalNfEv = rbind(records_A, incidental_records(true_data, nYrs=nYrs, nVisits=nVisits, p_short=mean(unlist(p_short)), focal=F)) #add nonfocal records
    #records$E2n_IncidentalNfIn = rbind(records_A, incidental_records(true_data, nYrs=nYrs, nVisits=nVisits, p_short=p_short, focal=F)) #add nonfocal records
    #records$E1f_IncidentalFcEv = rbind(records_A, incidental_records(true_data, nYrs=nYrs, nVisits=nVisits, p_short=mean(unlist(p_short)), focal=T)) #add nonfocal records
    #records$E2f_IncidentalFcIn = rbind(records_A, incidental_records(true_data, nYrs=nYrs, nVisits=nVisits, p_short=p_short, focal=T)) #add nonfocal records
  }
  if(grepl('F', Scenarios)){
    # Added 4/2/13: A new scenario in which a nonfocal species (indexed by cntrlsp) declines dramatically (50% over the period)
    #records$F_NfDecline <- generate_records(nYrs, pSVS, true_data, decline=c(decline,0.5), which.decline=c(1,cntrlsp))
    #records$F_NfDecline10 <- generate_records(nYrs, pSVS, true_data, decline=c(decline,0.5), which.decline=0.1)
    records$F_NfDecline <- generate_records(nYrs=nYrs, pSVS=pSVS, true_data=true_data, decline=c(decline,0.3), which.decline=0.5, mv=mv, vrs=vrs, stoch=stoch) # half the species are VU
  }
  if(combos){# now add the combination scenarios
    records$B2C1_IncVsShtEv = shorten_lists(records$B2_IncrnVisit, p_short=mean(unlist(p_short)))
    records$B2C2_IncVsShtIn = shorten_lists(records$B2_IncrnVisit, p_short=p_short)
    #records$B2D1_IncVsSelEv = bias_focal(records$B2_IncrnVisit, disadvantage=pDetMod)
    #records$C1D1_ShtEvSelEv = bias_focal(records$C1_pShortLEven, disadvantage=pDetMod)
    #records$C2D1_ShtInSelEv = bias_focal(records$C2_pShortLIncr, disadvantage=pDetMod)
    #records$B2C1D1_IncVShtSelEv = bias_focal(records$B2C1_IncVsShtEv,disadvantage=pDetMod) #should kill LL model (new 6/11)
  }
  # for validity, we there's an alternative version of B2, done differently (for reality check). Doesn't work for power test
  #if(decline==0) records$B2_IncrnVisit_a = generate_records(nYrs, pSVS=c(pSVS/2,pSVS/(2*nYrs)), true_data)
  
  attr(records,"true_data") <- true_data
  attr(records,"simpars") <- list(nSpecies=nSpecies,nSites=nSites,nYrs=nYrs,pSVS=pSVS,pFocal=pFocal,p_short=p_short,
                                  pDetMod=pDetMod,decline=decline,visit_rich_sites=vrs, stochastic=stoch, max_vis=mv)
  if(save_data) save_records(records, id)
  
  return(records)
}

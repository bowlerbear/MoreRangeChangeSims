############################################################################
# Script - to run generate_datasets_function.r to get n datasets           #
# with the same paratmeters to use to test the models with simulated data  #
############################################################################

rm(list = ls())

# set working directory to where paper functions etc saved
scripts_dir<- ('C:/Users/db40fysa/Nextcloud/sMon-Analyses/MoreRangeChangeSims/R/Alternative_data_simulation/RangeChangeSims')

# output directory
output_dir <- 'C:/Users/db40fysa/Nextcloud/sMon-Analyses/Scenario_data'

# set desired criteria
rec_int <- "Low"  # set the recording intensity (see values below)
decline <- 0.5    # set a decline for the focal species

# scenario details
scenario <- paste(rec_int, "_", decline, sep = "")

# create a directory to save these datasets in for this scenario
dir.create(paste(output_dir, scenario, sep = "/"))

# name this directory
dataset_dir <- paste(output_dir, "/", scenario, sep = "")

# set additional parameters
drop <- FALSE  # This parameter determines where there should be a sudden drop in species occupancy 
               # (rather than a gradual decline), this was for used to test the model

# These set the parameter that determines the recording intensity
# how these values are derived is described in Nick's sims paper (based on UK datasets)
if(rec_int == "Medium") pSVS <- 0.07
if(rec_int == "Low") pSVS <- 0.05
if(rec_int == "Super_Low") pSVS <- 0.01


# set requirements of the outputted dataset
nspecies <- 25
nsites <- 1000
nyears <- 40


#define occurence probability
Occ = 0.2 

#define detection probability
DetP = 0.5


n <- 5 # number of datasets required


# source in all the functions required, these are taken from Nick's sims paper
source(paste0(scripts_dir, "/Sim_functions.r"))

# These are also from Nick's sims paper, but I have edited these to record true occupancy over
# time and to generate a dataset using the control scenario only.  Framework for additional 
# scenarios are still in there so can be adapted for our needs if required.
source(paste0(scripts_dir, "/generaterecords_edit.r")) # edited function so that true occupancy over time is produced as an attribute
source(paste0(scripts_dir, "/generate_all_scenarios_simsedit.r"))

# creates n datasets
for(i in 1:n){
  
  ## now setting occ and det probabilities of the focal species to 0.2
  # reviewer comments, need to make dataset more comparable to the ants
  
  # here occupancy and detection set to 0.2.  This generated low recording intensity datasets similar to the 
  # UK ant data which I was comparing things to in my model testing paper
  
  recs <- generate_all_scenarios(nSites = nsites, nSpecies = nspecies,
                                 pFocal = list(Occ = Occ, DetP = DetP), 
                                 nYrs = nyears, pSVS = pSVS,
                                 mv = 10, Scenarios = "A", 
                                 p_short = list(init = 0.6, final = 0.9),
                                 pDetMod = 0.2, decline = decline, drop = drop)
  
  save(recs, file = paste(dataset_dir, "/", rec_int, "int_records_", decline, "_", i, ".rdata", sep = ""))
  
}


head(recs$A_EvenRcrdng)  # recs$A_EvenRcrding is the dataset producted from the input parameters and the even recording scenario (A)

###Simulation Study for occupancy models

This repo contains the code for simulations to explore the robustness of occupancy detection models to different types of human data collection behaviours.

Specifically we aim to:
- assess the effect of different biases on occupancy trend estimates
- develop diagnostic for assessing whether a given bias is likely present within a given real-world dataset
- extend the occupancy-models to account for the bias
- identify additional pieces of metadata that are necessary to fully account for the bias in the model.

###Simulation steps
The core model file is called "sims_workflow.Rmd".

The steps are the analysis are:
(1) generate an idealized community based on a known set of parameters
  - parameters are input in a csv sheet
  - apply the getSims() function, which is a wrapper for getParams(), getOccHistory and
    getRepeatVisits()
    
(2) degrade the dataset according to scenarios below
  - these are coded in the files "scenario_...."
  
(3) fit a standard occupancy model and compare true parameters in (1) with those estimated by     the model. We also compare a simple naive analysis that does not account for imperfect        detection. 
    - this is done by applying the function getModels, which is a wrapper function for            getSpartaFormat(), getBugsData() and runModel()
    
(4) Various functions for summarising and plotting the data are found in the helperFunctions.R script
  Specially, we compare the ability of the occupancy detection model to:
    - estimate true occurrence proportion of a focal species
    - estimate true occurrence change (difference between first and last year) of a focal           species
    
###Scenarios under consideration so far (additional ideas welcome)
(a) Atlas schemes
  - extension of grids (into lower quality areas)
  - greater effort (higher detection probability, more visits)
  
(b) Habiat quality effects
  - decline in the frequency of visits with habitat degradation
  - more visits initially in more higher quality sites
  
(c) Model mispecification
  - spatial mismatch between scale of dynamics and scale of "grid" in analysis
  - heteroegeneity of detection probability
  - auto-correlation among consecutive visits (trap shyness)

(d) Species/Community properties
  - regional pool size
  - specialist species (of species-poor vs species-rich sites)
# MoreRangeChangeSims

### Simulation Study for occupancy models

This repo contains the code for simulations to explore occupancy model mis-fits in more detail.
The aim is to explore the consequences of different types of sampling behaviour (by citizen scientists) for occupancy-detection models. Specifically we aim to:
- assess the effect of different biases on occupancy estimates estimated by the standard occupancy-detection model (mostly along the Sparta formulation).
- develop diagnostic for assessing whether a given bias is likely present within a given dataset
- extend the occupancy-models to account for the bias
- identify additional pieces of metadata that are necessary to fully account for the bias in the model.


#Scenarios
(1) Observer behaviour: atlas schemes
  - extension of grids (into lower quality areas)
  - greater effort (higher detection probability, more visits)
  
(2) Observer behaviour: habiatt quality effects
  - decline in the frequency of visits with habitat degradation
  - more visits initially in more higher quality sites
  
(3) Model mispecification
  - spatial mismatch
  - heteroegeneity of detection probability
  - auto-correlation among consecutive visits (trap shyness)

(4) Species/Community properties
  - regional pool size
  - specialist species (of species-poor vs species-rich sites)
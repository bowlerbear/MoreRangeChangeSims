# This carries out some simple simulations and fitting


source("R/SimOcc.R")

ProbOcc <- seq(from=0.2, to=0.9, length=100)
ProbObs <- seq(from=0.9, to=0.2, length=100)

FirstSims <- SimOcc(PrOcc=ProbOcc, PrObs=ProbOcc, NVisits = 10)

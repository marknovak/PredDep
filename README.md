# PredDep
Simulation code for "Quantifying predator dependence in the functional response of generalist predators."

Start with the "....-SimRunMe.R" script.  This file sources the other scripts, which:
(1) simulate some data ("...-SimData.R"),
(2) load the likelihood functions for the alternative functional response models ("...-Models.R"),
(3) use the analytical expressions to obtain starting values ("...-StartVals.R"), and
(4) fit the models ("...-FitModels.R").

The rest of the "...-SimRunMe.R" script just plots some simple figures for comparison purposes..

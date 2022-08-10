# MECFeedbackInTcWorld

This repository contains the model files to reproduce the results in: The Moisture-Entrainment-Convection feedback can be sufficient to cause Spontaneous Tropical Cyclone Genesis by Argel Ramirez Reyes and Da Yang (in revision for GRL).

## Simulation setup
* The model is the system for atmospheric modeling, SAM, version 6.10.6 and it can be found in /model. 

To reproduce the results, you must first replace the files found in model/SAM6.10.6/SRC with the ones found in model/modified_src. These files contain my additions and modifications.
1. buoyancy.f90 was modified by Da Yang to be able to turn off the buoyancy effect of water vapor
1. domain.f90 is modified to create our domain
1. nudging.f90 was modified to be able to relax the water vapor field in the clear sky

Additionally, model/SAM6.10.6/SRC/SGS_TKE_ARR contains the modified subgrid scale to disable the mixing of water vapor.

* model/initial_soundings contain the initial condition used for this study
* model/namelists contain the scripts used to create the namelists and launch simulations


## Analysis tools

* src/ Contains julia source files to detect TCs and to compute quantities as inflow and updraft ratio

* notebooks contain plotting scripts to reproduce the figures

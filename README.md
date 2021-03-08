# AMPAR-Trafficking-Model

AMPAR-Trafficking-Model contains the Python code used in the paper "The biophysical basis underlying the maintenance of early phase long-term potentiation" (BioRxiv: https://doi.org/10.1101/2020.08.12.247726). Code documentation can be found under https://moritzb90.github.io/AMPAR-Trafficking-Model/.

* A package "ampartrafficking" is provided containing four modules:
  - rate_model: Contains functions and classes for the rate model of AMPAR-trafficking (inlcuding the mean-field approximation of cooperative binding/unbinding rates).
  - stochastic_model: Contains functions and classes for the stochastic, spatial model of cooperative receptor binding.
  - parameter_sampling: Contains functions for random parameter sampling and subsequent simulation of the model. Depends on rate_model.
  - frap: Contains functions to carry out fluorescence recovery after photobleaching (FRAP) simulations. Depends on stochastic_model.
* The folders "Mean-Field-Model" and "Stochastic-Model" further contain python scripts to reproduce the figures shown in the paper. These scripts import modules from the ampartrafficking package. Each script is named after the respective figure as appearing in the Paper.
* The folders sphinx and docs contain documentation files.

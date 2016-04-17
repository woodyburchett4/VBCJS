# VBCJS

This repository contains a couple different variational Bayesian algorithms for fitting the Cormack-Jolly-Seber model with individual, time-varying, continuous covariates, as well as an MCMC implementation of the same model.

Each VBCJS algorithm has it's own subdirectory.  The simulations.R file will run each algorithm using simulated data and visually compare them to the MCMC implementation.  The analyze_csw_subset.R file will do the same on the first 10 years of cliff swallows data.
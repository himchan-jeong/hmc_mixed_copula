# Codes for `Fully and empirical Bayes approaches to estimating copula-based models for bivariate mixed outcomes using Hamiltonian Monte Carlo'

Note that supplementary codes are sourced to run the main codes.

## Main codes

**sim.r**: Code for simulation and Full/Empirical Bayes analysis of data

**burn_analysis.r**: Code for the analysis of burn injury data (data can be accessed via the link in the code)

## Supplementary codes

**copulalinks.r**: different transform of copula parameter for each copula family

**corrmat.r**: incorporating possible correlation in the covariates

**decodeCombo.r**: decode combination code to get simulation settings

**model_selection.r**: calculate DIC and LPML for model selection

**simulate_data.r**: source code for simulation of covariates and bivariate responses

**EB.stan**: RStan model used for empirical Bayesian analysis

**FB.stan**: RStan model used for full Bayesian analysis

**Cppfns.cpp**: copula functions written in C++ for effective optimization

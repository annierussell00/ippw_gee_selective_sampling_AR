# ippw_gee_selective_sampling_AR

This repository contains the R code used for the two simulation studies associated with research in "Generalising evidence under exposure-history-dependent sampling with an application in psoriatic arthritis".

The simulations illustrate the use of inverse probability of participation weighted generalised estimating equations (IPPW-GEE) under exposure-history-dependent selective sampling.

Contents
01_simulation_A.R — Code for Simulation A.
02_simulation_B.R — Code for Simulation B.

Software

The analyses were conducted using R version 4.4.1.

Data

Simulation dataset A is generated within the supplied R code.
Simulation_B_data.csv — simulated population data used to define the target estimand.
SEQUENCE study patient data not included in this repository due to patient confidentiality.

Reproducibility

The code includes the data-generating mechanisms, selective sampling procedures, calculation of participation weights, model fitting, and evaluation of simulation performance.


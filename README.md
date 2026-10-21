# ippw_gee_selective_sampling_AR

This repository contains the R code used for the two simulation studies associated with research in "Correcting selection bias in sequential treatment studies: population-standardised estimation using weighted generalised estimating equations".

The simulations illustrate the use of inverse probability of participation weighted generalised estimating equations (IPPW-GEE) under exposure-history-dependent selective sampling.

Contents
Simulation_A.R — Code for Simulation A.
Simulation_B.R — Code for Simulation B.
Simulation_B_Data_Generation is included for provenance. It documents how Simulation_B_data.csv was generated from restricted RNHRD/SEQUENCE inputs — but cannot be executed without access to those inputs. Simulation_B.R is the main reproducible analysis and only needs Simulation_B_data.csv and Sample_data.csv as inputs
Simulation_B_data.csv - Generated data used in simulation B

Software

The analyses were conducted using R version 4.4.1.

Data

Simulation dataset A is generated within the supplied R code.
Simulation_B_data.csv — synthetic population data used to define the target estimand. Elements of the data is callibrated to the SEQUENCE study data set.
The SEQUENCE study and RNHRD data are not included in this repository due to patient confidentiality.

Reproducibility

The code includes the data-generating mechanisms, selective sampling procedures, calculation of participation weights, model fitting, and evaluation of simulation performance.


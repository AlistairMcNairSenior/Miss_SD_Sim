
This repo was created by AM Senior to contain R code for the simulation described in Nakagawa et al. A practical and readily implementable method for handling missing standard deviations in the meta-analysis of response ratios.

There are 6 files containing code for the simulation and two .Rdata file containing results from the simulation. In order of use the files are:

Header.R - contains 2 functions used for the simulation.

MLMA_HPC.R - code to run the simulation for MLMA. Designed to run as an array (1000 instances) job on the HPC here in Sydney.

REMA_HPC.R - code to run the simulation for REMA. Designed to run as an array (1000 instances) job on the HPC here in Sydney.

Aggregate_MLMA.R - code to aggregate the results from the individual instances of the MLMA simulation on the array, calculate the biases from the raw data, and calculate the summary stats for each set of parameters. Produces the results data file 'agg_results_mlma.Rdata'. Note the raw data is not kept in the GitHub repo as there is so much of it

Aggregate_REMA.R - code to aggregate the results from the individual instances of the REMA simulation on the array, calculate the biases from the raw data, and calculate the summary stats for each set of parameters. Produces the results data file 'agg_results_rema.Rdata'. Note the raw data is not kept in the GitHub repo as there is so much of it

Plot_MLMA.R - plots the results of the mlma simulation (agg_results_mlma.Rdata).

Plot_REMA.R - plots the results of the rema simulation (agg_results_rema.Rdata). 


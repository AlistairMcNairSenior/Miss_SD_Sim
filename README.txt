
This repo was created by AM Senior to contain R code for the simulation described in Nakagawa et al. A practical and readily implementable method for handling missing standard deviations in the meta-analysis of response ratios.

There are 5 files containing code and one .Rdata file containing results from the simulation. In order of use the files are:

Header.R - contains 2 functions used for the simulation.

MLMA_HPC.R - code to run the simulation for MLMA. Designed to run as an array (1000 instances) job on the HPC here in Sydney.

REMA_HPC.R - code to run the simulation for REMA. Designed to run as an array (1000 instances) job on the HPC here in Sydney.

Aggregate.R - code to aggregate the results from the individual instances on the array, calculate the biases from the raw data, and calculate the summary stats for each set of parameters. Produces the results data file 'agg_results.Rdata'.

Plot_MLMA.R - plots the results of the simulation.
 

# Simulation to Support Nakagawa et al. A practical and readily implementable method for handling missing standard deviation in the meta-analysis of response ratios. 

# This script was used to aggregate multiple instances from the array job on the HPC. It also then calculates the bias, coverage and aggregate statsitics used in the plotting

# Clean up the R Environment 
rm(list=ls())

# Set the working directory
wd<-"/Users/alistairsenior/OneDrive - The University of Sydney (Staff)/Nakagawa_Ecology_Missing_SD" # Note this directory is not in the git repo. Contains a massive amount of raw data. Hence this file will not run from the repo.
setwd(wd)

# Load the relevant libraries, and header
library(metafor)
library(plyr)
library(doSNOW)

# Check the files
files<-dir("HPC_data/Miss_Sim_REMA")[-c(1:5)]

load(paste0("HPC_data/Miss_Sim_REMA/", files[1]))
res_unlist<-results
for(i in 2:length(files)){
	load(paste0("HPC_data/Miss_Sim_REMA/", files[i]))
	for(j in 1:length(res_unlist)){
		res_unlist[[j]]<-rbind(res_unlist[[j]], results[[j]])
	}
}
results<-res_unlist

# Load the parameters 
load("HPC_data/Miss_Sim_REMA/Parameters_REMA.Rdata")

# Long format all the results
for(i in 1:length(results)){
	long_i<-cbind(results[[i]], parameters[i,])
	if(i == 1){
		long<-long_i
	}else{
		long<-rbind(long, long_i)	
	}
}
long_mat<-as.matrix(long)

# Reformat from long to wide for the methods
long_full<-as.data.frame(rbind(long_mat[,c(1:3, 16:27)], long_mat[,c(4:6, 16:27)], long_mat[,c(7:9, 16:27)], long_mat[,c(10:12, 16:27)], long_mat[,c(13:15, 16:27)]))
long_full$Method<-c(rep("Complete Data", dim(long_mat)[1]), rep("Method 1.1", dim(long_mat)[1]), rep("Method 1.2", dim(long_mat)[1]), rep("Method 2", dim(long_mat)[1]), rep("Method 3", dim(long_mat)[1]))

# Calculate the bias and coverage
long_full$bias<-long_full$whole_ests - long_full$lnRR
dfs<-long_full$k_study-1
tag<-which(long_full$icc_study == 0)
dfs[tag]<-long_full$k_study[tag] * long_full$k_effect_mu[tag] - 1
ts<-qt(0.975, df=dfs)
long_full$coverage<-(abs(long_full$bias) - ts*long_full$whole_SE) <= 0
long_full$bias_tau2<-long_full$whole_Tau2 - long_full$tau2
long_full$bias_tau2_lnR<-log(long_full$whole_Tau2 / long_full$tau2)

# Add a unique code for each parameter set
long_full$code<-paste0(long_full[,4], "_", long_full[,5])
for(i in 6:15){
	long_full$code<-paste0(long_full$code, "_", long_full[,i])
	print(i)
}

# Save that
save(long_full, file="HPC_data/long_REMA.Rdata")

# Now get the aggregate stats
long_full$code2<-paste0(long_full$code, long_full$Method)
agg_stats<-ddply(long_full, .(code2), summarise, code=code[1], Method=Method[1], mean_bias=mean(bias), median_bias=median(bias), min_bias=min(bias), max_bias=max(bias), q10_bias=quantile(bias, 0.1), q90_bias=quantile(bias, 0.9), mean_bias_tau2=mean(bias_tau2), median_bias_tau2=median(bias_tau2), min_bias_tau2=min(bias_tau2), max_bias_tau2=max(bias_tau2), q10_bias_tau2=quantile(bias_tau2, 0.1), q90_bias_tau2=quantile(bias_tau2, 0.9), mean_bias_tau2_lnR=mean(bias_tau2_lnR), median_bias_tau2_lnR=median(bias_tau2_lnR), min_bias_tau2_lnR=min(bias_tau2_lnR), max_bias_tau2_lnR=max(bias_tau2_lnR), q10_bias_tau2_lnR=quantile(bias_tau2_lnR, 0.1), q90_bias_tau2_lnR=quantile(bias_tau2_lnR, 0.9), coverage=mean(coverage))

# Combine with the parameter set by code
parameters$code<-paste0(parameters[,1], "_", parameters[,2])
for(i in 3:12){
	parameters$code<-paste0(parameters$code, "_", parameters[,i])
	print(i)
}
agg_stats<-cbind(agg_stats, parameters[match(agg_stats$code, parameters$code),])
agg_results<-agg_stats[,-c(1:2)]

# Save that
save(agg_results, file="Miss_SD_Sim/agg_results_rema.Rdata")


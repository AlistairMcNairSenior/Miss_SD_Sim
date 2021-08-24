
# Simulation to Support Nakagawa et al. A practical and readily implementable method for handling missing standard deviation in the meta-analysis of response ratios. 

# This script includes a random effect for study ID (i.e. multilevel meta-analysis)

# First written by AM Senior @ The University of Sydney, 25/06/2021.

# Clean up the R Environment 
rm(list=ls())

# Set the working directory
wd<-"/Users/alistairsenior/Dropbox (Sydney Uni)/Nakagawa_Ecology_Missing_SD/Full_Simulation/Miss_Sim"
setwd(wd)

# Load the relevant libraries, and header
library(metafor)
library(plyr)
library(doSNOW)

# Check the files
files<-dir()[-c(1:7)]

load(files[1])
res_unlist<-results
for(i in 1:length(files)){
	load(files[i])
	for(j in 1:length(res_unlist)){
		res_unlist[[j]]<-rbind(res_unlist[[j]], results[[j]])
	}
}
results<-res_unlist

# Load the parameters 
load("Parameters_MLMA.Rdata")

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
long_full<-as.data.frame(rbind(long_mat[,-c(5:20)], long_mat[,-c(1:4,9:20)], long_mat[,-c(1:8,13:20)], long_mat[,-c(1:12,17:20)], long_mat[,-c(1:16)]))
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
long_full$bias_ICC<-long_full$whole_ICC - long_full$icc_study

# Add a unique code for each parameter set
long_full$code<-paste0(long_full[,5], "_", long_full[,6])
for(i in 7:16){
	long_full$code<-paste0(long_full$code, "_", long_full[,i])
	print(i)
}

# Save that
save(long_full, file="long_MLMA.Rdata")

# Now get the aggregate stats
long_full$code2<-paste0(long_full$code, long_full$Method)
agg_stats<-ddply(long_full, .(code2), summarise, code=code[1], Method=Method[1], mean_bias=mean(bias), median_bias=median(bias), min_bias=min(bias), max_bias=max(bias), q10_bias=quantile(bias, 0.1), q90_bias=quantile(bias, 0.9), mean_bias_tau2=mean(bias_tau2), median_bias_tau2=median(bias_tau2), min_bias_tau2=min(bias_tau2), max_bias_tau2=max(bias_tau2), q10_bias_tau2=quantile(bias_tau2, 0.1), q90_bias_tau2=quantile(bias_tau2, 0.9), mean_bias_tau2_lnR=mean(bias_tau2_lnR), median_bias_tau2_lnR=median(bias_tau2_lnR), min_bias_tau2_lnR=min(bias_tau2_lnR), max_bias_tau2_lnR=max(bias_tau2_lnR), q10_bias_tau2_lnR=quantile(bias_tau2_lnR, 0.1), q90_bias_tau2_lnR=quantile(bias_tau2_lnR, 0.9),  mean_bias_ICC=mean(bias_ICC), median_ICC=median(bias_ICC), min_bias_ICC=min(bias_ICC), max_ICC=max(bias_ICC), q10_bias_ICC=quantile(bias_ICC, 0.1), q90_bias_ICC=quantile(bias_ICC, 0.9), coverage=mean(coverage))

# Combine with the parameter set by code
parameters$code<-paste0(parameters[,1], "_", parameters[,2])
for(i in 3:12){
	parameters$code<-paste0(parameters$code, "_", parameters[,i])
	print(i)
}
agg_stats<-cbind(agg_stats, parameters[match(agg_stats$code, parameters$code),])
agg_results<-agg_stats[,-c(1:2)]

# Save that
save(agg_results, file="agg_results.Rdata")


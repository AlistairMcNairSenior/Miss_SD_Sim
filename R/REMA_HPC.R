
# Simulation to Support Nakagawa et al. A practical and readily implementable method for handling missing standard deviation in the meta-analysis of response ratios. 

# This script runs the simulation for the MLMA (i.e. multilevel meta-analysis). It is designed to run as an array job on the HPC and Sydney UNI. There are only 10 reps in here, but I ran on an array 0-999

# Clean up the R Environment 
rm(list=ls())

# Incoming arguments from bash
args<-commandArgs()

# Where are we working
directory<-"/project/RDS-FSC-EvolNutStrats-RW/Miss_Sim"
setwd(directory)

# Load the relevant libraries, and header
library(metafor)
library(plyr)
source("Header.R")

# Get the incoming info on the PBS array index
index<-args[6]

# In this simulation I will include non-independence. The the analysis will be MLMA.

###################################################
################## Parameters #####################
###################################################

# Specify the overall effect magnitude - lets assume a lnRR of 0.3 - this is effectively a 35% increase in the mean
lnRR<-0.3

# Specify the number of effect sizes per study 1 with SD 0 for REMA
k_effect_mu<-3
k_effect_sd<-2.4

# Specify total tau2 of lnRR^2. This sets 0 at 1 SD below the mean effect meaning there will be a set proportion of negative effects. A large degree of heteorgeneity
tau2<-c(0.003^2, lnRR^2)

# Specify the % tau2 100 - gives REMA data
icc_study<-c(1)

# Note I assume the mean in the control group is 100 - to give +ve means for lnRR
mu_control<-100

# Specify the mean sample size, and its variance - to be simulated from a double poisson distribution
n_study_mu<-c(5, 30)
n_study_sd<-c(2.7, 6.7)
# Note we dont want all combos - I'll edit out the low SD with low mean and vice versa

# Specify the mean study SD in outcome, I will set to 15% of mu_control; i.e. the CV is 0.15
sd_study_mu<-mu_control * 0.15

# But I think we need to test a few degrees of heterogeneity in the variation
sd_study_sd<-sd_study_mu * seq(0, 0.5, 0.25)

# We will assume a small, medium and large meta-analysis
k_study<-c(12, 30, 100)

# We will drop the SD from, numbers_lost proportion of the data
proportion_lost<-seq(0.05, 0.55, 0.1)

# Number of reps for this array of the simulation - it is run on an array of 1000 cores, therefore 10 reps per core
n_reps<-10

# Get all the parameter combinations to run
parameters<-expand.grid(lnRR=lnRR, k_effect_mu=k_effect_mu, k_effect_sd=k_effect_sd, tau2=tau2, icc_study=icc_study, mu_control=mu_control, n_study_mu=n_study_mu, n_study_sd=n_study_sd, sd_study_mu=sd_study_mu, sd_study_sd=sd_study_sd, k_study=k_study, proportion_lost=proportion_lost)

# Remove the parameters with the combinations of low sd n and high mean sd n etc
drop<-c(which(parameters$n_study_mu == 5 & parameters$n_study_sd == 6.7), which(parameters$n_study_mu == 30 & parameters$n_study_sd == 2.7))
parameters<-parameters[-drop,]
rownames(parameters)<-seq(nrow(parameters))

# Create a template to hold the simulation results. For each of the replications, we will perform the analysis by the three methods and record the est, the SE and the tau, then compare to above
template<-data.frame(whole_ests=rep(NA, n_reps), whole_SE=NA, whole_Tau2=NA, m1.1_ests=NA, m1.1_SE=NA, m1.1_Tau2=NA, m1.2_ests=NA, m1.2_SE=NA, m1.2_Tau2=NA, m2_ests=NA, m2_SE=NA, m2_Tau2=NA, m3_ests=NA, m3_SE=NA, m3_Tau2=NA)

# Output list to fill in there will be an object for each set of parameters in the batch - that object will be the filled in template
results<-list()

# Loop for the each parameter set
for(p in 1:nrow(parameters)){

	# Get the pth set of parameters
	parameters_p<-parameters[p,]
				
	# Create a copy of the template to fill in for this set of parameters
	template_p<-template
	
	# We will do it reps times, each time deleting n numbers_lost of the SDs and reanalysing by the 3 methods
	for(n in 1:n_reps){
		
		# Keep doing until we get the whole nth row of template_p filled in - repeats simulation in the event that any model has convergence issues
		while(sum(is.na(template_p[n,])) > 0){
		
			# Create an nth data set
			data_n<-sim_data(lnRR=parameters_p$lnRR, k_study=parameters_p$k_study, k_effect_mu=parameters_p$k_effect_mu, k_effect_sd=parameters_p$k_effect_sd, tau2=parameters_p$tau2, icc_study=parameters_p$icc_study, n_study_mu=parameters_p$n_study_mu, n_study_sd=parameters_p$n_study_sd, sd_study_mu=parameters_p$sd_study_mu, sd_study_sd=parameters_p$sd_study_sd)
			
			# Calculate the CVs and effect sizes for the complete data and analyse
			data_n$Control.CV<-data_n$Control.SD / data_n$Control.Mean
			data_n$Treatment.CV<-data_n$Treatment.SD / data_n$Treatment.Mean
			data_n<-cbind(data_n, my_lnRR(Control.Mean=data_n$Control.Mean, Treatment.Mean=data_n$Treatment.Mean, Control.CV=data_n$Control.CV, Treatment.CV=data_n$Treatment.CV, Control.n=data_n$Control.n, Treatment.n=data_n$Treatment.n))
			
			# Fit the model to the whole data
			model_w<-try(rma.mv(yi = yi, V = vi, random=list(~1|Effect), data=data_n), silent=T)
			
			# Save the results if the model fitted
			if(class(model_w)[1] == "rma.mv"){
				template_p$whole_ests[n]<-model_w$b
				template_p$whole_SE[n]<-model_w$se
				template_p$whole_Tau2[n]<-model_w$sigma2[1]
			}

			# Find the row IDs you want to have missing SDs and set SD, CV and effect sizes as NA
			drop<-sample(seq(1, nrow(data_n), 1), round(nrow(data_n) * parameters_p$proportion_lost))
			data_n$Control.SD[drop]<-NA
			data_n$Treatment.SD[drop]<-NA
			data_n$Control.CV[drop]<-NA
			data_n$Treatment.CV[drop]<-NA
			data_n$yi[drop]<-NA
			data_n$vi[drop]<-NA
								
			# Here we calculate the weighted mean CV^2 from those data available - this is pooled by study first
			pooled<-ddply(data_n, .(Study), summarise, mean_Control.CV=mean(Control.CV^2, na.rm=T), mean_Treatment.CV=mean(Treatment.CV^2, na.rm=T), mean_Control.n=mean(Control.n), mean_Treatment.n=mean(Treatment.n)) 			
			mcv21<-weighted.mean(pooled$mean_Treatment.CV, pooled$mean_Treatment.n, na.rm=TRUE)
			mcv22<-weighted.mean(pooled$mean_Control.CV, pooled$mean_Control.n, na.rm=TRUE)
										
			# Method 1.1
			# Use the weighted mean cv where SDs are missing.
			# Then proceed to fit a meta-analysis as normal.
			
			# Now estimate a new yi and vi using the average CV^2s where data are missing
			data_n$yi[drop]<-log(data_n$Treatment.Mean[drop] / data_n$Control.Mean[drop]) + 0.5 * (mcv21 / data_n$Treatment.n[drop] - mcv22 / data_n$Control.n[drop])
			data_n$vi[drop]<-mcv21/data_n$Treatment.n[drop] + mcv22/data_n$Control.n[drop] + mcv21^2/(2 * data_n$Treatment.n[drop]^2) + mcv22^2/(2 * data_n$Control.n[drop]^2)
			
			# Fit the model
			model1.1<-try(rma.mv(yi = yi, V = vi, random=list(~1|Effect), data=data_n), silent=T)
			
			# Save the results if the model fitted
			if(class(model1.1)[1] == "rma.mv"){
				template_p$m1.1_ests[n]<-model1.1$b
				template_p$m1.1_SE[n]<-model1.1$se
				template_p$m1.1_Tau2[n]<-model1.1$sigma2[1]
			}
			
			# Method 1.2
			# Use the weighted mean cv everywhere
			# Then proceed to fit a meta-analysis as normal.
			
			# Now estimate a new yi and vi using the average CV^2s everywhere
			data_n$yi_2<-log(data_n$Treatment.Mean / data_n$Control.Mean) + 0.5 * (mcv21 / data_n$Treatment.n - mcv22 / data_n$Control.n)			
			data_n$vi_2<-mcv21/data_n$Treatment.n + mcv22/data_n$Control.n + mcv21^2/(2 * data_n$Treatment.n^2) + mcv22^2/(2 * data_n$Control.n^2)
			
			# Fit the model
			model1.2<-try(rma.mv(yi = yi_2, V = vi_2, random=list(~1|Effect), data=data_n), silent=T)
			
			# Save the results if the model fitted
			if(class(model1.2)[1] == "rma.mv"){
				template_p$m1.2_ests[n]<-model1.2$b
				template_p$m1.2_SE[n]<-model1.2$se
				template_p$m1.2_Tau2[n]<-model1.2$sigma2[1]
			}
			
			# Method 2
			# Now we fit a weighted regression(ish) rather than a meta-analysis in metafor
				
			# Create a matrix with v_tilda in the diagonal (note this is NOT the mix v_tilda and v_i)
			Vf<-diag(mcv21/data_n$Treatment.n + mcv22/data_n$Control.n + mcv21^2/(2 * data_n$Treatment.n^2) + mcv22^2/(2 * data_n$Control.n^2))
			row.names(Vf)<-data_n$Effect
			colnames(Vf)<-data_n$Effect
			
			# Also apparently need to copy the IDs across to second predictor
			data_n$ID2<-data_n$Effect
				
			# Now fit a weighted regression
			model2<-try(rma.mv(yi = yi, V = 0, random=list(~1|Effect, ~1|ID2), data=data_n, R=list(ID2=Vf), Rscale=F), silent=T)
			
			# Save the results if the model fitted
			if(class(model2)[1] == "rma.mv"){
				template_p$m2_ests[n]<-model2$b
				template_p$m2_SE[n]<-model2$se
				template_p$m2_Tau2[n]<-model2$sigma2[1]
			}
			
			# Method 3
			# OK now to combine the methods - we set vi to 0 for missing SDs, create a new V matrix of 0s where data are known and the vi from the mean CV^2 for missings, then fit together in one model
			
			# Switch the inverse n in the matrix with the vi - note this is v.tilda where SDs are missing and vi where they are not, as per method 1
			diag(Vf)<-data_n$vi
			
			# Add 0s to V matrix where SDs are known
			diag(Vf)[-drop]<-0
			
			# Now add 0s back to the vi where SDs are missing to remove v.tilda from where SD is missing
			data_n$vi[drop]<-0
			
			# Try fit the model
			model3<-try(rma.mv(yi = yi, V = vi, random=list(~1|Effect, ~1|ID2), data=data_n, R=list(ID2=Vf), Rscale=F), silent=T)
			
			# Save the results if the model fitted
			if(class(model3)[1] == "rma.mv"){
				template_p$m3_ests[n]<-model3$b
				template_p$m3_SE[n]<-model3$se
				template_p$m3_Tau2[n]<-model3$sigma2[1]
			}
			
		}# Closing the while loop
			
	}# Closing the reps
	
	# Save the results from this n_reps set of parameters
	results[[p]]<-template_p

}# Closing the parameters loop

# Save the results from the simulation
save(results, file=paste0("Simulation_REMA_", index, ".Rdata"))
save(parameters, file="Parameters_REMA.Rdata")




# Simulation to Support Nakagawa et al. A practical and readily implementable method for handling missing standard deviation in the meta-analysis of response ratios. 

# This script includes a random effect for study ID (i.e. multilevel meta-analysis)

# First written by AM Senior @ The University of Sydney, 25/06/2021.

############################################
############### sim_data ###################
############################################

# The function is designed to simulate data at a multi-level with study and within-study
# The same function will all for a more standard design by setting k_study as 1, and k_effects the desired nu,ber of studies, and icc_study = 0, giving all the heterogeneity at the effect size level

# lnRR = the overall mean effect size
# k_study = the number of studies
# k_effect_mu = the mean number of effect sizes per study
# k_effect_sd = the sd in number of effect sizes per study
# tau2 = the total heterogeneity
# icc_study = the proportion (intraclass correlation) of tau2 coming from the study level effect 
# n_study_mu = the mean sample size studies in the dataset
# n_study_sd = the sd of sample size of studies in the dataset
# sd_study_mu is the mean sd of data within effect size for each study
# sd_study_sd is the sd of sd of data within effect size for each study
# return_true_effects allows the user to see the simulated effect at the level of each study and effect size in the outputted data

# NOTE SDs and ns vary at the level of study rather than effect size (or treatment group). In my experiuence this is the most realistic scenario

sim_data<-function(lnRR, k_study, k_effect_mu, k_effect_sd, tau2, icc_study, n_study_mu, n_study_sd, sd_study_mu, sd_study_sd, return_true_effects=F){
	
	# Load the packages
	require(gamlss.dist)
	
	# A bit of reparameterization
	
	# Make sure the various sds are not <= 0: for the gamma and DPO distributions
	if(sd_study_sd <= 0){
		sd_study_sd<-10^-8
	}
	if(k_effect_sd <= 0){
		k_effect_sd<-10^-8
	}
	if(n_study_sd <= 0){
		n_study_sd<-10^-8
	}
	
	# Calculate the between and within-study heterogeneity from the tau2 and icc study
	study_sigma2<-tau2 * icc_study
	resid_sigma2<-tau2 * (1 - icc_study)
	
	# For the gamma distirbution for the variance the a (shape) and b (scale) are approximated as
	a<-sd_study_mu^2 / sd_study_sd^2
	b<-sd_study_sd^2 / sd_study_mu
	
	# For the double poisson distirbution the sigma is parameterised as - I will specify the mean as n_study_mu - 3, then add 3 to avoid 0 sample sizes - same for k_effects. Note the addition of 0.0001 avoids 0 means
	sigma_n<-n_study_sd^2 / (n_study_mu - 3 + 10^-8)
	sigma_k<-k_effect_sd^2 / (k_effect_mu - 1 + 10^-8)
 	
 	# Loop for the k studies
	for(i in 1:k_study){
		
		# Get the study specific-effect for study i
		ES_i<-rnorm(1, lnRR, sqrt(study_sigma2))
		
		# Find study-specific number of effect sizes, sample size and variance for study i
		k_effect_i<-rDPO(1, mu=(k_effect_mu - 1 + 10^-8), sigma=sigma_k) + 1
		n_i<-rDPO(1, mu=(n_study_mu - 3 + 10^-8), sigma=sigma_n) + 3
		sd_i<-rgamma(1, shape=a, scale=b)
		
		# Loop for the k effects
		for(j in 1:k_effect_i){
			
			# Get the within-study effects for study i effect size j
			ES_ij<-rnorm(1, ES_i, sqrt(resid_sigma2))	
		
			# Simulate the control group for study i, effect size j
			control_ij<-rnorm(n=n_i, mean=mu_control, sd=sd_i)
			
			# Simulate the treatment group for study i, effect size j
			treatment_ij<-rnorm(n_i, mean(mu_control * exp(ES_ij)), sd=sd_i)
		
			# Now calculate the mean, sd and give n in a nice way for downstream code
			data_ij<-data.frame(Effect = paste(i, j, sep="_"), Study = i, Control.Mean = mean(control_ij), Control.SD = sd(control_ij), Control.n = length(control_ij), Treatment.Mean = mean(treatment_ij), Treatment.SD = sd(treatment_ij), Treatment.n = length(treatment_ij))
			
			# If you want the true simulated effects returned, add those in
			if(return_true_effects){
				data_ij$ES_i<-ES_i
				data_ij$ES_ij<-ES_ij
				data_ij$sd_i<-sd_i
			}
			
			# Keep a dataframe of all the data
			if(i == 1 & j == 1){
				data<-data_ij
			}else{
				data<-rbind(data, data_ij)	
			}	
		}
	}
	
	# Return the data
	return(data)
	
}

############################################
################ my_lnRR ###################
############################################

# Function to calcualte the lnRR based on CV as given in the manuscript

my_lnRR<-function(Control.Mean, Treatment.Mean, Control.CV, Treatment.CV, Control.n, Treatment.n){
	
	# The Effect size
	yi<-log(Treatment.Mean / Control.Mean) + 0.5 * (Treatment.CV^2 / Treatment.n - Control.CV^2 / Control.n)

	# The sampling variance
	vi<-(Treatment.CV^2 / Treatment.n) + (Control.CV^2 / Control.n) + (Treatment.CV^4 / (2 * Treatment.n^2)) + (Control.CV^4 / (2 * Control.n^2))
	
	# Return the data
	return(data.frame(yi=yi, vi=vi))
	
}


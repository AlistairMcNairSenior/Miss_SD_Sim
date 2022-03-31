
# Simulation to Support Nakagawa et al. A practical and readily implementable method for handling missing standard deviation in the meta-analysis of response ratios. 

# This script includes a random effect for study ID (i.e. multilevel meta-analysis)

# First written by AM Senior @ The University of Sydney, 25/06/2021.

# Clean up the R Environment 
rm(list=ls())

# Set the working directory
wd<-"/Users/alistairsenior/OneDrive - The University of Sydney (Staff)/Nakagawa_Ecology_Missing_SD/Miss_SD_Sim" # MacMini Home and work iMac
setwd(wd)

# Load the relevant libraries, and header
library(ggplot2)
library(plyr)
library(gridExtra)
library(ggbeeswarm)
library(ggcorrplot)

# Load the results
load("Rdata/agg_results_rema.Rdata")

########################################################################
############## Averaging each parameter set ############################
########################################################################

# Analysis of Bias in the lnRR

head(agg_results)
agg_results$range_bias<-agg_results$max_bias - agg_results$min_bias

# Rename 1.1 and 1.2 to 1.A and 1.B
agg_results$Method[which(agg_results$Method == "Method 1.1")]<-"Method 1A"
agg_results$Method[which(agg_results$Method == "Method 1.2")]<-"Method 1B"
agg_results$Method[which(agg_results$Method == "Complete Data")]<-"Full Data"

A<-ggplot(data=agg_results, aes(x=median_bias, y=Method, col=Method, fill=Method)) +
	geom_violin() + 
	geom_beeswarm(groupOnX=F, size=0.05, cex=0.25, col="black") + 
	xlim(-0.005, 0.008) +
	xlab("Bias lnRR") + ylab("") + theme_bw() + 
	theme(legend.position="none", axis.title.x=element_text(size=15), axis.text.y=element_text(size=15), plot.title=element_text(size=15)) +
	geom_vline(xintercept=0, col="dark grey", size=0.5)


long<-agg_results[,c("Method", "median_bias", "code")]
wide<-reshape(long, direction="wide", idvar="code", timevar="Method")
cor.mat<-cor(wide[,-1])
rownames(cor.mat)<-c("Full", "Method 1A", "Method 1B", "Method 2", "Method 3")
colnames(cor.mat)<-c("Full", "Method 1A", "Method 1B", "Method 2", "Method 3")

B<-ggcorrplot(cor.mat, lab=T, show.legend=FALSE, type="lower") + theme(plot.title=element_text(size=15))

# Difference in the bias
df<-data.frame(delta=abs(wide[,"median_bias.Method 1A"]) - abs(wide[,"median_bias.Method 1B"]))
C<-ggplot(df, aes(x=delta)) + geom_histogram(bins=20) + theme_bw() + geom_vline(xintercept=0, col="red") +
	xlab("Diff. |Bias| Meths 1A & 1B") + ylab("Frequency")

D<-ggplot(data=agg_results, aes(x=log(range_bias, base=10), y=Method, col=Method, fill=Method)) +
	geom_violin() + 
	geom_beeswarm(groupOnX=FALSE, size=0.05, cex=0.25, col="black") + 
	xlab("log10 Range of Bias lnRR") + ylab("") + theme_bw() + 
	theme(legend.position="none", axis.title.x=element_text(size=15), axis.text.y=element_text(size=15), plot.title=element_text(size=15))

# Get the data for Methods 1.A & 1.B
dat1.1<-agg_results[which(agg_results$Method=="Method 1A"),]
dat1.2<-agg_results[which(agg_results$Method=="Method 1B"),]

E<-ggplot(dat1.2, aes(x=log(range_bias, base=10), y=as.factor(sd_study_sd), fill=as.factor(n_study_mu))) +
	geom_violin() +
	geom_point(size=-1, col="grey") +
	theme_bw() +
	xlab("log10 Range of Bias lnRR") +
	ylab(expression(sigma[s])) + 
	theme(axis.title.x=element_text(size=15), axis.text.y=element_text(size=15), axis.title.y=element_text(size=15), legend.position=c(0.8, 0.2), legend.text=element_text(size=15), legend.title=element_text(size=15), plot.title=element_text(size=15)) + guides(fill=guide_legend(title="Study Sample Size")) + xlim(-1.75, 1.25)	



pdf("MS/fig/Bias_lnRR_REMA.pdf", height=10, width=12)

grid.arrange(A+labs(title="A."), B+labs(title="B."), C+labs(title="C."), D+labs(title="D."), E+labs(title="E."), layout_matrix=rbind(c(1,1,3,4),
													c(5,5,6,6)))	
													
dev.off()

# Analysis of bias in coverage

A<-ggplot(data=agg_results, aes(x=coverage, y=Method, col=Method, fill=Method)) +
	geom_violin() + 
	geom_beeswarm(groupOnX=F, size=0.05, cex=0.25, col="black") + 
	xlim(0.85, 1) +
	xlab("Coverage (95% CI)") + ylab("") + theme_bw() + 
	theme(legend.position="none", axis.title.x=element_text(size=15), axis.text.y=element_text(size=15), plot.title=element_text(size=15)) +
	geom_vline(xintercept=0.95, col="dark grey", size=0.5)
	
B<-ggplot(data=dat1.1, aes(x=as.factor(tau2), y=coverage)) +
	geom_violin(fill="lightgoldenrod3") + 
	theme_bw() + 
	ylim(0.9, 1) +
	xlab(expression(italic(T)^2)) + ylab("Coverage (95% CI)") +
	theme(axis.text.x=element_text(size=12), axis.title.x=element_text(size=15), axis.text.y=element_text(size=12), axis.title.y=element_text(size=15), plot.title=element_text(size=15), plot.subtitle=element_text(size=15)) +
	geom_hline(yintercept=0.95, col="dark grey") +
	scale_x_discrete(labels=c(expression(9e-06~~or~~over(italic(T), lnRR)==0.01), expression(0.09~~or~~over(italic(T), lnRR)==1))) + labs(subtitle="Method 1A")

C<-ggplot(data=dat1.2, aes(x=as.factor(tau2), y=coverage)) +
	geom_violin(fill="aquamarine3") + 
	theme_bw() + 
	ylim(0.9, 1) +
	xlab(expression(italic(T)^2)) + ylab("Coverage (95% CI)") +
	theme(axis.text.x=element_text(size=12), axis.title.x=element_text(size=15), axis.text.y=element_text(size=12), axis.title.y=element_text(size=15), plot.title=element_text(size=15), plot.subtitle=element_text(size=15)) +
	geom_hline(yintercept=0.95, col="dark grey") +
	scale_x_discrete(labels=c(expression(9e-06~~or~~over(italic(T), lnRR)==0.01), expression(0.09~~or~~over(italic(T), lnRR)==1))) + labs(subtitle="Method 1B")

pdf("MS/fig/Coverage_REMA.pdf", height=7, width=18)

grid.arrange(A+labs(title="A."), B+labs(title="B."), C+labs(title="C."), layout_matrix=array(c(1,2,3), c(1,3)))
													
dev.off()
	
# Analysis of bias in Tau2
	
A<-ggplot(data=agg_results, aes(x=median_bias_tau2_lnR, y=Method, col=Method, fill=Method)) +
	geom_violin() + 
	geom_beeswarm(groupOnX=F, size=0.05, cex=0.25, col="black") + 
	xlim(-15, 6) +
	xlab("Bias in Heterogeneity (log Ratio)") + ylab("") + theme_bw() + 
	theme(legend.position="none", axis.title.x=element_text(size=15), axis.text.x=element_text(size=12), axis.text.y=element_text(size=15), plot.title=element_text(size=15)) +
	geom_vline(xintercept=0, col="dark grey", size=0.5)

B<-ggplot(data=agg_results, aes(x=as.factor(tau2), y= median_bias_tau2_lnR, fill=Method)) +
	geom_boxplot() + theme_bw() + 
	xlab(expression(italic(T)^2)) + ylab("Bias in Heterogeneity (log Ratio)") +
	theme(legend.position="none", axis.text.x=element_text(size=12), axis.title.x=element_text(size=15), axis.text.y=element_text(size=12), axis.title.y=element_text(size=15), plot.title=element_text(size=15)) +
	geom_hline(yintercept=0, col="dark grey") +
	scale_x_discrete(labels=c(expression(9e-06~~or~~over(italic(T), lnRR)==0.01), expression(0.09~~or~~over(italic(T), lnRR)==1)))
	

pdf("MS/fig/Bias_Het_REMA.pdf", height=5, width=12)

grid.arrange(A+labs(title="A."), B+labs(title="B."), layout_matrix=rbind(c(1,2)))	
													
dev.off()	
	

# Simulation to Support Nakagawa et al. A practical and readily implementable method for handling missing standard deviation in the meta-analysis of response ratios. 

# This script includes a random effect for study ID (i.e. multilevel meta-analysis)

# First written by AM Senior @ The University of Sydney, 25/06/2021.

# Clean up the R Environment 
rm(list=ls())

# Set the working directory
wd<-"/Users/alistairsenior/OneDrive - The University of Sydney (Staff)/Nakagawa_Ecology_Missing_SD/Miss_SD_Sim" # MacMini Home
setwd(wd)

# Load the relevant libraries, and header
library(ggplot2)
library(plyr)
library(gridExtra)
library(ggbeeswarm)
library(ggcorrplot)
library(ggthemes)

# Load the results
load("Rdata/agg_results_mlma.Rdata")

# Params
title_size<-20

########################################################################
############## Averaging each parameter set ############################
########################################################################

# Analysis of Bias in the lnRR

head(agg_results)
agg_results$range_bias<-agg_results$max_bias - agg_results$min_bias

# Rename 1.1 and 1.2 to 1.A and 1.B, and complete to full
# agg_results$Method_old<-agg_results$Method
# agg_results$Method[which(agg_results$Method == "Method 1.1")]<-"Missing Cases"
# agg_results$Method[which(agg_results$Method == "Method 1.2")]<-"All Cases"
# agg_results$Method[which(agg_results$Method == "Method 2")]<-"Multiplicative"
# agg_results$Method[which(agg_results$Method == "Method 3")]<-"Hybrid"
# agg_results$Method[which(agg_results$Method == "Complete Data")]<-"Full Data"
method_cols<-colorblind_pal()(8)[c(2,3,4,7,8)]
other_cols<-colorblind_pal()(8)[c(1,5)]

A<-ggplot(data=agg_results, aes(y=median_bias, x=Method, col=Method, fill=Method)) + 
	geom_violin() + 
	geom_beeswarm(groupOnX=T, size=0.6, cex=0.2, col="black", alpha=1/3, shape=1) + 
	ylim(-0.0025, 0.008) +
	ylab("Bias lnRR") + xlab("") + theme_bw() + 
	theme(legend.position="none", axis.title.y=element_text(size=15), axis.text.y=element_text(size=15), axis.text.x=element_text(size=15), plot.title=element_text(size=title_size)) +
	geom_vline(xintercept=0, col="dark grey", size=0.5) + 
	scale_x_discrete(labels=c("Full Data", "Miss.
Cases", "All
Cases", "Multiplic.", "Hybrid")) + scale_color_manual(values=method_cols) + scale_fill_manual(values=method_cols)


long<-agg_results[,c("Method", "median_bias", "code")]
wide<-reshape(long, direction="wide", idvar="code", timevar="Method")
cor.mat<-cor(wide[,-1])
rownames(cor.mat)<-c("Full Data", "Miss. Cases", "All Cases", "Multiplicative", "Hybrid")
colnames(cor.mat)<-c("Full Data", "Miss. Cases", "All Cases", "Multiplicative", "Hybrid")

B<-ggcorrplot(cor.mat, lab=T, show.legend=FALSE, type="lower") + theme(plot.title=element_text(size=title_size), axis.text.x=element_text(size=15),  axis.text.y=element_text(size=15))

# Difference in the bias
df<-data.frame(delta=abs(wide[,"median_bias.Missing Cases"]) - abs(wide[,"median_bias.All Cases"]))
C<-ggplot(df, aes(x=delta)) + geom_histogram(bins=20) + theme_bw() + geom_vline(xintercept=0, col="red") +
	xlab("Diff. |Bias| Miss. & All Cases") + ylab("Frequency") + 
	theme(axis.title.y=element_text(size=15), axis.text.y=element_text(size=15), axis.text.x=element_text(size=12), axis.title.x=element_text(size=12), plot.title=element_text(size=title_size))

D<-ggplot(data=agg_results, aes(y=log(range_bias, base=10), x=Method, col=Method, fill=Method)) +
	geom_violin() + 
	geom_beeswarm(groupOnX=T, size=2, cex=0.5, col="black", alpha=1/3, shape=1) + 
	ylab("log10 Range of Bias lnRR") + xlab("") + theme_bw() + 
	theme(legend.position="none", axis.title.y=element_text(size=15), axis.text.y=element_text(size=15), axis.text.x=element_text(size=15), plot.title=element_text(size=title_size)) + 
	scale_x_discrete(labels=c("Full Data", "Miss. 
Cases", "All 
Cases", "Multiplic.", "Hybrid")) + scale_color_manual(values=method_cols) + scale_fill_manual(values=method_cols)


# Get the data for Methods 1.A & 1.B
dat1.1<-agg_results[which(agg_results$Method=="Method 1.1"),]
dat1.2<-agg_results[which(agg_results$Method=="Method 1.2"),]

E<-ggplot(dat1.2, aes(y=log(range_bias, base=10), x=as.factor(sd_study_sd), fill=as.factor(n_study_mu))) +
	geom_violin() +
	geom_point(size=-1, col="grey") +
	theme_bw() +
	ylab("log10 Range of Bias lnRR") +
	xlab(expression(sigma[s])) + 
	theme(axis.title.x=element_text(size=20), axis.text.y=element_text(size=15), axis.title.y=element_text(size=15), axis.text.x=element_text(size=15), legend.position=c(0.2, 0.8), legend.text=element_text(size=15), legend.title=element_text(size=15), plot.title=element_text(size=title_size)) + guides(fill=guide_legend(title="Study Sample Size")) + scale_color_manual(values=other_cols) + scale_fill_manual(values=other_cols)	

pdf("MS/fig/Bias_lnRR.pdf", height=10, width=12)

grid.arrange(A+labs(title="A."), B+labs(title="B."), C+labs(title="C."), D+labs(title="D."), E+labs(title="E."), layout_matrix=rbind(c(1,1,3,4),
													c(5,5,6,6)))	
													
dev.off()

# Analysis of bias in coverage

A<-ggplot(data=agg_results, aes(x=coverage, y=Method, col=Method, fill=Method)) +
	geom_violin() + 
	geom_beeswarm(groupOnX=F, size=0.05, cex=0.25, col="black") + 
	xlim(0.9, 1) +
	xlab("Coverage (95% CI)") + ylab("") + theme_bw() + 
	theme(legend.position="none", axis.title.x=element_text(size=15), axis.text.y=element_text(size=15), plot.title=element_text(size=15)) +
	geom_vline(xintercept=0.95, col="dark grey", size=0.5)
	
B<-ggplot(data=dat1.1, aes(x=as.factor(tau2), y=coverage, fill=as.factor(icc_study))) +
	geom_violin() + 
	theme_bw() + 
	ylim(0.9, 1) +
	xlab(expression(italic(T)^2)) + ylab("Coverage (95% CI)") +
	theme(legend.position=c(0.85, 0.875), axis.text.x=element_text(size=12), axis.title.x=element_text(size=15), axis.text.y=element_text(size=12), axis.title.y=element_text(size=15), legend.title=element_text(size=15), plot.title=element_text(size=15), plot.subtitle=element_text(size=15), legend.text=element_text(size=15)) +
	geom_hline(yintercept=0.95, col="dark grey") + guides(fill=guide_legend(title="Simulated ICC")) +
	scale_x_discrete(labels=c(expression(9e-06~~or~~over(italic(T), lnRR)==0.01), expression(0.09~~or~~over(italic(T), lnRR)==1))) + labs(subtitle="Method 1A")

C<-ggplot(data=dat1.2, aes(x=as.factor(tau2), y=coverage, fill=as.factor(icc_study))) +
	geom_violin() + 
	theme_bw() + 
	xlab(expression(italic(T)^2)) + ylab("Coverage (95% CI)") +
	ylim(0.9, 1) +
	theme(legend.position=c(0.85, 0.875), axis.text.x=element_text(size=12), axis.title.x=element_text(size=15), axis.text.y=element_text(size=12), axis.title.y=element_text(size=15), legend.title=element_text(size=15), plot.title=element_text(size=15), plot.subtitle=element_text(size=15), legend.text=element_text(size=15)) +
	geom_hline(yintercept=0.95, col="dark grey") + guides(fill=guide_legend(title="Simulated ICC")) +
	scale_x_discrete(labels=c(expression(9e-06~~or~~over(italic(T), lnRR)==0.01), expression(0.09~~or~~over(italic(T), lnRR)==1))) + labs(subtitle="Method 1B")

pdf("MS/fig/Coverage.pdf", height=7, width=18)

grid.arrange(A+labs(title="A."), B+labs(title="B."), C+labs(title="C."), layout_matrix=array(c(1,2,3), c(1,3)))
													
dev.off()
	
# Analysis of bias in Tau2
	
A<-ggplot(data=agg_results, aes(x=median_bias_tau2_lnR, y=Method, col=Method, fill=Method)) +
	geom_violin() + 
	geom_beeswarm(groupOnX=F, size=0.05, cex=0.25, col="black") + 
	xlim(-1, 6) +
	xlab("Bias in Heterogeneity (log Ratio)") + ylab("") + theme_bw() + 
	theme(legend.position="none", axis.title.x=element_text(size=15), axis.text.x=element_text(size=12), axis.text.y=element_text(size=15), plot.title=element_text(size=15)) +
	geom_vline(xintercept=0, col="dark grey", size=0.5)

B<-ggplot(data=agg_results, aes(x=as.factor(tau2), y=median_bias_tau2_lnR, fill=Method)) +
	geom_boxplot() + theme_bw() + 
	xlab(expression(italic(T)^2)) + ylab("Bias in Heterogeneity (log Ratio)") +
	theme(legend.position="none", axis.text.x=element_text(size=12), axis.title.x=element_text(size=15), axis.text.y=element_text(size=12), axis.title.y=element_text(size=15), plot.title=element_text(size=15)) +
	geom_hline(yintercept=0, col="dark grey") +
	scale_x_discrete(labels=c(expression(9e-06~~or~~over(italic(T), lnRR)==0.01), expression(0.09~~or~~over(italic(T), lnRR)==1)))

C<-ggplot(data=agg_results, aes(x=median_ICC, y=Method, col=Method, fill=Method)) +
	geom_violin() + 
	geom_beeswarm(groupOnX=F, size=0.05, cex=0.25, col="black") + 
	xlim(-0.5, 0.75) +
	xlab("Bias in ICC") + ylab("") + theme_bw() + 
	theme(legend.position="none", axis.title.x=element_text(size=15), axis.text.x=element_text(size=12), axis.text.y=element_text(size=15), plot.title=element_text(size=15)) +
	geom_vline(xintercept=0, col="dark grey", size=0.5)
		
D<-ggplot(data=dat1.2, aes(x=as.factor(tau2), y=median_ICC, fill=as.factor(icc_study))) +
	geom_violin() + 
	theme_bw() + 
	xlab(expression(italic(T)^2)) + ylab("Bias in ICC") +
	theme(legend.position=c(0.85, 0.85), axis.text.x=element_text(size=12), axis.title.x=element_text(size=15), axis.text.y=element_text(size=12), axis.title.y=element_text(size=15), legend.title=element_text(size=15), plot.title=element_text(size=15), legend.text=element_text(size=15)) +
	geom_hline(yintercept=0, col="dark grey") + guides(fill=guide_legend(title="Simulated ICC")) +
	scale_x_discrete(labels=c(expression(9e-06~~or~~over(italic(T), lnRR)==0.01), expression(0.09~~or~~over(italic(T), lnRR)==1)))	
	

pdf("MS/fig/Bias_Het.pdf", height=10, width=12)

grid.arrange(A+labs(title="A."), B+labs(title="B."), C+labs(title="C."), D+labs(title="D."), layout_matrix=rbind(c(1,2),
					c(3,4)))	
													
dev.off()	
	
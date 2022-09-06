
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

# Note for coding purposes in the Method column
# "Method 1.1" = "Missing Cases"
# "Method 1.2" = "All Cases"
# "Method 2" = "Multiplicative"
# "Method 3" = "Hybrid"
# "Complete Data" = "Full Data"
method_cols<-colorblind_pal()(8)[c(2,3,4,7,8)]
other_cols<-colorblind_pal()(8)[c(1,5)]

A<-ggplot(data=agg_results, aes(y=median_bias, x=Method, col=Method, fill=Method)) + 
	geom_violin() + 
	geom_beeswarm(groupOnX=T, size=0.6, cex=0.2, col="black", alpha=1/3, shape=1) + 
	ylim(-0.0025, 0.008) +
	ylab("Bias lnRR") + xlab("") + theme_bw() + 
	theme(legend.position="none", axis.title.y=element_text(size=15), axis.text.y=element_text(size=15), axis.text.x=element_text(size=15), plot.title=element_text(size=title_size)) +
	geom_hline(yintercept=0, col="dark grey", size=0.5) + 
	scale_x_discrete(labels=c("Full Data", "Missing
Cases", "All
Cases", "Multiplicative", "Hybrid")) + scale_color_manual(values=method_cols) + scale_fill_manual(values=method_cols)


long<-agg_results[,c("Method", "median_bias", "code")]
wide<-reshape(long, direction="wide", idvar="code", timevar="Method")
cor.mat<-cor(wide[,-1])
rownames(cor.mat)<-c("Full Data", "Missing Cases", "All Cases", "Multiplicative", "Hybrid")
colnames(cor.mat)<-c("Full Data", "Missing Cases", "All Cases", "Multiplicative", "Hybrid")

B<-ggcorrplot(cor.mat, lab=T, show.legend=FALSE, type="lower") + theme(plot.title=element_text(size=title_size), axis.text.x=element_text(size=15),  axis.text.y=element_text(size=15))

# Difference in the bias
df<-data.frame(delta=abs(wide[,"median_bias.Method 1.1"]) - abs(wide[,"median_bias.Method 1.2"]))
C<-ggplot(df, aes(x=delta)) + geom_histogram(bins=20) + theme_bw() + geom_vline(xintercept=0, col="red") +
	xlab("Difference in |Bias| 
Missing - All Cases") + ylab("Frequency") + 
	theme(axis.title.y=element_text(size=15), axis.text.y=element_text(size=15), axis.text.x=element_text(size=15), axis.title.x=element_text(size=15), plot.title=element_text(size=title_size))

D<-ggplot(data=agg_results, aes(y=log(range_bias, base=10), x=Method, col=Method, fill=Method)) +
	geom_violin() + 
	geom_beeswarm(groupOnX=T, size=2, cex=0.5, col="black", alpha=1/3, shape=1) + 
	ylab("log10 Range of Bias lnRR") + xlab("") + theme_bw() + 
	theme(legend.position="none", axis.title.y=element_text(size=15), axis.text.y=element_text(size=15), axis.text.x=element_text(size=15), plot.title=element_text(size=title_size)) + 
	scale_x_discrete(labels=c("Full Data", "Missing
Cases", "All 
Cases", "Multiplicative", "Hybrid")) + scale_color_manual(values=method_cols) + scale_fill_manual(values=method_cols)


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

A<-ggplot(data=agg_results, aes(y=coverage, x=Method, col=Method, fill=Method)) +
	geom_violin() + 
	geom_beeswarm(groupOnX=T, size=2, cex=0.5, col="black", alpha=1/3, shape=1) + 
	ylim(0.9, 1) +
	ylab("Coverage (95% CI)") + xlab("") + theme_bw() + 
	theme(legend.position="none", axis.title.y=element_text(size=15), axis.text.x=element_text(size=15), axis.text.y=element_text(size=15), plot.title=element_text(size=title_size)) +
	geom_hline(yintercept=0.95, col="dark grey", size=0.5) +
	scale_x_discrete(labels=c("Full Data", "Missing
Cases", "All 
Cases", "Multiplicative", "Hybrid")) + scale_color_manual(values=method_cols) + scale_fill_manual(values=method_cols)

	
B<-ggplot(data=dat1.1, aes(x=as.factor(tau2), y=coverage, fill=as.factor(icc_study))) +
	geom_violin() + 
	theme_bw() + 
	ylim(0.9, 1) +
	xlab("") + ylab("Coverage (95% CI)") +
	theme(legend.position=c(0.85, 0.875), axis.text.x=element_text(size=15), axis.title.x=element_text(size=15), axis.text.y=element_text(size=15), axis.title.y=element_text(size=15), legend.title=element_text(size=15), plot.title=element_text(size=title_size), plot.subtitle=element_text(size=20), legend.text=element_text(size=15)) +
	geom_hline(yintercept=0.95, col="dark grey") + guides(fill=guide_legend(title="Simulated ICC")) +
	scale_x_discrete(labels=c("Low Heterogeneity", "High Heterogeneity")) + scale_color_manual(values=other_cols) + scale_fill_manual(values=other_cols) +
	annotate("text", 0.9, 1, label="Missing Cases", size=8)

C<-ggplot(data=dat1.2, aes(x=as.factor(tau2), y=coverage, fill=as.factor(icc_study))) +
	geom_violin() + 
	theme_bw() + 
	xlab("") + ylab("Coverage (95% CI)") +
	ylim(0.9, 1) +
	theme(legend.position=c(0.85, 0.875), axis.text.x=element_text(size=15), axis.title.x=element_text(size=15), axis.text.y=element_text(size=15), axis.title.y=element_text(size=15), legend.title=element_text(size=15), plot.title=element_text(size=title_size), plot.subtitle=element_text(size=20), legend.text=element_text(size=15)) +
	geom_hline(yintercept=0.95, col="dark grey") + guides(fill=guide_legend(title="Simulated ICC")) +
	scale_x_discrete(labels=c("Low Heterogeneity", "High Heterogeneity")) + scale_color_manual(values=other_cols) + scale_fill_manual(values=other_cols) +
	annotate("text", 0.8, 1, label="All Cases", size=8)


pdf("MS/fig/Coverage.pdf", height=7, width=18)

grid.arrange(A+labs(title="A."), B+labs(title="B."), C+labs(title="C."), layout_matrix=array(c(1,2,3), c(1,3)))
													
dev.off()
	
# Analysis of bias in Tau2
	
A<-ggplot(data=agg_results, aes(y=median_bias_tau2_lnR, x=Method, col=Method, fill=Method)) +
	geom_violin() + 
	geom_beeswarm(groupOnX=T, size=2, cex=0.175, col="black", alpha=1/3, shape=1) + 
	ylim(-1, 6) +
	ylab("Bias in Heterogeneity (log Ratio)") + xlab("") + theme_bw() + 
	theme(legend.position="none", axis.title.x=element_text(size=15), axis.text.x=element_text(size=15), axis.text.y=element_text(size=15), plot.title=element_text(size=title_size), axis.title.y=element_text(size=15), ) +
	geom_hline(yintercept=0, col="dark grey", size=0.5) + scale_x_discrete(labels=c("Full Data", "Missing
Cases", "All 
Cases", "Multiplicative", "Hybrid")) + scale_color_manual(values=method_cols) + scale_fill_manual(values=method_cols)

B<-ggplot(data=agg_results, aes(x=as.factor(tau2), y=median_bias_tau2_lnR, fill=Method)) +
	geom_boxplot() + theme_bw() + 
	xlab("") + ylab("Bias in Heterogeneity (log Ratio)") +
	theme(legend.position="none", axis.text.x=element_text(size=15), axis.title.x=element_text(size=15), axis.title.y=element_text(size=15), axis.text.y=element_text(size=15), plot.title=element_text(size=title_size)) +
	geom_hline(yintercept=0, col="dark grey") +
	scale_x_discrete(labels=c("Low Heterogeneity", "High Heterogeneity")) + scale_color_manual(values=method_cols) + scale_fill_manual(values=method_cols)

C<-ggplot(data=agg_results, aes(y=median_ICC, x=Method, col=Method, fill=Method)) +
	geom_violin(adjust=10) + 
	geom_beeswarm(groupOnX=T, size=2, cex=0.175, col="black", alpha=1/3, shape=1) +
	ylim(-0.5, 0.75) +
	ylab("Bias in ICC") + xlab("") + theme_bw() + 
	theme(legend.position="none", axis.title.y=element_text(size=15), axis.text.x=element_text(size=15), axis.text.y=element_text(size=15), plot.title=element_text(size=title_size)) +
	geom_hline(yintercept=0, col="dark grey", size=0.5) + scale_x_discrete(labels=c("Full Data", "Missing
Cases", "All 
Cases", "Multiplicative", "Hybrid")) + scale_color_manual(values=method_cols) + scale_fill_manual(values=method_cols)
		
D<-ggplot(data=dat1.1, aes(x=as.factor(tau2), y=median_ICC, fill=as.factor(icc_study))) +
	geom_violin() + 
	theme_bw() + 
	xlab("") + ylab("Bias in ICC") +
	theme(legend.position=c(0.84, 0.85), axis.text.x=element_text(size=15), axis.title.x=element_text(size=15), axis.text.y=element_text(size=15), axis.title.y=element_text(size=15), legend.title=element_text(size=15), plot.title=element_text(size=title_size), legend.text=element_text(size=15)) +
	geom_hline(yintercept=0, col="dark grey") + guides(fill=guide_legend(title="Simulated ICC")) +
	scale_x_discrete(labels=c("Low Heterogeneity", "High Heterogeneity")) + scale_color_manual(values=other_cols) + scale_fill_manual(values=other_cols) +
	annotate("text", 0.9, 0.5, label="Missing Cases", size=8)	
	
E<-ggplot(data=dat1.2, aes(x=as.factor(tau2), y=median_ICC, fill=as.factor(icc_study))) +
	geom_violin() + 
	theme_bw() + 
	xlab("") + ylab("Bias in ICC") +
	theme(legend.position=c(0.84, 0.85), axis.text.x=element_text(size=15), axis.title.x=element_text(size=15), axis.text.y=element_text(size=15), axis.title.y=element_text(size=15), legend.title=element_text(size=15), plot.title=element_text(size=title_size), legend.text=element_text(size=15)) +
	geom_hline(yintercept=0, col="dark grey") + guides(fill=guide_legend(title="Simulated ICC")) +
	scale_x_discrete(labels=c("Low Heterogeneity", "High Heterogeneity")) + scale_color_manual(values=other_cols) + scale_fill_manual(values=other_cols) +
	annotate("text", 0.8, 0.5, label="All Cases", size=8)	
	
pdf("MS/fig/Bias_Het.pdf", height=14, width=12)

grid.arrange(A+labs(title="A."), B+labs(title="B."), C+labs(title="C."), D+labs(title="D."), E+labs(title="E."), layout_matrix=rbind(c(1,1,1,2,2,2),
					c(6,3,3,3,3,6),
					c(4,4,4,5,5,5)))	
													
dev.off()

# % of Missing Data Plot
agg_results$prop_jitter<-(agg_results$proportion_lost + runif(nrow(agg_results), 0, 0.016)) * 100
A<-ggplot(data=agg_results, aes(x=prop_jitter, y=median_bias, col=Method, fill=Method)) + geom_point(size=0.75, alpha=1/2)	+
	geom_smooth(alpha=1/3, method="gam", se=F) +
	theme_bw() + geom_hline(yintercept=0, col="dark grey", size=0.5) + scale_color_manual(values=method_cols, labels=c("Full Data", "Missing Cases", "All Cases", "Multiplicative", "Hybrid")) + scale_fill_manual(values=method_cols, labels=c("Full Data", "Missing Cases", "All Cases", "Multiplicative", "Hybrid")) +
	xlab("% Studies with Missing SDs") + ylab("Bias lnRR") +
	theme(legend.position=c(0.88, 0.83), axis.text.x=element_text(size=15), axis.title.x=element_text(size=15), axis.text.y=element_text(size=15), axis.title.y=element_text(size=15), legend.title=element_text(size=0), plot.title=element_text(size=title_size), legend.text=element_text(size=14), legend.key=element_rect(fill="transparent"), legend.background=element_rect(fill="transparent")) +
	ylim(-0.003, 0.009)


B<-ggplot(data=agg_results, aes(x=prop_jitter, y=coverage, col=Method, fill=Method)) + geom_point(size=0.75, alpha=1/2)	+
	geom_smooth(alpha=1/3, method="gam", se=F) +
	theme_bw() + geom_hline(yintercept=0.95, col="dark grey", size=0.5) + scale_color_manual(values=method_cols, labels=c("Full Data", "Missing Cases", "All Cases", "Multiplicative", "Hybrid")) + scale_fill_manual(values=method_cols, labels=c("Full Data", "Missing Cases", "All Cases", "Multiplicative", "Hybrid")) +
	xlab("% Studies with Missing SDs") + ylab("Coverage (95% CI)") +
	theme(legend.position="none", axis.text.x=element_text(size=15), axis.title.x=element_text(size=15), axis.text.y=element_text(size=15), axis.title.y=element_text(size=15), plot.title=element_text(size=title_size)) +
	ylim(0.925, 1)
B


C<-ggplot(data=agg_results, aes(x=prop_jitter, y=median_bias_tau2_lnR, col=Method, fill=Method)) + geom_point(size=0.75, alpha=1/2)	+
	geom_smooth(alpha=1/3, method="gam", se=F) +
	theme_bw() + geom_hline(yintercept=0, col="dark grey", size=0.5) + scale_color_manual(values=method_cols, labels=c("Full Data", "Missing Cases", "All Cases", "Multiplicative", "Hybrid")) + scale_fill_manual(values=method_cols, labels=c("Full Data", "Missing Cases", "All Cases", "Multiplicative", "Hybrid")) +
	xlab("% Studies with Missing SDs") + ylab("Bias in Heterogeneity (log Ratio)") +
	theme(legend.position="none", axis.text.x=element_text(size=15), axis.title.x=element_text(size=15), axis.text.y=element_text(size=15), axis.title.y=element_text(size=15), plot.title=element_text(size=title_size)) +
	ylim(-0.5, 8)
	
pdf("MS/fig/Prop_Missing.pdf", height=15, width=10)

grid.arrange(A+labs(title="A."), B+labs(title="B."), C+labs(title="C."), layout_matrix=rbind(1,2,3))	
													
dev.off()






# Simulation to Support Nakagawa et al. A practical and readily implementable method for handling missing standard deviation in the meta-analysis of response ratios. 

# This script includes is for plotting the results of the aggregated simulation

# Clean up the R Environment 
rm(list=ls())

# Set the working directory
wd<-"/Users/alistairsenior/Dropbox (Sydney Uni)/Nakagawa_Ecology_Missing_SD/Miss_SD_Sim" # MacMini Home
setwd(wd)

# Load the relevant libraries, and header
library(ggplot2)
library(plyr)
library(gridExtra)
library(ggbeeswarm)
library(ggcorrplot)

# Load the results
load("agg_results.Rdata")

# Create a folder for the plots
if(file.exists("Plots_MLMA")){
	unlink("Plots_MLMA", recursive=T)
}
dir.create("Plots_MLMA")
setwd("Plots_MLMA")

########################################################################
############## Averaging each parameter set ############################
########################################################################

# Analysis of Bias in the lnRR

head(agg_results)
agg_results$range_bias<-agg_results$max_bias - agg_results$min_bias

A<-ggplot(data=agg_results, aes(x=mean_bias, y=Method, col=Method, fill=Method)) +
	geom_violin() + 
	geom_beeswarm(groupOnX=F, size=0.05, cex=0.25, col="black") + 
	xlim(-0.003, 0.008) +
	xlab("Bias lnRR") + ylab("") + theme_bw() + 
	theme(legend.position="none", axis.title.x=element_text(size=15), axis.text.y=element_text(size=15), plot.title=element_text(size=15)) +
	geom_vline(xintercept=0, col="dark grey", size=0.5)


# # B<-ggplot(data=agg_results, (aes(x=as.factor(proportion_lost), y=mean_bias, fill=Method, col=Method))) + 
	# geom_violin() +
	# theme_bw() +
	# geom_hline(yintercept=0, col="grey") +
	# xlab("Proportion Missing") + ylab("Bias lnRR") +
	# theme(legend.position="none", axis.title.x=element_text(size=15), axis.title.y=element_text(size=15))


long<-agg_results[,c("Method", "mean_bias", "code")]
wide<-reshape(long, direction="wide", idvar="code", timevar="Method")
cor.mat<-cor(wide[,-1])
rownames(cor.mat)<-c("Complete", "Method 1.1", "Method 1.2", "Method 2", "Method 3")
colnames(cor.mat)<-c("Complete", "Method 1.1", "Method 1.2", "Method 2", "Method 3")

B<-ggcorrplot(cor.mat, lab=T, show.legend=FALSE, type="lower") + theme(plot.title=element_text(size=15))

# Difference in the bias
df<-data.frame(delta=abs(wide[,"mean_bias.Method 1.1"]) - abs(wide[,"mean_bias.Method 1.2"]))
C<-ggplot(df, aes(x=delta)) + geom_histogram(bins=20) + theme_bw() + geom_vline(xintercept=0, col="red") +
	xlab("Diff. |Bias| Meths 1.1 & 1.2") + ylab("Frequency")

D<-ggplot(data=agg_results, aes(x=log(range_bias, base=10), y=Method, col=Method, fill=Method)) +
	geom_violin() + 
	geom_beeswarm(groupOnX=FALSE, size=0.05, cex=0.25, col="black") + 
	xlab("log10 Range of Bias lnRR") + ylab("") + theme_bw() + 
	theme(legend.position="none", axis.title.x=element_text(size=15), axis.text.y=element_text(size=15), plot.title=element_text(size=15))

# Get the data for Methods 1.1 & 1.2
dat1.1<-agg_results[which(agg_results$Method=="Method 1.1"),]
dat1.2<-agg_results[which(agg_results$Method=="Method 1.2"),]

E<-ggplot(dat1.2, aes(x=log(range_bias, base=10), y=as.factor(sd_study_sd), fill=as.factor(n_study_mu))) +
	geom_violin() +
	geom_point(size=-1, col="grey") +
	theme_bw() +
	xlab("log10 Range of Bias lnRR") +
	ylab(expression(sigma[s])) + 
	theme(axis.title.x=element_text(size=15), axis.text.y=element_text(size=15), axis.title.y=element_text(size=15), legend.position=c(0.8, 0.2), legend.text=element_text(size=15), legend.title=element_text(size=15), plot.title=element_text(size=15)) + guides(fill=guide_legend(title="Study Sample Size"))	



pdf("Bias_lnRR.pdf", height=10, width=12)

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
	scale_x_discrete(labels=c(expression(9e-06~~or~~over(italic(T), lnRR)==0.01), expression(0.09~~or~~over(italic(T), lnRR)==1))) + labs(subtitle="Method 1.1")

C<-ggplot(data=dat1.2, aes(x=as.factor(tau2), y=coverage, fill=as.factor(icc_study))) +
	geom_violin() + 
	theme_bw() + 
	xlab(expression(italic(T)^2)) + ylab("Coverage (95% CI)") +
	ylim(0.9, 1) +
	theme(legend.position=c(0.85, 0.875), axis.text.x=element_text(size=12), axis.title.x=element_text(size=15), axis.text.y=element_text(size=12), axis.title.y=element_text(size=15), legend.title=element_text(size=15), plot.title=element_text(size=15), plot.subtitle=element_text(size=15), legend.text=element_text(size=15)) +
	geom_hline(yintercept=0.95, col="dark grey") + guides(fill=guide_legend(title="Simulated ICC")) +
	scale_x_discrete(labels=c(expression(9e-06~~or~~over(italic(T), lnRR)==0.01), expression(0.09~~or~~over(italic(T), lnRR)==1))) + labs(subtitle="Method 1.2")

pdf("Coverage.pdf", height=7, width=18)

grid.arrange(A+labs(title="A."), B+labs(title="B."), C+labs(title="C."), layout_matrix=array(c(1,2,3), c(1,3)))
													
dev.off()
	
# Analysis of bias in Tau2
	
A<-ggplot(data=agg_results, aes(x=mean_bias_tau2_lnR, y=Method, col=Method, fill=Method)) +
	geom_violin() + 
	geom_beeswarm(groupOnX=F, size=0.05, cex=0.25, col="black") + 
	xlim(-6, 6) +
	xlab("Bias in Heterogeneity (log Ratio)") + ylab("") + theme_bw() + 
	theme(legend.position="none", axis.title.x=element_text(size=15), axis.text.x=element_text(size=12), axis.text.y=element_text(size=15), plot.title=element_text(size=15)) +
	geom_vline(xintercept=0, col="dark grey", size=0.5)

B<-ggplot(data=agg_results, aes(x=as.factor(tau2), y= mean_bias_tau2_lnR, fill=Method)) +
	geom_boxplot() + theme_bw() + 
	xlab(expression(italic(T)^2)) + ylab("Bias in Heterogeneity (log Ratio)") +
	theme(legend.position="none", axis.text.x=element_text(size=12), axis.title.x=element_text(size=15), axis.text.y=element_text(size=12), axis.title.y=element_text(size=15), plot.title=element_text(size=15)) +
	geom_hline(yintercept=0, col="dark grey") +
	scale_x_discrete(labels=c(expression(9e-06~~or~~over(italic(T), lnRR)==0.01), expression(0.09~~or~~over(italic(T), lnRR)==1)))

C<-ggplot(data=agg_results, aes(x=mean_bias_ICC, y=Method, col=Method, fill=Method)) +
	geom_violin() + 
	geom_beeswarm(groupOnX=F, size=0.05, cex=0.25, col="black") + 
	xlim(-0.5, 0.75) +
	xlab("Bias in ICC") + ylab("") + theme_bw() + 
	theme(legend.position="none", axis.title.x=element_text(size=15), axis.text.x=element_text(size=12), axis.text.y=element_text(size=15), plot.title=element_text(size=15)) +
	geom_vline(xintercept=0, col="dark grey", size=0.5)
		
D<-ggplot(data=dat1.2, aes(x=as.factor(tau2), y= mean_bias_ICC, fill=as.factor(icc_study))) +
	geom_violin() + 
	theme_bw() + 
	xlab(expression(italic(T)^2)) + ylab("Bias in ICC") +
	theme(legend.position=c(0.85, 0.85), axis.text.x=element_text(size=12), axis.title.x=element_text(size=15), axis.text.y=element_text(size=12), axis.title.y=element_text(size=15), legend.title=element_text(size=15), plot.title=element_text(size=15), legend.text=element_text(size=15)) +
	geom_hline(yintercept=0, col="dark grey") + guides(fill=guide_legend(title="Simulated ICC")) +
	scale_x_discrete(labels=c(expression(9e-06~~or~~over(italic(T), lnRR)==0.01), expression(0.09~~or~~over(italic(T), lnRR)==1)))	
	

pdf("Bias_Het.pdf", height=10, width=12)

grid.arrange(A+labs(title="A."), B+labs(title="B."), C+labs(title="C."), D+labs(title="D."), layout_matrix=rbind(c(1,2),
					c(3,4)))	
													
dev.off()	
	
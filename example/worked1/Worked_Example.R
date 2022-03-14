
# Clear out the R-Environment
rm(list=ls())

# Load packages
library(metafor)
library(plyr)

# Function to calculate the effect sizes both ways
my_calc_es<-function(c_mean, c_sd, c_n, t_mean, t_sd, t_n, study_id){

	# Load package locally
	require(plyr)

	# bind all together
	data<-data.frame(c_mean=c_mean, c_sd=c_sd, c_n=c_n, t_mean=t_mean, t_sd=t_sd, t_n=t_n, study_id=study_id)

	# Calculate the CVs
	data$c_cv<-data$c_sd/data$c_mean
	data$t_cv<-data$t_sd/data$t_mean

	# Get the weighted CV^2
	# First pool at the level of the study, weighted
	pooled<-ddply(data, .(study_id), summarise, c_cv2=weighted.mean(c_cv^2, c_n, na.rm=T), t_cv2=weighted.mean(t_cv^2, t_n, na.rm=T), c_n=mean(c_n, na.rm=T), t_n=mean(t_n, na.rm=T))

	# Now get the weighted cv2s
	mean_c_cv2<-weighted.mean(pooled$c_cv2, pooled$c_n, na.rm=T)
	mean_t_cv2<-weighted.mean(pooled$t_cv2, pooled$t_n, na.rm=T)

	# Calculate the effect sizes and variances different ways
	yi_miss<-log(data$t_mean / data$c_mean) + 0.5 * (data$t_cv^2 / data$t_n - data$c_cv^2 / data$c_n)
	vi_miss<-data$t_cv^2 / data$t_n + data$c_cv^2 / data$c_n + data$t_cv^4 / (2 * data$t_n^2) + data$c_cv^4 / (2 * data$c_n^2)

	yi_1B<-log(data$t_mean / data$c_mean) + 0.5 * (mean_t_cv2 / data$t_n - mean_c_cv2 / data$c_n)
	vi_1B<-mean_t_cv2 / data$t_n + mean_c_cv2 / data$c_n + mean_t_cv2^2 / (2 * data$t_n^2) + mean_c_cv2^2 / (2 * data$c_n^2)

	yi_1A<-yi_miss
	vi_1A<-vi_miss

	missing<-which(is.na(data$c_sd) == T)
	yi_1A[missing]<-yi_1B[missing]
	vi_1A[missing]<-vi_1B[missing]


	# Caluclate the matrix for method 2 - note this is just method 1B in the diagonal
	Vf<-diag(vi_1B)
	row.names(Vf)<-seq(1, length(yi_1A), 1)
	colnames(Vf)<-seq(1, length(yi_1A), 1)


	# Package up and return as list - first object is effect sizes, second is matrix for method 2
	output<-list()
	output[[1]]<-data.frame(ES_ID = seq(1, length(yi_1A), 1), yi_miss = yi_miss, vi_miss = vi_miss, yi_1A = yi_1A, vi_1A = vi_1A, yi_1B = yi_1B, vi_1B = vi_1B)
	output[[2]]<-Vf
	return(output)
}


# Load the data
data<-read.csv("./example/worked1/data_test.csv")
head(data)

# Lets follow their lead and split by grazing type as per the paper. Here I will only re-analyse the CG group
data<-data[which(data$Comparison == "CG-SRG"),]

# Double check some stuff

# How many effect sizes
dim(data)[1]
# 173

# How many studies
length(unique(data$Study))
# 67 Studies

# How many grazing organisms
length(unique(data$Stock.type))
# 4, which are
unique(data$Stock.type)

# Is there a common control issue
data$control_ID<-paste0(data$Study, "_", data$Common.control)
length(unique(data$control_ID))
length((data$control_ID))
# Yes some studies have the issue

# Which are they
check<-ddply(data, .(Study), summarise, n_controls=length(unique(control_ID)), n_effects=length(control_ID))
shared<-check[which(check$n_controls != check$n_effects),]
shared
dim(shared)[1] / dim(check)[1]
# 13.4% have a shared control issue to some degree

# Double check they really are shared control - should be no variance within shared control
check<-ddply(data, .(control_ID), summarise, length(CN), sd(CN, na.rm=T), sd(CM, na.rm=T))
check[which(check[,2] > 1),]
# Yep looks good

# Can also check number of effect sizes with shared control issue here
sum(check[which(check[,2] > 1),2]) / dim(data)[1]
# 13.9% of effect sizes

# How much missing SD data?
plot(is.na(data$CSD), is.na(data$TSD))
# If control is missing, so is treatment - always here

sum(is.na(data$CSD)) / dim(data)[1]
# 35.8% of effects sizes have missing sds

# Function to calculate the effect sizes the four different ways

# Originally Macdonald et al handle the missing SD data by calculating the mean CV (for each grazing type, CG-SRG and UG-SRG) for the present data, then use this with the known mean to impute the SD. Essentially impute as a function of the mean assuming a linear assocation between mean and SD. This is very similar to our method 1A, but without weighted average and not on CV^2

# Now calculate effect sizes
data<-cbind(data, my_calc_es(c_mean=data$CM, c_sd=data$CSD, c_n=data$CN, t_mean=data$TM, t_sd=data$TSD, t_n=data$TN, study_id=data$Study)[[1]])
data$ES_ID<-as.factor(seq(1, nrow(data), 1))

# Model 1A: MLMA of missing SD effect sizes using pooled CV2
MLMA1A<-rma.mv(yi = yi_1A, V = vi_1A, random=list(~1|Study, ~1|ES_ID), data=data)
summary(MLMA1A)

# Model 1B: MLMA of all effect sizes using pooled CV2
MLMA1B<-rma.mv(yi = yi_1B, V = vi_1B, random=list(~1|Study, ~1|ES_ID), data=data)
summary(MLMA1B)

# Model 2: weighted regression of all effect sizes using pooled CV2
# Note here and below I analyse here with effect sizes estimated as the mix of known CV and mean CV - doesn't have to be though
Vf<-my_calc_es(c_mean=data$CM, c_sd=data$CSD, c_n=data$CN, t_mean=data$TM, t_sd=data$TSD, t_n=data$TN, study_id=data$Study)[[2]]
data$ES_ID2<-rownames(Vf)
MLMA2<-rma.mv(yi = yi_1A, V = 0, random=list(~1|Study, ~1|ES_ID, ~1|ES_ID2), data=data, R=list(ES_ID2=Vf), Rscale=F)
summary(MLMA2)

# Model 3: add 0s in the matrix where SDs are not missing, and to vi where they are
missing<-which(is.na(data$yi_miss) == T)
diag(Vf)[-missing]<-0
data$vi_3<-data$vi_miss
data$vi_3[missing]<-0
MLMA3<-rma.mv(yi = yi_1A, V = vi_3, random=list(~1|Study, ~1|ES_ID, ~1|ES_ID2), data=data, R=list(ES_ID2=Vf), Rscale=F)
summary(MLMA3)

# Package up as table
results<-data.frame(method=c("1A", "1B", "2", "3"), b=c(MLMA1A$b, MLMA1B$b, MLMA2$b, MLMA3$b), ci.lb=c(MLMA1A$ci.lb, MLMA1B$ci.lb, MLMA2$ci.lb, MLMA3$ci.lb), ci.ub=c(MLMA1A$ci.ub, MLMA1B$ci.ub, MLMA2$ci.ub, MLMA3$ci.ub), tau2=c(sum(MLMA1A$sigma2), sum(MLMA1B$sigma2), sum(MLMA2$sigma2), sum(MLMA3$sigma2)))
write.table(results, file="Example.results.csv", sep=",", row.names=F, col.names=names(results))


############################
# Load libraries
###########################
rm(list=ls())
pacman::p_load(ape, metafor, MCMCglmm, tidyverse, readxl)
source("./R/func.R")

############################
# Loading and sorting data out
###########################
#load data
dat  <- data.frame(read_excel("./example/worked2/ele13245-sup-0003-metas3.xlsx"))

# Fix column classes
dat$Experimental_standard_deviation <- as.numeric(dat$Experimental_standard_deviation)
dat$Control_standard_deviation <- as.numeric(dat$Control_standard_deviation)

#remove rows with missing data
    dat <- na.omit(dat)
dat$obs <- 1:nrow(dat) # Needed for metafor to make model equivalent to MCMCglmm

#load phylogenetic tree
tree<- read.tree("./example/worked2/ele13245-sup-0008-phylogenys8.tre")

#######################
    # Example 2
#######################
# Intercept only model. for Abundance fitness ONLY. The authors state: "We first estimated the average impact of insect interaction on each fitness component by fitting linear models in which the intercept was the only fixed effect." Table 1 provides the fitness components. Abundance has the most data.
    a2 <- dat %>% filter(Fitness_component == "Abundance") # Doesn't quite match Table 1, but does match figure 1, and even exclusion N = 311, which matches.

    # Here, we want to use lnRR so we need to calculate this
    a2 <- escalc( m1i = Control_mean,
                  m2i = Experimental_mean,
                 sd1i = Control_standard_deviation,
                 sd2i = Experimental_standard_deviation,
                  n1i = Control_sample_size,
                  n2i = Experimental_sample_size,
                 measure = "ROM", var.names=c("yi_lnrr","vi_lnrr"),
                 append = TRUE, data = a2)
    # Some NA's because lnRR, not too many though 311 -> 306
    a2 <- na.omit(a2)

    # Prune tree
    #align focal species with branches of phylogenetic tree
    tree2_meta<-drop.tip(tree,
                         tree$tip.label[which(tree$tip.label%in%a2$Focal_insect==FALSE)])
    check.species<-function(x) {any(x==tree2_meta$tip.label)}
    a2 <- a2[sapply(a2[,"Focal_insect"],check.species),]

    tree_checks(data = a2, tree = tree2_meta, species_name_col = "Focal_insect", type = "checks")

    phylo <- vcv(tree2_meta, corr = TRUE)

################################################
    # Generate missing data
################################################
    # Create missingness at the study level. It's more likely a study doesn't present SD than a subset of effects within a study.
    # Important to note that, how much of an impact missing data will have really depends on the sample of studies,
    # whether it's random or not.
    set.seed(675)
    stdies <- sample(unique(a2$Author), size = 0.2*(length(unique(a2$Author))))
    a2missSD_stdy <- a2
    a2missSD_stdy[which(a2missSD_stdy$Author %in% stdies), c("Experimental_standard_deviation", "Control_standard_deviation")] <- NA

    # Seems to be some issue data when calculating sampling variances with CV, which suggest some raw data maybe wrong.
    a2missSD_stdy <- a2missSD_stdy[!a2missSD_stdy$ID%in%c("833", "1176"),]

    # Now, assume you needto exclude data with missing SD because you can't calculate effect size and sampling variance
    complete_case_MV <- na.omit(a2missSD_stdy)

    # Prune tree.
    tree3_meta<-drop.tip(tree,
                         tree$tip.label[which(tree$tip.label%in%complete_case_MV$Focal_insect==FALSE)])
    check.species<-function(x) {any(x==tree3_meta$tip.label)}
    complete_case_MV <- complete_case_MV[sapply(complete_case_MV[,"Focal_insect"],check.species),]

    #Check that it was done right
    tree_checks(data = complete_case_MV, tree = tree3_meta,
                species_name_col = "Focal_insect", type = "checks")

    # Make new VCV matrix
    phylo2 <- vcv(tree3_meta, corr = TRUE)

################################################
# Whole/full data model
################################################

    whole_mv <- rma.mv(yi_lnrr ~ 1, V = vi_lnrr,
                          random=list(~1|Group, ~1|Year, ~1|Focal_insect, ~1|obs),
                          R = list(Focal_insect = phylo), data = a2)
    whole_mv_res <- get_est(whole_mv)


################################################
# Complete case analysis
################################################
    # Fit complete case analysis. Note that data is currently missing at random.
    complete_case_mv <- rma.mv(yi_lnrr ~ 1, V = vi_lnrr,
                          random=list(~1|Group, ~1|Year, ~1|Focal_insect, ~1|obs),
                          R = list(Focal_insect = phylo2), data = complete_case_MV)

    complete_case_mv_res <- get_est(complete_case_mv)

################################################
    # METHOD 1A
################################################
    # Spake and Doncaster Method that takes CV across studies
    # FIrst calculate CV on missing dataset. Note missing data will be ignored
    a2missSD_stdy <- a2missSD_stdy %>%
                 mutate(     cv_Control = na_if(Control_mean / Control_standard_deviation, Inf),
                        cv_Experimental = na_if(Experimental_mean / Experimental_standard_deviation, Inf))

     # Now calculate the average between study CV, which will replace missing values.
     # Note that mean N used
    a2missSD_stdy <- cv_avg(cv_Control, Control_sample_size, group = Author,
                            name = "1", data = a2missSD_stdy)
    a2missSD_stdy <- cv_avg(cv_Experimental, Experimental_sample_size, group = Author,
                            name = "2", data = a2missSD_stdy)

    # Now using wighted mean CV in replacement for where CV's are missing.
    # Note that function above already CV^2 so need to do that on original CV
    a2missSD_stdy <- a2missSD_stdy %>%
                      mutate(cv_cont_new = if_else(is.na(cv_Control),      b_CV_1, cv_Control^2),
                             cv_expt_new = if_else(is.na(cv_Experimental), b_CV_2, cv_Experimental^2))

    # Now calculate new vi, called vi_DS_lnrr. Note that CV is alreday ^2
    a2missSD_stdy <- a2missSD_stdy %>%
                      mutate(vi_DS_lnrr = (cv_cont_new / Control_sample_size) + (cv_expt_new / Experimental_sample_size))

    # Fit model with new sampling variance
    method_1A_mv <- rma.mv(yi_lnrr ~ 1, V = vi_DS_lnrr,
                               random=list(~1|Group, ~1|Year, ~1|Focal_insect, ~1|obs),
                               R = list(Focal_insect = phylo), data = a2missSD_stdy)

    method_1A_mv_res <- get_est(method_1A_mv)
################################################
    # METHOD 1B
################################################
    lnrr_laj <- function(m1, m2, cv1_2, cv2_2, n1, n2){
            log(m1 / m2) + 0.5*((cv1_2 / n1) - (cv2_2 / n2))
    }

    v_lnrr_laj <- function(cv1_2, cv2_2, n1, n2){
          ((cv1_2) / n1) + ((cv2_2) / n2) +
        ((cv1_2)^2 / (2*n1)^2) + ((cv2_2)^2 / (2*n2)^2)
    }

    # Now calculate new yi andvi, called lnrr_laj & v_lnrr_laj, respectively.
    a2missSD_stdy <- a2missSD_stdy %>%
      mutate(lnrr_laj = lnrr_laj(m1 = Control_mean, m2 = Experimental_mean, cv1 = cv_cont_new, cv2 = cv_expt_new,
                                 n1= Control_sample_size, n2 = Experimental_sample_size),
           v_lnrr_laj = v_lnrr_laj(cv1 = cv_cont_new, n1= Control_sample_size,
                                   cv2 = cv_expt_new, n2 = Experimental_sample_size))

    # Fit model with new sampling variance
    method_1B_mv <- rma.mv(lnrr_laj ~ 1, V = v_lnrr_laj,
                           random=list(~1|Group, ~1|Year, ~1|Focal_insect, ~1|obs),
                           R = list(Focal_insect = phylo), data = a2missSD_stdy)

    method_1B_mv_res <- get_est(method_1B_mv)
################################################
    # METHOD 2
################################################
                 V_es <- diag(a2missSD_stdy$v_lnrr_laj)
      row.names(V_es) <- a2missSD_stdy$obs
   a2missSD_stdy$obs2 <- rownames(V_es)

    method2_mv <-rma.mv(lnrr_laj ~ 1, V = 0, random=list(~1|Group, ~1|Year, ~1|Focal_insect, ~1|obs, ~1|obs2),
                        data=a2missSD_stdy, R=list(obs2=V_es), Rscale=F)
    method2_mv_res <- get_est(method2_mv)

################################################
    # METHOD 3
################################################
    # Find where missing SD's are in the data
      missing_dat <- which(is.na(a2missSD_stdy$Control_standard_deviation) & is.na(a2missSD_stdy$Experimental_standard_deviation))

    # Set the effect sizes not missing data to zero in the V_es matrix
                        V_es2 <- diag(a2missSD_stdy$v_lnrr_laj)
    diag(V_es2)[-missing_dat] <- 0
             row.names(V_es2) <- a2missSD_stdy$obs
           a2missSD_stdy$obs2 <- rownames(V_es2)

    # Set the v_lnrr_laj to 0
                 a2missSD_stdy$v_lnrr_laj_m3 <- a2missSD_stdy$v_lnrr_laj
    a2missSD_stdy$v_lnrr_laj_m3[missing_dat] <- 0

    method3_mv <-rma.mv(lnrr_laj ~ 1, V = v_lnrr_laj_m3, random=list(~1|Group, ~1|Year, ~1|Focal_insect, ~1|obs, ~1|obs2),
                        data=a2missSD_stdy, R=list(obs2=V_es2), Rscale=F)
    method3_mv_res <- get_est(method3_mv)

################################################
    # Results Table
################################################
    results <- rbind(whole_mv_res, complete_case_mv_res, method_1A_mv_res, method_1B_mv_res, method2_mv_res, method3_mv_res)
    row.names(results) <- c("Whole Data", "Complete Case", "Method 1A", "Method 1B", "Method 2", "Method 3")
    write.csv(results, "./example/worked2/results2.csv", row.names = TRUE)

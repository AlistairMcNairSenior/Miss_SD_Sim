################################################
    # Load libraries
################################################
rm(list=ls())
pacman::p_load(ape, metafor, MCMCglmm, tidyverse, readxl)
source("./R/func.R")

################################################
    # Loading and sorting data out
################################################
#load data
dat  <- data.frame(read_excel("./example/worked1/ele13245-sup-0003-metas3.xlsx"))

# Fix column classes
dat$Experimental_standard_deviation <- as.numeric(dat$Experimental_standard_deviation)
dat$Control_standard_deviation <- as.numeric(dat$Control_standard_deviation)

#remove rows with missing data
dat$obs <- 1:nrow(dat) # Needed for metafor to make model equivalent to MCMCglmm

# Clean up the data to drop unnecessary columns
dat <- dat[,-which(colnames(dat)%in%c("vi", "yi_g", "Focal_insect_resolved_name", "Competing_insect_resolved_name", "Focal_insect_diet_breadth", "Competing_insect_diet_breadth", "Phylogenetic_distance_between_focal_insect_and_competing_insect", "Spatial_separation", "temporal_separation"))]

#load phylogenetic tree
tree<- read.tree("./example/worked1/ele13245-sup-0008-phylogenys8.tre")

################################################
    # Greary's Test, Calculate effect sizes
    # and prune tree
################################################
  # Intercept only model. for Abundance fitness ONLY. The authors state: "We first estimated the average impact of insect interaction on each fitness component by fitting linear models in which the intercept was the only fixed effect." Table 1 provides the fitness components. Abundance has the most data.
    a2 <- dat %>% filter(Fitness_component == "Abundance")

  # First calculate CV on missing dataset. Note missing data will be ignored
           a2 <- a2 %>%
                  mutate(cv_Control = na_if(Control_standard_deviation / Control_mean, Inf),
                    cv_Experimental = na_if(Experimental_standard_deviation / Experimental_mean, Inf))


    # Add in Laj lnRR correction. Note we need to ^2 cv's here
    a2 <-  a2 %>%
  mutate(lnrr_laj_orig = na_if(lnrr_laj(m1 = Control_mean, m2 = Experimental_mean,
                                      cv1_2 = cv_Control^2, cv2_2  = cv_Experimental^2,
                              n1= Control_sample_size, n2 = Experimental_sample_size), Inf),
        v_lnrr_laj_orig = na_if(v_lnrr_laj(cv1_2 = cv_Control^2, n1= Control_sample_size,
                              cv2_2 = cv_Experimental^2, n2 = Experimental_sample_size), Inf))
################################################
    # Geary's test
################################################
    # Function to calculate Geary's "number"
    geary <- function(mean, sd, n){
      (1 / (sd / mean)) * ((4*n)^(3/2) / (1 + 4*n))
    }

    # Geary's test; assumption of normality assumed to be approximately correct when values are >= 3.
    a2 <- a2 %>%
      mutate(geary_control = geary(Control_mean, Control_standard_deviation, Control_sample_size),
             geary_trt = geary(Experimental_mean, Experimental_standard_deviation, Experimental_sample_size),
             geary_test = ifelse(geary_control >= 3 & geary_trt >= 3, "pass", "fail"))

    # Exclude data failing test
       a2  <- a2  %>%
      filter(geary_test == "pass")

################################################
    # Phylogenetic test
################################################
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
    a2missSD_stdy[which(a2missSD_stdy$Author %in% stdies),
                  c("Experimental_standard_deviation", "Control_standard_deviation")] <- NA

    # First calculate CV on missing dataset. Note missing data will be ignored
    a2missSD_stdy <- a2missSD_stdy %>%
      mutate(cv_Control = na_if(Control_standard_deviation / Control_mean, Inf),
             cv_Experimental = na_if(Experimental_standard_deviation / Experimental_mean, Inf))


    # Now, assume you need to exclude data with missing SD because you can't calculate effect size and sampling variance
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


    # Now calculate the average between study CV, which will replace missing values.
    # Note that mean N used
    a2missSD_stdy <- cv_avg(x = Control_mean, sd = Control_standard_deviation,
                            n = Control_sample_size, group = Author, label = "1", cv2=FALSE,
                             data = a2missSD_stdy)
    a2missSD_stdy <- cv_avg(x = Experimental_mean, sd = Experimental_standard_deviation,
                            n = Experimental_sample_size, group = Author,
                            label = "2", cv2=FALSE, data = a2missSD_stdy)

    # Now using weighted mean CV in replacement for where CV's are missing.
    # Note that function above already CV^2 so need to do that on original CV
    a2missSD_stdy <- a2missSD_stdy %>%
      mutate(cv2_cont_new = if_else(is.na(cv_Control),      b_CV2_1, cv_Control^2),
             cv2_expt_new = if_else(is.na(cv_Experimental), b_CV2_2, cv_Experimental^2))


    # Calculate the new unbiased lnRR

    # Now calculate new yi and vi, called lnrr_laj & v_lnrr_laj, respectively.Note that the functions take CV^2
    a2missSD_stdy <- a2missSD_stdy %>%
      mutate(lnrr_laj = lnrr_laj(m1 = Control_mean, m2 = Experimental_mean,
                                 cv1 = cv2_cont_new, cv2 = cv2_expt_new,
                                 n1= Control_sample_size, n2 = Experimental_sample_size),
             v_lnrr_laj = v_lnrr_laj(cv1 = cv2_cont_new, n1= Control_sample_size,
                                     cv2 = cv2_expt_new, n2 = Experimental_sample_size))


################################################
    # Whole/full data model
################################################

    whole_mv <- rma.mv(lnrr_laj_orig ~ 1, V = v_lnrr_laj_orig,
                          random=list(~1|Group, ~1|Year, ~1|Focal_insect, ~1|obs),
                          R = list(Focal_insect = phylo), data = a2)
    whole_mv_res <- get_est(whole_mv)

################################################
    # Complete case analysis
################################################
    # Fit complete case analysis. Note that data is currently missing at random.
    complete_case_mv <- rma.mv(lnrr_laj_orig ~ 1, V = v_lnrr_laj_orig,
                          random=list(~1|Group, ~1|Year, ~1|Focal_insect, ~1|obs),
                          R = list(Focal_insect = phylo2), data = complete_case_MV)

    complete_case_mv_res <- get_est(complete_case_mv)

################################################
    # METHOD 1A
################################################

    # Fit model with new sampling variance
    method_1A_mv <- rma.mv(lnrr_laj ~ 1, V = v_lnrr_laj,
                               random=list(~1|Group, ~1|Year, ~1|Focal_insect, ~1|obs),
                               R = list(Focal_insect = phylo), data = a2missSD_stdy)

    method_1A_mv_res <- get_est(method_1A_mv)

################################################
    # METHOD 1B
################################################
    a2missSD_stdy <- a2missSD_stdy %>%
      mutate(v_lnrr_laj_1B = v_lnrr_laj(cv1 = b_CV2_1, n1= Control_sample_size,
                                        cv2 = b_CV2_2, n2 = Experimental_sample_size))

    # Fit model with new sampling variance and point estimate
    method_1B_mv <- rma.mv(lnrr_laj ~ 1, V = v_lnrr_laj_1B,
                           random=list(~1|Group, ~1|Year, ~1|Focal_insect, ~1|obs),
                           R = list(Focal_insect = phylo), data = a2missSD_stdy)

    method_1B_mv_res <- get_est(method_1B_mv)

################################################
    # METHOD 2
################################################
                 V_es <- diag(a2missSD_stdy$v_lnrr_laj_1B)
      row.names(V_es) <- a2missSD_stdy$obs
   a2missSD_stdy$obs2 <- rownames(V_es)

    method2_mv <-rma.mv(lnrr_laj ~ 1, V = 0, random=list(~1|Group, ~1|Year,
                                                         ~1|Focal_insect, ~1|obs, ~1|obs2),
                        data=a2missSD_stdy, R=list(obs2=V_es, Focal_insect = phylo),
                        Rscale=F)
    method2_mv_res <- get_est(method2_mv)

################################################
    # METHOD 3
################################################
    # Find where missing SD's are in the data
      missing_dat <- which(is.na(a2missSD_stdy$Control_standard_deviation) & is.na(a2missSD_stdy$Experimental_standard_deviation))

    # Set the effect sizes not missing data to zero in the V_es matrix
                        V_es2 <- diag(a2missSD_stdy$v_lnrr_laj_1B)
    diag(V_es2)[-missing_dat] <- 0
             row.names(V_es2) <- a2missSD_stdy$obs
           a2missSD_stdy$obs2 <- rownames(V_es2)

    # Set the v_lnrr_laj to 0
                 a2missSD_stdy$v_lnrr_laj_m3 <- a2missSD_stdy$v_lnrr_laj
    a2missSD_stdy$v_lnrr_laj_m3[missing_dat] <- 0

    method3_mv <-rma.mv(lnrr_laj ~ 1, V = v_lnrr_laj_m3, random=list(~1|Group, ~1|Year, ~1|Focal_insect, ~1|obs, ~1|obs2),
                        data=a2missSD_stdy, R=list(obs2=V_es2, Focal_insect = phylo), Rscale=F)
    method3_mv_res <- get_est(method3_mv)

################################################
    # Results Table
################################################
    results <- rbind(whole_mv_res, complete_case_mv_res, method_1A_mv_res, method_1B_mv_res, method2_mv_res, method3_mv_res)
    row.names(results) <- c("Whole Data", "Complete Case", "Method 1A", "Method 1B", "Method 2", "Method 3")
    write.csv(round(results, 3), "./example/worked1/results2.csv", row.names = TRUE)

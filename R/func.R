#' @title cv_avg
#' @description Calculates the weighted average CV^2 within a study and the weighted average CV^2 across a study
#' @param x Mean of an experimental group
#' @param sd Standard deviation of an experimental group
#' @param n The sample size of an experimental group
#' @param group Study, grouping or cluster variable one wishes to calculate the within and between weighted CV^2. In meta-analysis this will most likely be 'study'.
#' @param data The dataframe containing the mean, sd, n and grouping variables
#' @param label A character string specifying the label one wishes to attach to columns to identify the treatment. Otherwise, if not specified it will default to the variable name for x
#' @param sub_b A logical indicating whether the between study CV^2 (b_CV2) should be appended to the data only ('TRUE') or whether both within study CV^2 (w_CV2), mean sample size (n_mean) and between study CV^2 (b_CV2) should all be appended to the data only ('FALSE')
#' @example \dontrun{
#' # test data for cv_avg function
#' library(tidyverse)
#' set.seed(76)
#' x1 = rnorm(16, 6, 1)
#' x2 = rnorm(16, 6, 1)
#' test_dat <- data.frame(stdy = rep(c(1,2,3,4), each = 4),x1 = x1,sd1 = exp(log(x1)*1.5 + rnorm(16, 0,
#' sd = 0.10)),n1 = rpois(16, 15),x2 = x2,sd2 = exp(log(x2)*1.5 + rnorm(16, 0, sd = 0.10)),n2 = rpois(16, 15))
#' rm(list = c("x1", "x2"))
#' # # Now generate some missing data
#' t2 <- gen.miss(test_dat, "sd1", "sd2", 6)
#' t2_cv <- cv_avg(x = x1, sd = sd1, n = n1, stdy, data =  t2, sub_b = TRUE)
#' t2_cv <- cv_avg(x2, sd2, n2, stdy, label = "2", data =  t2_cv)
#' # Check calculations are correct. All match what is expected
#' test <- t2_cv %>%  filter(stdy == "1")
#' # Within
#' sum(test$n1) #59
#' sum(test$n2) #56
#' weighted.mean((test$sd1 / test$x1)^2, test$n1, na.rm = T)
#' weighted.mean((test$sd2 / test$x2)^2, test$n2, na.rm = T)
#' # Between
#' wCV1 = unique(t2_cv$w_CV2_1)
#' w_nt1 = c(59,58,58,50)
#' weighted.mean(wCV1, w_nt1)
#' wCV2 = unique(t2_cv$w_CV2_2)
#' w_nt2 = c(56, 56, 72, 63)
#' weighted.mean(wCV2, w_nt2)
#' }

cv_avg <- function(x, sd, n, group, data, label = NULL, sub_b = TRUE){

  # Check if the name is specified or not. If not, then assign it the name of the mean, x, variable input in the function. https://stackoverflow.com/questions/60644445/converting-tidyeval-arguments-to-string
  if(is.null(label)){
    label <- purrr::map_chr(enquos(x), rlang::as_label)
  }

  # Calculate between study CV. Take weighted mean CV within study, and then take a weighted mean across studies of the within study CV. Weighted based on sample size and pooled sample size.
  b_grp_cv_data <- data                                             %>%
    dplyr::group_by({{group}})                            %>%
    dplyr::mutate(   w_CV2 = weighted.mean(na_if(({{sd}} / {{x}})^2, Inf), {{n}},
                                           na.rm = TRUE),
                     n_mean = mean({{n}}, na.rm = TRUE))   %>%
    dplyr::ungroup(.)                                     %>%
    dplyr::mutate(b_CV2 = weighted.mean(w_CV2, n_mean, na.rm = TRUE), .keep = "used")

  # Make sure that label of the calculated columns is distinct from any other columns
  names(b_grp_cv_data) <- paste0(names(b_grp_cv_data), "_", label)

  # Append these calculated columns back to the original data and return the full dataset.
  if(sub_b){
    b_grp_cv_data <- b_grp_cv_data %>% dplyr::select(grep("b_", names(b_grp_cv_data)))
    dat_new <- cbind(data, b_grp_cv_data)
  } else {
    dat_new <- cbind(data, b_grp_cv_data)
  }

  return(data.frame(dat_new))
}


#' @title get_est
#' @description Extracts estimates from rma.mv and rma model objects
#' @param model The rma.mv model object

get_est <- function(model){
  est <- coef(model)
  ci.lb <- model$ci.lb
  ci.ub <- model$ci.ub
  se <- model$se
  return(data.frame(Est. = est, SE=se, "95% LCI" = ci.lb, "95% UCI" = ci.ub, check.names = FALSE))
}

#' @title lnrr_laj
#' @description Calculates log response ratio based on Taylor expansion from Jajeunesse 2011
#' @param m1 Mean of treatment group 1
#' @param m2 Mean of treatment group 2
#' @param cv1_2 Coefficient of variation squared (CV^2) for treatment group 1
#' @param cv2_2 Coefficient of variation squared (CV^2) for treatment group 2
#' @param n1 Sample size for treatment group 1
#' @param n2 Sample size for treatment group 2
#'
lnrr_laj <- function(m1, m2, cv1_2, cv2_2, n1, n2){
  log(m1 / m2) + 0.5*((cv1_2 / n1) - (cv2_2 / n2))
}

#' @title v_lnrr_laj
#' @description Calculates the sampling variance for log response ratio based on second order Taylor expansion proposed by Lajeunesse 2011
#' @param cv1_2 Coefficient of variation squared (CV^2) for treatment group 1
#' @param cv2_2 Coefficient of variation squared (CV^2) for treatment group 2
#' @param n1 Sample size for treatment group 1
#' @param n2 Sample size for treatment group 2
v_lnrr_laj <- function(cv1_2, cv2_2, n1, n2){
  ((cv1_2) / n1) + ((cv2_2) / n2) +
    ((cv1_2)^2 / (2*n1)^2) + ((cv2_2)^2 / (2*n2)^2)
}


#' @title tree_checks
#' @description A function for checking species names in a phylogenetic tree and dataframe to ensure they match.
#' @param data The dataframe containing the mean, sd, n and grouping variables along with species name
#' @param tree A phylogenetic tree
#' @param species_name_col A character string specifying the column name describing the species in the dataframe
#' @param type Type of method to be implemented. 'checks' will compare teh data and tree and identify 1) how many species are in the tree and data, 2) which species are in the data but not tree and 3) which species are in the tree but not the data. 'prune' can be used to prune the phylogenetic tree down to the species that exist within the dataframe.

tree_checks <- function(data, tree, species_name_col, type = c("checks", "prune")){
  type = match.arg(type)
  # How many unique species exist in data and tree
  Numbers <- matrix(nrow = 2, ncol = 1)
  Numbers[1,1] <- length(unique(data[,species_name_col]))
  Numbers[2,1] <- length(tree$tip.label)
  rownames(Numbers)<- c("Species in data:", "Species in tree:")
  # Missing species or species not spelt correct
  species_list1= setdiff(sort(tree$tip.label), sort(unique(data[,species_name_col])))
  species_list2= setdiff(sort(unique(data[,species_name_col])), sort(tree$tip.label) )
  if(type == "checks"){
    return(list(SpeciesNumbers = data.frame(Numbers),
                Species_InTree_But_NotData=species_list1,
                Species_InData_But_NotTree=species_list2))
  }
  if(type == "prune"){
    if(length(species_list2) >=1) stop("Sorry, you can only prune a tree when you have no taxa existing in the data that are not in the tree")
    return(ape::drop.tip(tree, species_list1))
  }
}

#' @title gen.miss
#' @description Generates random missing data in two columns of a meta-analytic dataset
#' @param data The dataframe containing the mean, sd, n and grouping variables
#' @param missVar A character string specifying the first column one wishes to have missing data. Example, "sd1" column.
#' @param missCol2 A character string specifying the second column one wishes to have missing data. Example, "sd2" column.
#' @param n_miss The number of missing data (NA) one wishes to introduce to the missVar and missCol2
gen.miss <- function(data, missVar, missCol2, n_miss){
  data[sample(rownames(data), n_miss), missVar] <- NA
  data[is.na(data[,missVar]), missCol2] <- NA
  return(data)
}

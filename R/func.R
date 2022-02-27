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


gen.miss <- function(data, missVar, missCol2, n_miss){
  data[sample(rownames(data), n_miss), missVar] <- NA
  data[is.na(data[,missVar]), missCol2] <- NA
  return(data)
}

gen.miss2 <- function(data, missVar, n_miss){
  data[sample(rownames(data), n_miss, replace = TRUE), missVar] <- NA
  return(data)
}

cv_avg <- function(cv, n, group, data, name){
  # Calculate between study CV (or whatever). Take weighted mean CV within study, and then take a weighted mean across studies of the within study CV. Weighted based on sample size and pooled sample size.
    b_grp_cv_data <- data                                      %>%
                group_by({{group}})                            %>%
                mutate(   w_CV = weighted.mean({{cv}}, {{n}},
                                               na.rm = TRUE),
                       n_total = sum({{n}}, na.rm = TRUE))     %>% # Pooled N or mean N?
                ungroup(.)                                     %>%
                mutate(b_CV = weighted.mean(w_CV, n_total), .keep = "used")

  # Make sure that label of the calculated columns is distinct from any other columns
    names(b_grp_cv_data) <- paste0(names(b_grp_cv_data), "_", name)

  # Append these calculated columns back to the original data and return the full dataset.
    dat_new <- cbind(data, b_grp_cv_data)

    return(data.frame(dat_new))
}


# # test data for cv_avg function
# library(tidyverse)
# set.seed(76)
# x1 = rnorm(16, 6, 1)
# x2 = rnorm(16, 6, 1)
# test_dat <- data.frame(stdy = rep(c(1,2,3,4), each = 4),
#                           x1 = x1,
#                          sd1 = exp(log(x1)*1.5 + rnorm(16, 0, sd = 0.10)),
#                           n1 = rpois(16, 15),
#                           x2 = x2,
#                          sd2 = exp(log(x2)*1.5 + rnorm(16, 0, sd = 0.10)),
#                           n2 = rpois(16, 15))
# test_dat <- test_dat %>% mutate(cv1 = x1 / sd1,
#                                 cv2 = x2 / sd2)
# plot(sd ~ x, data = test_dat)
# plot(log(sd) ~ log(x), data = test_dat)
#
# # Now generate some missing data
# t2 <- gen.miss(test_dat, "sd1", "cv1", 6)
# t2 <- gen.miss(test_dat, "sd2", "cv2", 6)
#
# t2_cv <- cv_avg(cv1, n1, stdy, name = "1", data =  t2)
# t2_cv <- cv_avg(cv2, n2, stdy, name = "2", data =  t2_cv)
#
# # Check calculations are correct. All match what is expected
# test <- t2_cv %>%  filter(stdy == "1")
#
# # Within
# sum(test$n1) #59
# sum(test$n2) #56
# weighted.mean(test$cv1, test$n1, na.rm = T) #0.4183278
# weighted.mean(test$cv2, test$n2, na.rm = T) #0.4579275
#
# # Between
# wCV1 = unique(t2_cv$w_CV_1)
# w_nt1 = c(59,58,58,50)
# weighted.mean(wCV1, w_nt1) #0.4207684
#
# wCV2 = unique(t2_cv$w_CV_2)
# w_nt2 = c(56, 56, 72, 63)
# weighted.mean(wCV2, w_nt2) # 0.4000197


get_est <- function(model){
  est <- coef(model)
  ci.lb <- model$ci.lb
  ci.ub <- model$ci.ub
  return(data.frame(Est. = est, "95% LCI" = ci.lb, "95% UCI" = ci.ub, check.names = FALSE))
}


#' createMiss: Function for creating missingness in a data frames variables
#' @param data the data frame to create missing values for variables. Note the last column should be an individual ID column
#' @param p.miss the percentage of missing data for each variable.
#' @param mcar if TRUE missingness is completely random for each variable.
#` @param W if TRUE missingness is at the within individual level (i.e. individuals are missing a percentage of values).
#` @param w.miss percentage of missingness at the within individual level (i.e. percentage of missing values within each individual). Used when mcar = 'TRUE', this will over-write p.miss.

createMiss <- function(data, grp = "ID", missVar = c("x1", "x2", "x3"), col = 1:(ncol(data) - 1), p.miss = rep(0.2, ncol(data) - 1),  type = c("w", "all", "nrand")){

  totCol <- 1:ncol(data)
  cols <- totCol[-match(c(grp), colnames(data))]

  #if(length(cols) != length(p.miss)){
  #stop("The percentage of missingness is not defined for some columns OR you have entered too many missingness values")
  #}

  type = match.arg(type)

  if(type  == "all"){
    # Will ignore individuals, but create a set amount of missingness throughout the entire dataset. Can specify different missingness for each variable.
    n_miss <- p.miss*nrow(data)
    data <- gen.miss(data, grp, n_miss)
    return(data)
  }

  if(type == "w"){
    # Split data by ID
    dtmp <- split(data, data[,grp])

    #Generate number of rows to create missing data for each ind
    n_miss <- p.miss*(nrow(data) / length(unique(data[,grp])))

    # Create missing data within each individual. Restricted to all individuals having same amount of data at this point. But you can change the missingness per variable per individual.
    dat <- lapply(dtmp, function(x) gen.miss2(x, missVar, n_miss))
    return(plyr::ldply(dat, .id = NULL))
  }


  if(type == "nrand"){
    #pattern1 <- c(rep(TRUE, 3), rep(FALSE,5), rep(TRUE, 2))
    pattern1 <- c(rep(FALSE, 3), rep(TRUE,5), rep(FALSE, 2))
    pattern2 <- c(rep(TRUE, 3), rep(FALSE,5), rep(TRUE, 2))


    # Split data by ID
    # Sample IDs within treatments
    ids <- as.vector(unlist(as.vector(lapply(data %>% split(trt, drop = TRUE), function(x) sample(unique(x$ind), 5)))))

    # Split by ID
    dtmp <- split(data, data[,grp])

    # create missing data patterns
    for(i in names(dtmp)){
      if(i %in% ids){
        dtmp[[i]][pattern1, missVar] <- NA
      } else{
        dtmp[[i]][pattern2, missVar] <- NA
      }
    }

    return(plyr::ldply(dtmp, .id = NULL))
  }

}

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


cv_avg <- function(cv, n, group, data){
  b_grp_cv_data <- data                                      %>%
              group_by({{group}})                            %>%
              mutate(   w_CV = weighted.mean({{cv}}, {{n}},
                                             na.rm = TRUE),
                     n_total = sum({{n}}, na.rm = TRUE))     %>%
              ungroup(.)                                     %>%
              mutate(b_CV = weighted.mean(w_CV, n_total))

  return(data.frame(b_grp_cv_data))

}


# test data for cv_avg function
library(tidyverse)
set.seed(76)
x1 = rnorm(16, 6, 1)
test_dat <- data.frame(stdy = rep(c(1,2,3,4), each = 4),
                          x = x1,
                         sd = exp(log(x1)*1.5 + rnorm(16, 0, sd = 0.10)),
                          n = rpois(16, 15))
test_dat <- test_dat %>% mutate(cv = x / sd)
plot(sd ~ x, data = test_dat)
plot(log(sd) ~ log(x), data = test_dat)

t <- cv_avg(cv, n, stdy, data =  test_dat)

# Now generate some missing data
t2 <- gen.miss(test_dat, "sd", "cv", 6)

t2_cv <- cv_avg(cv, n, stdy, data =  t2)

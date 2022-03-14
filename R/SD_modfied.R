# test for modified SD

#install.packages("pacman")
pacman::p_load(tidyverse, metafor, here, osfr, ape, phytools)

# additional

lnrr_laj2 <- function(m1, m2, sd1, sd2, n1, n2){
  sd1c <- sd1*exp(1/(2*n1-2))
  sd2c <- sd2*exp(1/(2*n2-2))
  lnrr_laj2<- log(m1 / m2) + 0.5*((sd1c^2 /(n1*m1^2)) - (sd2c^2 /(n2*m2^2)))
  lnrr_laj2
}

v_lnrr_laj2 <- function(m1, m2, sd1, sd2, n1, n2){
  sd1c <- sd1*exp(1/(2*n1-2))
  sd2c <- sd2*exp(1/(2*n2-2))
  v_lnrr_laj2 <- (sd1c^2 /(n1*m1^2)) + (sd2c^2 /(n2*m2^2)) +
    (sd1c^4 /((2*n1^2)*(m1^4))) + (sd1c^4 /((2*n1^2)*(m1^4)))
  v_lnrr_laj2
}

cv_avg2 <- function(x, sd, n, group, data, label = NULL, sub_b = TRUE){

  # Check if the name is specified or not. If not, then assign it the name of the mean, x, variable input in the function. https://stackoverflow.com/questions/60644445/converting-tidyeval-arguments-to-string
  if(is.null(label)){
    label <- purrr::map_chr(enquos(x), rlang::as_label)
  }

  # Calculate between study CV. Take weighted mean CV within study, and then take a weighted mean across studies of the within study CV. Weighted based on sample size and pooled sample size.
  b_grp_cv_data <- data                                             %>%
    dplyr::group_by({{group}})                            %>%
    dplyr::mutate(sdc = {{sd}}*exp(1/(2*({{n}}-1))),
      w_CV2 = weighted.mean(na_if((sdc/{{x}})^2, Inf), {{n}},
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



# Useful functions for calculating CV^2 within and between studies
osfr::osf_retrieve_file("https://osf.io/sqr4w/") %>% osfr::osf_download(conflicts = "overwrite")
source("./func.R")

# Download the data file from OSF
osfr::osf_retrieve_file("https://osf.io/evysw/") %>% osfr::osf_download(conflicts = "overwrite")

# Download the phylogeny from OSF
osfr::osf_retrieve_file("https://osf.io/t5kh4/") %>% osfr::osf_download(conflicts = "overwrite")

# Load the data file and tree file
data1 <- read.csv(here::here("example1.csv"))
tree <- read.tree(here::here("phylo_tree"))
phylo <- vcv(tree, corr = TRUE)

#

data1a <- data1

# Calculate the average between study CV, which will replace missing values.
data1 <- cv_avg(x = Control_mean, sd = Control_standard_deviation,
                n = Control_sample_size, group = Author, label = "1",
                data = data1)
data1 <- cv_avg(x = Experimental_mean, sd = Experimental_standard_deviation,
                n = Experimental_sample_size, group = Author,
                label = "2", data = data1)


data1a <- cv_avg2(x = Control_mean, sd = Control_standard_deviation,
                n = Control_sample_size, group = Author, label = "1",
                data = data1a)
data1a <- cv_avg2(x = Experimental_mean, sd = Experimental_standard_deviation,
                n = Experimental_sample_size, group = Author,
                label = "2", data = data1a)


# Use weighted mean CV in replacement for where CV's are missing. Otherwise, calculate CV^2 of data that is known.
data1 <- data1 %>%
  mutate(cv2_cont_new = if_else(is.na(Control_standard_deviation),      b_CV2_1, cv_Control^2),
         cv2_expt_new = if_else(is.na(Experimental_standard_deviation), b_CV2_2, cv_Experimental^2))

data1a <- data1a %>%
  mutate(cv2_cont_new = if_else(is.na(Control_standard_deviation),      b_CV2_1, cv_Control^2),
         cv2_expt_new = if_else(is.na(Experimental_standard_deviation), b_CV2_2, cv_Experimental^2))

# Now calculate new yi and vi, called lnrr_laj & v_lnrr_laj, respectively. This uses either the between individual CV^2 when missing or normal CV^2 when not missing.
data1 <- data1 %>%
  mutate(lnrr_laj = lnrr_laj(m1 = Control_mean, m2 = Experimental_mean,
                             cv1_2 = cv2_cont_new, cv2_2 = cv2_expt_new,
                             n1= Control_sample_size, n2 = Experimental_sample_size),
         v_lnrr_laj = v_lnrr_laj(cv1_2 = cv2_cont_new, n1= Control_sample_size,
                                 cv2_2 = cv2_expt_new, n2 = Experimental_sample_size))

data1a <- data1a %>%
  mutate(lnrr_laj2 = if_else(is.na(Control_standard_deviation),
                            lnrr_laj(m1 = Control_mean, m2 = Experimental_mean,
                             cv1_2 = cv2_cont_new, cv2_2 = cv2_expt_new,
                             n1= Control_sample_size, n2 = Experimental_sample_size),
                            lnrr_laj2(m1 = Control_mean, m2 = Experimental_mean,
                                     sd1 = Control_standard_deviation, sd2 = Experimental_standard_deviation,
                                     n1= Control_sample_size, n2 = Experimental_sample_size)),
         v_lnrr_laj2 = if_else(is.na(Control_standard_deviation),
                               v_lnrr_laj(cv1_2 = cv2_cont_new, n1= Control_sample_size,
                                 cv2_2 = cv2_expt_new, n2 = Experimental_sample_size),
                               v_lnrr_laj2(m1 = Control_mean, m2 = Experimental_mean,
                                           sd1 = Control_standard_deviation, sd2 = Experimental_standard_deviation,
                                          n1= Control_sample_size, n2 = Experimental_sample_size)
                               )
         )


plot(log(data1$v_lnrr_laj), log(data1a$v_lnrr_laj2))
cor(log(data1$v_lnrr_laj), log(data1a$v_lnrr_laj2))
plot(data1$lnrr_laj, data1a$lnrr_laj2)
cor.test(data1$lnrr_laj, data1a$lnrr_laj2, method = "spearman")
cor.test(data1$lnrr_laj, data1a$lnrr_laj2, method = "pearson")
hist(log(data1$v_lnrr_laj))
hist(log(data1a$v_lnrr_laj2))

# 1A
method_1A_bird <- rma.mv(lnrr_laj ~ 1, V = v_lnrr_laj,
                         random=list(~1|Group, ~1|Year, ~1|Focal_insect, ~1|obs),
                         R = list(Focal_insect = phylo), data = data1)

method_1A_bird_res <- get_est(method_1A_bird)

# 1B

# Now calculate new v_lnrr_laj_1B using between study CV^2 for all observations.
data1 <- data1 %>%
  mutate(v_lnrr_laj_1B = v_lnrr_laj(cv1 = b_CV2_1, n1= Control_sample_size,
                                    cv2 = b_CV2_2, n2 = Experimental_sample_size))

# Fit model with new sampling variance and point estimate
method_1B_bird <- rma.mv(lnrr_laj ~ 1, V = v_lnrr_laj_1B,
                         random=list(~1|Group, ~1|Year, ~1|Focal_insect, ~1|obs),
                         R = list(Focal_insect = phylo), data = data1)

method_1B_bird_res <- get_est(method_1B_bird)

# 1C
method_1C_bird <- rma.mv(lnrr_laj2 ~ 1, V = v_lnrr_laj2,
                         random=list(~1|Group, ~1|Year, ~1|Focal_insect, ~1|obs),
                         R = list(Focal_insect = phylo), data = data1)
method_1C_bird_res <- get_est(method_1C_bird)

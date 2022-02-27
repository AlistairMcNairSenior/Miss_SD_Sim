############################
# Load libraries
###########################
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

tree_checks(data = dat, tree = tree, species_name_col = "Focal_insect", type = "checks")

##################
# Checks with authors code
##################
#set up for phylogeny model

#create copy of data file to use with tree
a1 <- dat

#align focal species with branches of phylogenetic tree
tree1<-drop.tip(tree,tree$tip.label[which(tree$tip.label%in%a1$Focal_insect==FALSE)])
check.species<-function(x) {any(x==tree1$tip.label)}
a1 <- a1[sapply(a1[,"Focal_insect"],check.species),]

# Quick check to make sure that pruning has been done correctly.
tree_checks(data = a1, tree = tree1, species_name_col = "Focal_insect", type = "checks")

INtree <- inverseA(tree1)  # Metafor takes a correlation matrix
INtree2 <- vcv(tree1, corr = TRUE)  # Metafor takes a correlation matrix

#phylogeny model with fixed effects
    rerun=FALSE
    if(rerun){
      #priors
      prior3<-list(R=list(V=1,nu=0.002),
                   G=list(G1=list(V=1,nu=0.002),
                          G2=list(V=1,nu=0.002),
                          G3=list(V=1,nu=0.002)))
      prior6<-list(R=list(V=1,nu=0.002),
                   G=list(G1=list(V=1,nu=0.002),
                          G2=list(V=1,nu=0.002),
                          G3=list(V=1,nu=0.002),
                          G4=list(V=1,nu=0.002),
                          G5=list(V=1,nu=0.002),
                          G6=list(V=1,nu=0.002)))

      phylogenymodel<-MCMCglmm(yi_g ~ Spatial_separation + temporal_separation + F_Guild + C_Guild + F_Guild:C_Guild,
                             random=~Group+Year+Focal_insect, data=a1, mev=a1$vi, ginverse=list(Focal_insect=INtree$Ainv),
                             prior=prior3, nitt=1000000, thin=100, burnin=500000)

    saveRDS(phylogenymodel, "./example/worked2/phylogenymodel")} else {phylogenymodel <- readRDS("./example/worked2/phylogenymodel")}


    ## Metafor models. First, check it matches autors code.
    phylogenymodel_mv <- rma.mv(yi_g ~ Spatial_separation + temporal_separation + F_Guild + C_Guild + F_Guild:C_Guild, V = vi, random=list(~1|Group, ~1|Year, ~1|Focal_insect, ~1|obs), R = list(Focal_insect = INtree2), data = a1)


#######################
    # Example 2
#######################
# Intercept only model. for Abundance fitness ONLY. The authors state: "We first estimated the average impact of insect interaction on each fitness component by fitting linear models in which the intercept was the only fixed effect." Table 1 provides the fitness components. Abundance has the most data.
    a2 <- dat %>% filter(Fitness_component == "Abundance") # Doesn't quite match Table 1, but does match figure 1, and even exclusion N = 311, which matches.

    # Prune tree
    #align focal species with branches of phylogenetic tree
    tree2_meta<-drop.tip(tree,
                         tree$tip.label[which(tree$tip.label%in%a2$Focal_insect==FALSE)])
    check.species<-function(x) {any(x==tree2_meta$tip.label)}
    a2 <- a2[sapply(a2[,"Focal_insect"],check.species),]

    tree_checks(data = a2, tree = tree2_meta, species_name_col = "Focal_insect", type = "checks")

    phylo <- vcv(tree2_meta, corr = TRUE)


#####################
# Complete data model
#####################

    complete_mv <- rma.mv(yi_g ~ 1, V = vi,
                          random=list(~1|Group, ~1|Year, ~1|Focal_insect, ~1|obs),
                          R = list(Focal_insect = phylo), data = a2)
    complete_mv_res <- get_est(complete_mv)
#####################
# Generate missing data & complete case analysis
#####################
    # Generate missing data in SD's at the effect size level. If you're missing one SD then your missing for control and experimental
    a2_missSD <- gen.miss(a2, missVar = "Experimental_standard_deviation",
                           missCol2 = "Control_standard_deviation",
                           n_miss = 0.2*nrow(a2))

    # Now, assume you needto exclude data with missing SD because you can't calculate effect size and sampling variance
    complete_case_MV <- na.omit(a2_missSD)

    # Prune tree.
    tree3_meta<-drop.tip(tree,
                         tree$tip.label[which(tree$tip.label%in%complete_case_MV$Focal_insect==FALSE)])
    check.species<-function(x) {any(x==tree3_meta$tip.label)}
    complete_case_MV <- complete_case_MV[sapply(complete_case_MV[,"Focal_insect"],check.species),]

    #Check that it was done right
    tree_checks(data = complete_case_MV, tree = tree3_meta, species_name_col = "Focal_insect", type = "checks")

    # Make new VCV matrix
    phylo2 <- vcv(tree3_meta, corr = TRUE)

    # Fit complete case analysis. Note that data is currently missing at random.
    complete_case_mv <- rma.mv(yi_g ~ 1, V = vi,
                          random=list(~1|Group, ~1|Year, ~1|Focal_insect, ~1|obs),
                          R = list(Focal_insect = phylo2), data = complete_case_MV)

    complete_case_mv_res <- get_est(complete_case_mv)

#run libraries
pacman::p_load(ape, metafor, MCMCglmm, tidyverse, readxl)

#load data
dat  <- read_excel("./example/worked2/ele13245-sup-0003-metas3.xlsx")

#remove rows with missing data
dat <- na.omit(dat)

#load phylogenetic tree
tree<- read.tree("./example/worked2/ele13245-sup-0008-phylogenys8.tre")

#priors
prior3<-list(R=list(V=1,nu=0.002),G=list(G1=list(V=1,nu=0.002),G2=list(V=1,nu=0.002),G3=list(V=1,nu=0.002)))
prior6<-list(R=list(V=1,nu=0.002),G=list(G1=list(V=1,nu=0.002),G2=list(V=1,nu=0.002),G3=list(V=1,nu=0.002),G4=list(V=1,nu=0.002),G5=list(V=1,nu=0.002),G6=list(V=1,nu=0.002)))

#taxonomy model with fixed effects
taxonomymodel <- MCMCglmm(yi~ density + species + spatial_separation + temporal_separation + F_Guild + C_Guild + F_Guild:C_Guild,
                            random=~Group+Year+Focal_insect+F_Order+F_Family+F_Genus, mev=dat$vi, data=dat,prior=prior6,
                            family="gaussian", nitt=1000000, thin=100, burnin=100000)


#set up for phylogeny model

#create copy of data file to use with tree
a1 <- dat

#align focal species with branches of phylogenetic tree
tree1<-drop.tip(tree,tree$tip.label[which(tree$tip.label%in%a1$`Focal insect`==FALSE)])
check.species<-function(x) {any(x==tree1$tip.label)}
a1 <- a1[sapply(a1[,"Focal insect"],check.species),]
INtree <- inverseA(tree1,nodes="TIPS")

#phylogeny model with fixed effects
phylogenymodel<-MCMCglmm(yi ~ density + species + spatial_separation + temporal_separation + F_Guild + C_Guild + F_Guild:C_Guild,
                         random=~Group+Year+Focal_insect, data=a1, mev=a1$vi, ginverse=list(Focal_insect=INtree$Ainv),
                         prior=prior3, nitt=1000000, thin=100, burnin=500000)



#create subset to only model abundance data
dat <- subset(dat, Measurement_type == "Abundance")

#taxonomy model with fixed effects
taxonomymodelabundance <- MCMCglmm(yi~ density + species + spatial_separation + temporal_separation + F_Guild + C_Guild + F_Guild:C_Guild,
                          random=~Group+Year+Focal_insect+F_Order+F_Family+F_Genus, mev=dat$vi, data=dat,prior=prior6,
                          family="gaussian", nitt=1000000, thin=100, burnin=100000)


#set up for phylogeny model

#create copy of data file to use with tree
a1 <- dat

#align focal species with branches of phylogenetic tree
tree1<-drop.tip(tree,tree$tip.label[which(tree$tip.label%in%a1$Focal_insect==FALSE)])
check.species<-function(x) {any(x==tree1$tip.label)}
a1 <- a1[sapply(a1[,"Focal_insect"],check.species),]
INtree <- inverseA(tree1,nodes="TIPS")

#phylogeny model with fixed effects
phylogenymodelabundance<-MCMCglmm(yi ~ density + species + spatial_separation + temporal_separation + F_Guild + C_Guild + F_Guild:C_Guild,
                         random=~Group+Year+Focal_insect, data=a1, mev=a1$vi, ginverse=list(Focal_insect=INtree$Ainv),
                         prior=prior3, nitt=1000000, thin=100, burnin=500000)






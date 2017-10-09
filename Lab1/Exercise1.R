#!/usr/bin/Rscript
#Author: "Emma Fox (eafox@ucla.edu)"

# SET WD #
setwd("UCLA/Evo/Lab1/")
#remove if using in other location

# PACKAGES #
library(phytools)

# LOAD DATA #
darter.tree<-read.tree("etheostoma_percina_chrono.tre")

# PLOT TREE #
plotTree(darter.tree,ftype="i",fsize=0.4,type="fan",lwd=1)

# LTT PLOT #
obj<-ltt(darter.tree,log.lineages=FALSE)
obj

# EXAMINE DATA SET #
darter.tree
#bifurcated tree would have 200 internal nodes (n_tips-1)
is.binary(darter.tree)
#another check of whether fully bifurcating

# FULLY BIFURCATE TREE #
darter.tree<-multi2di(darter.tree)
darter.tree
is.binary(darter.tree)
#"randomly resolved the internal nodes that were multifurcating"

# STORE LTT #
obj<-ltt(darter.tree,plot = FALSE)
obj

# PLOT FULLY BIFURCATED DATA SET #
plot(obj,log.lineages=FALSE,main="LTT plot for darters")

# OVERLAY TREE ONTO LTT PLOT #
plot(obj,show.tree=TRUE,log.lineages=FALSE,main="LTT plot for darters")

# PLOT LOG OF LINEAGE NUMBER AGAINST TIME #
plot(obj,log.lineages=FALSE,log="y",main="LTT plot for darters",
     ylim=c(2,Ntip(darter.tree)))

# ADD LINE SHOWING PREDICTION UNDER PURE BIRTH MODEL #
h<-max(nodeHeights(darter.tree))
x<-seq(0,h,by=h/100)
b<-(log(Ntip(darter.tree))-log(2))/h
lines(x,2*exp(b*x),col="red",lty="dashed",lwd=2)

# ADD SIMULATED LINEAGE NUMBERS THROUGH TIME #
trees<-pbtree(b=b,n=Ntip(darter.tree),t=h,nsim=100,method="direct",
              quiet=TRUE)
obj<-ltt(trees,plot=FALSE)
plot(obj,col="grey",main="LTT of darters compared to simulated LTTs")
lines(c(0,h),log(c(2,Ntip(darter.tree))),lty="dashed",lwd=2,col="red")

# ADD ORIGINAL DATA #
ltt(darter.tree,add=TRUE,lwd=2)

# 95% CONFIDENCE INTERVAL FOR LTT FROM A SET OF TREES #
ltt95(trees,log=TRUE)
title(main="LTT of darters compared to simulated LTTs")
ltt(darter.tree,add=TRUE,log.lineages=FALSE,col="red",lwd=2)







#!/usr/bin/Rscript
#Author: "Emma Fox (eafox@ucla.edu)"

# SET WD #
#setwd("UCLA/Evo/Lab1/")
#remove if using in other location

# PACKAGES #
library(phytools)

# LOAD DATA #
darter.tree<-read.tree("etheostoma_percina_chrono.tre")

# FULLY BIFURCATE TREE #
darter.tree<-multi2di(darter.tree)
darter.tree
is.binary(darter.tree)

# GENERATE SIMULATED TREES #
h<-max(nodeHeights(darter.tree))
x<-seq(0,h,by=h/100)
b<-(log(Ntip(darter.tree))-log(2))/h
trees<-pbtree(b=b,n=Ntip(darter.tree),t=h,nsim=100,method="direct",
              quiet=TRUE)

# PLOT HISTOGRAM OF SIMULATED GAMMA VALUES #
g<-sapply(trees,function(x) ltt(x,plot=FALSE)$gamma)
hist(g,main=expression(paste("Distribution of ",gamma," from simulation")))

# RETURN SUMMARY STATS #
mean(g)
var(g)

# TESTING HYPOTHESES #
obj<-ltt(darter.tree,plot=FALSE)
print(obj)

# SIMULATING A COALESCENT TREE #
coal.tree<-rcoal(n=100)
plotTree(coal.tree,ftype="off")
coal.obj<-ltt(coal.tree,plot=FALSE)
print(coal.obj)

# PRINT COAL LTT PLOT #
obj<-ltt(coal.tree,log.lineages=FALSE,log="y")

# SAVE GAMMA VALUE #
darter.gamma <- obj$gamma
darter.gamma

# COMPARING GAMMA BETWEEN PURE BIRTH AND COALESCENT SIMULATED TREES #
trees<-pbtree(n=100,nsim=200,scale=max(nodeHeights(coal.tree)))
ltt95(trees,log=TRUE)
title(main="Simulated coalescent trees compared to pure-birth LTTs")
ltt(coal.tree,add=TRUE,log.lineages=FALSE,col="red",lwd=2,lty="dashed")
#coalescent is exponential even when logged

#######################
# INCOMPLETE SAMPLING #
#######################
#Authors were missing 15 species
#Going to see what affect that has by generating 216 species and subsampling down to 201
#and examining each gamma. If ours is extreme, can conclude not caused by incomplete sampling

# SIMULATE PURE BIRTH TREE #
library(geiger)
age <- 25.91862
richness <- 216
darterbirth =  (log(richness) - log(2))/age
darterbirth

# SIMULATE INCOMPLETE SAMPLING #
richness <- 216
missing <- 15
num_simulations<-200 
g1_null<-numeric(num_simulations) #g1_null will hold the simulated gamma values
for(i in 1:num_simulations) {
  sim.bdtree(darterbirth, d=0, stop = "taxa", n=richness)->sim_tree 
  drop.random(sim_tree, missing)->prune # prune down to the # of taxa in the phylogeny
  gammaStat(prune)->g1_null[i]
}

# PLOT HISTOGRAM OF NULL DISTRIBUTION #
hist(g1_null)

# ADD ARROW TO SHOW OBSERVED GAMMA # 
arrows(darter.gamma, 40, darter.gamma, 0, col="red", lwd=2) 

# COMPARE OBSERVED GAMMA TO SIMULATED VALUES #
smallerNull<-g1_null<=darter.gamma
count<-sum(smallerNull)

# USE TO GENERATE P VALUE #
mccr_pval<-(count+1)/(num_simulations+1)
mccr_pval
#definitely just due to subsampling
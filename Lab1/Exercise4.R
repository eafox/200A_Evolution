#!/usr/bin/Rscript
#Author: "Emma Fox (eafox@ucla.edu)"

# SET WD #
#setwd("UCLA/Evo/Lab1/")
#remove if using in other location

# PACKAGES #
library(phytools)
library(geiger)

##############
# QUESTION 1 #
##############

# Load data
snake.tree<-read.tree("homalops.phy")

# Check tree is fully bifurcated
snake.tree
# Looks to be! One more check...
is.binary(snake.tree)

# Run ltt
snake.obj<-ltt(snake.tree,plot = FALSE)
snake.obj
snake.gamma <- snake.obj$gamma
snake.gamma

# GAMMA: -3.2411, p-val: 0.0012

##############
# QUESTION 2 #
##############
# The gamma value is conclusively negative which means there are more nodes/divergences earlier in the tree
# This would be reason for further study to determine if an adaptive radiation event has occurred

##############
# QUESTION 3 #
##############
# Determine if observed gamma is due to incomplete sampling

# Parameters for pure birth model
age <- 22
richness <- 34
snakebirth =  (log(richness) - log(2))/age
snakebirth

# Parameters for tree simulation with incomplete sampling
richness <- 216
missing <- 15
num_simulations<-200 

# For loop to simulate populations 
# COPIED FROM MATERIAL PROVIDED BY PROF. MIKE ALFARO
g1_null<-numeric(num_simulations) #g1_null will hold the simulated gamma values
for(i in 1:num_simulations) {
  sim.bdtree(snakebirth, d=0, stop = "taxa", n=richness)->sim_tree 
  drop.random(sim_tree, missing)->prune # prune down to the # of taxa in the phylogeny
  gammaStat(prune)->g1_null[i]
}

# Plot histogram of values and point to observed value
hist(g1_null)
arrows(snake.gamma, 40, snake.gamma, 0, col="red", lwd=2) 
#seems to be well into the tails. Looks possibily significant

# Generate p-value of observed gamma
smallerNull<-g1_null<=snake.gamma
count<-sum(smallerNull)
mccr_pval<-(count+1)/(num_simulations+1)
mccr_pval
# p-val: 0.00498

# On the basis of the MCCR test, I conclude the rate of speciation slowing as time increases

##############
# QUESTION 4 #
##############

# Fit birth-death model
fitbd <- birthdeath(snake.tree)
fitbd
#Print birth and death rates
bd(fitbd)

# Birth Rate: 0.0684 
# Death Rate: 0.0000

##############
# QUESTION 5 #
##############
#(1) Clade description and source
# Number of Tips:
# Number of Species in Clades:
# Citation:

# Perissodactyla (odd-toed ungulates and cetaceans) data set downloaded from:
# http://10ktrees.nunn-lab.org/Perissodactyla/dataset.html

# Load tree
prac.tree<-read.tree("")
# Check tree is fully bifurcated
prac.tree
# Looks to be! One more check...
is.binary(prac.tree)


#(2) Fit birth death model and return b and d
# Fit birth-death model
fitbd <- birthdeath(prac.tree)
fitbd
#Print birth and death rates
bd(fitbd)
# Birth Rate:  
# Death Rate: 

#(3) Perform MCCR test and analyse
# Run ltt and get observed gamma
prac.obj<-ltt(prac.tree,plot = FALSE)
prac.obj
prac.gamma <- prac.obj$gamma
prac.gamma
# Simulate incomplete sampling
age <- 22
richness <- 34
pracbirth =  (log(richness) - log(2))/age
pracbirth
richness <- 216
missing <- 15
num_simulations<-200 
# For loop to simulate populations 
# COPIED FROM MATERIAL PROVIDED BY PROF. MIKE ALFARO
g1_null<-numeric(num_simulations) #g1_null will hold the simulated gamma values
for(i in 1:num_simulations) {
  sim.bdtree(pracbirth, d=0, stop = "taxa", n=richness)->sim_tree 
  drop.random(sim_tree, missing)->prune # prune down to the # of taxa in the phylogeny
  gammaStat(prune)->g1_null[i]
}
# Plot histogram of values and point to observed value
hist(g1_null)
arrows(prac.gamma, 40, prac.gamma, 0, col="red", lwd=2) 
#seems to be well into the tails. Looks possibily significant
# Generate p-value of observed gamma
smallerNull<-g1_null<=prac.gamma
count<-sum(smallerNull)
mccr_pval<-(count+1)/(num_simulations+1)
mccr_pval
# p-val:

# Interpratation 












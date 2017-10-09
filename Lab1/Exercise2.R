#!/usr/bin/Rscript
#Author: "Emma Fox (eafox@ucla.edu)"

# SET WD #
setwd("UCLA/Evo/Lab1/")
#remove if using in other location

# PACKAGES #
library(phytools)

# LOAD DATA #
darter.tree<-read.tree("etheostoma_percina_chrono.tre")

# FULLY BIFURCATE TREE #
darter.tree<-multi2di(darter.tree)
darter.tree
is.binary(darter.tree)

# FIT BIRTH-DEATH MODEL #
fitbd <- birthdeath(darter.tree)
fitbd

# RETURN BIRTH AND DEATH RATES #
bd(fitbd)

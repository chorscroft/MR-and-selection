require(systemfit)

################################################################
## Basic simulation, 0.05(0.01)
################################################################
setwd("~/Bristol/Simulations/basic")

dat<-read.table("output.txt",header = TRUE)
sf<-systemfit(fit~pheno,inst = ~geno,method = "2SLS",data=dat)
summary(sf)
require(systemfit)

################################################################
## Basic simulation, 0.05(0.01)
################################################################
setwd("~/Bristol/Simulations/children")

dat<-read.table("output.txt",header = TRUE)
par<-read.table("parent.txt",header = TRUE)

dat$children<-sapply(dat$index,function(x) sum(par$parent1==x | par$parent2==x))

sf<-systemfit(fit~pheno,inst = ~geno,method = "2SLS",data=dat)
summary(sf)

sf_children<-systemfit(children~pheno,inst = ~geno,method = "2SLS",data=dat)
summary(sf_children)
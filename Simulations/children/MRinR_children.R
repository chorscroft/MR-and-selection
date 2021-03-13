require(systemfit)

################################################################
## Basic simulation, 0.05
################################################################
setwd("~/Documents/MR-and-selection/Simulations/children")

dat<-read.table("dat.txt",header = TRUE)
table(dat$geno)
tapply(dat$pheno,dat$geno,mean)
tapply(dat$fit,dat$geno,mean)
tapply(dat$children,dat$geno,mean)

sf<-systemfit(fit~pheno,inst = ~geno,method = "2SLS",data=dat)
summary(sf)

sf_children<-systemfit(children~pheno,inst = ~geno,method = "2SLS",data=dat)
summary(sf_children)

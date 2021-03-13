require(systemfit)

################################################################
## Basic simulation, 0.05(0.01)
################################################################
setwd("~/Documents/MR-and-selection/Simulations/basic")

dat<-read.table("output.txt",header = TRUE)
table(dat$geno)
tapply(dat$pheno,dat$geno,mean)
tapply(dat$fit,dat$geno,mean)
plot(dat$pheno,dat$fit)
summary(lm(dat$fit~dat$pheno))

sf<-systemfit(fit~pheno,inst = ~geno,method = "2SLS",data=dat)
summary(sf)

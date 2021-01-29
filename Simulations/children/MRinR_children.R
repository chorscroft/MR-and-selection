require(systemfit)

################################################################
## Basic simulation, 0.05(0.01)
################################################################
setwd("~/Bristol/Simulations/children")

dat<-read.table("output.txt",header = TRUE)
par<-read.table("parent.txt",header = TRUE)

dat$children<-sapply(dat$index,function(x) sum(par$parent1==x | par$parent2==x))
dat$childrenOE<-dat$children/2

barplot(table(dat$children),space = 0)
lines(0:10+0.5,dpois(0:10,2)*10000,col="blue")

sf<-systemfit(fit~pheno,inst = ~geno,method = "2SLS",data=dat)
summary(sf)

sf_children<-systemfit(childrenOE~pheno,inst = ~geno,method = "2SLS",data=dat)
summary(sf_children)
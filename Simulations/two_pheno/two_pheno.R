require(systemfit)

################################################################
## Two phenotype simulation, 0.05, -0.03, additive
################################################################
setwd("~/Documents/MR-and-selection/Simulations/two_pheno/")

dat<-read.table("output.txt",header = TRUE)
table(dat$geno)
tapply(dat$pheno1,dat$geno,mean)
tapply(dat$pheno2,dat$geno,mean)
tapply(dat$fit,dat$geno,mean)

sf1<-systemfit(fit~pheno1+pheno2,inst = ~geno,method = "2SLS",data=dat)
summary(sf1)
sf2<-systemfit(fit~pheno1,inst = ~geno,method = "2SLS",data=dat)
summary(sf2)
sf3<-systemfit(fit~pheno2,inst = ~geno,method = "2SLS",data=dat)
summary(sf3)

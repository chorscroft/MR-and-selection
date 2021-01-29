require(systemfit)

################################################################
## Two phenotype simulation, 0.05(0.01), -0.03(0.01), additive
################################################################
setwd("~/Bristol/Simulations/two_pheno/")

dat2<-read.table("output.txt",header = TRUE)
sf1<-systemfit(fit~pheno1+pheno2,inst = ~geno,method = "2SLS",data=dat2)
summary(sf1)
sf2<-systemfit(fit~pheno2,inst = ~geno,method = "2SLS",data=dat2)
summary(sf2)
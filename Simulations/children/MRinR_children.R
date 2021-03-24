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

sf<-systemfit(log(fit)~pheno,inst = ~geno,method = "2SLS",data=dat)
summary(sf)

mod1<-summary(lm(dat$pheno~dat$geno))
bgx<-mod1$coefficients[2,1]
segx<-mod1$coefficients[2,2]
mod2<-summary(glm(dat$children~dat$geno,family=poisson("log")))
bgy<-mod2$coefficients[2,1]
segy<-mod2$coefficients[2,2]
mr_wald_ratio(bgx,bgy,segx,segy)

require(systemfit)
require(jsonlite)
config <- read_json("~/config.json")
mr_sel_path <- config$mr_sel_path

################################################################
## Basic simulation, 0.05(0.01)
################################################################

dat<-read.table(paste0(mr_sel_path,"/Simulations/basic/output.txt"),header = TRUE)
table(dat$geno)
tapply(dat$pheno,dat$geno,mean)
tapply(dat$fit,dat$geno,mean)
plot(dat$pheno,dat$fit)
summary(lm(dat$fit~dat$pheno))

sf<-systemfit(fit~pheno,inst = ~geno,method = "2SLS",data=dat)
summary(sf)

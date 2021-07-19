require(systemfit)
require(jsonlite)
config <- read_json("~/config.json")
mr_sel_path <- config$mr_sel_path

################################################################
## Two phenotype simulation, 0.05, -0.03, additive
################################################################

dat<-read.table(paste0(mr_sel_path,"/Simulations/two_pheno/output.txt"),header = TRUE)
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

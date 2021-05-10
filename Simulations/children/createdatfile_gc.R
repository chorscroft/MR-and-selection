#!/usr/bin/env Rscript
config <- read_json("~/config.json")
mr_sel_path <- config$mr_sel_path
workingdir<-paste0(mr_sel_path,"/Simulations/children/multigen/selStop_sim_1")
datfile<-"dat_0.txt"
datcfile<-"dat_1.txt"
datgcfile<-"dat_gc.txt"


## function to create dat_gc.txt file from two dat.txt files in the given directory
create_dat_file<-function(workingdir,datfile,datcfile,datgcfile="dat_gc.txt"){
  dat<-read.table(paste0(workingdir,"/",datfile),header = TRUE)
  dat_c<-read.table(paste0(workingdir,"/",datcfile),header = TRUE)
  dat$grandchildren<-0
  for (i in 1:nrow(dat)){
    dat$grandchildren[i]<-sum(dat_c$children[dat_c$parent1==dat$index[i]])+sum(dat_c$children[dat_c$parent2==dat$index[i]])
  }
  write.table(dat,paste0(workingdir,"/",datfile),quote=FALSE,row.names = FALSE)
}
## runs the function given a directory
args = as.list(commandArgs(TRUE))
names(args)<-c("workingdir","outfile","parfile","datfile")[1:length(args)]
do.call(create_dat_file,args)
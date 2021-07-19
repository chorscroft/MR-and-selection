#!/usr/bin/env Rscript

## function to create dat_gc.txt file from two dat.txt files in the given directory
create_dat_gc_file<-function(workingdir,datfile,datcfile,datgcfile="dat_gc.txt"){
  dat<-read.table(paste0(workingdir,"/",datfile),header = TRUE)
  dat_c<-read.table(paste0(workingdir,"/",datcfile),header = TRUE)
  dat$grandchildren<-0
  for (i in 1:nrow(dat)){
    dat$grandchildren[i]<-sum(dat_c$children[dat_c$parent1==dat$index[i]])+sum(dat_c$children[dat_c$parent2==dat$index[i]])
  }
  write.table(dat,paste0(workingdir,"/",datgcfile),quote=FALSE,row.names = FALSE)
}
## runs the function given a directory
args = as.list(commandArgs(TRUE))
names(args)<-c("workingdir","datfile","datcfile","datgcfile")[1:length(args)]
do.call(create_dat_gc_file,args)
#!/usr/bin/env Rscript

## function to create dat.txt file from an output.txt and parent.txt file in the given directory
create_dat_file<-function(workingdir,outfile="output.txt",parfile="parent.txt",datfile="dat.txt"){
  dat<-read.table(paste0(workingdir,"/",outfile),header = TRUE)
  par<-read.table(paste0(workingdir,"/",parfile),header = TRUE)
  temp<-as.data.frame(table(c(par$parent1,par$parent2)))
  colnames(temp)<-c("index","children")
  dat<-merge(dat,temp,all.x = TRUE)
  dat$children[is.na(dat$children)]<-0
  write.table(dat,paste0(workingdir,"/",datfile),quote=FALSE,row.names = FALSE)
}
## runs the function given a directory
args = as.list(commandArgs(TRUE))
names(args)<-c("workingdir","outfile","parfile","datfile")[1:length(args)]
do.call(create_dat_file,args)
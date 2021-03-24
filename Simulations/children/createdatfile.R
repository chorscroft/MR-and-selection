#!/usr/bin/env Rscript

## function to create dat.txt file from an output.txt and parent.txt file in the give directory
create_dat_file<-function(workingdir){
  dat<-read.table(paste0(workingdir,"/output.txt"),header = TRUE)
  par<-read.table(paste0(workingdir,"/parent.txt"),header = TRUE)
  temp<-as.data.frame(table(c(par$parent1,par$parent2)))
  colnames(temp)<-c("index","children")
  dat<-merge(dat,temp,all.x = TRUE)
  dat$children[is.na(dat$children)]<-0
  write.table(dat,paste0(workingdir,"/dat.txt"),quote=FALSE,row.names = FALSE)
}
## runs the function given a directory
create_dat_file(commandArgs(TRUE)[1])
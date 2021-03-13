#!/usr/bin/env Rscript

## function to create dat.txt file from an output.txt and parent.txt file in the give directory
create_dat_file<-function(workingdir){
  dat<-read.table(paste0(workingdir,"/output.txt"),header = TRUE)
  par<-read.table(paste0(workingdir,"/parent.txt"),header = TRUE)
  dat$children<-sapply(dat$index,function(x) sum(par$parent1==x)+sum(par$parent2==x))
  write.table(dat,paste0(workingdir,"/dat.txt"),quote=FALSE,row.names = FALSE)
}
## runs the function given a directory
create_dat_file(commandArgs(TRUE)[1])
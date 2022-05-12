library(simulateGP)

# Function to get child's genotype
get_child_g<-function(g1,g2){
  g3<-rep(0,length(g1))
  for (i in 1:length(g1)){
    g3[i]<-sample(c(0,1),1,T,prob=c(2-g1[i],g1[i]))+sample(c(0,1),1,T,prob=c(2-g2[i],g2[i]))
  }
  return(g3)
}

set.seed(1083219317)
nid<-1000
nsnps<-20
g<-list()
x<-list()
fitness<-list()
no_children<-list()
g[[1]]<-make_geno(nid, nsnps, 0.5)
eff<-choose_effects(nsnps,0.2)
x[[1]]<-make_phen(eff,g[[1]])
fitness[[1]]<-5+x[[1]]+rnorm(nid,0,0.1)


no_children<-list()
for(gen in 1:50){
  g[[gen+1]]<-matrix(0,nid,nsnps)
  no_children[[gen]]<-rep(0,nid)
  for (i in 1:nid){
    child_parents<-sample(1:nid,2,F,fitness[[gen]])
    g[[gen+1]][i,]<-get_child_g(g[[gen]][child_parents[1],],g[[gen]][child_parents[2],])
    no_children[[gen]][child_parents[1]]<-no_children[[gen]][child_parents[1]]+1
    no_children[[gen]][child_parents[2]]<-no_children[[gen]][child_parents[2]]+1
  }
  x[[gen+1]]<-make_phen(eff,g[[gen+1]])
  fitness[[gen+1]]<-5+x[[gen+1]]+rnorm(nid,0,0.1)
}

# change in x
delta_x<-sapply(x[1:50],mean)-sapply(x[2:51],mean)
cov_x_w<-sapply(1:50,function(gen)cov(x[[gen]],no_children[[gen]]))
plot(cov_x_w,delta_x)
       

# Parameters

beta1<-0.1
beta2<--0.05
year_of_change<-1980
decades<-5
distributions<-list(
  data.frame(age=c(10,20,30),prop=c(0.2,0.65,0.15)),
  data.frame(age=c(10,20,30),prop=c(0.15,0.7,0.15)),
  data.frame(age=c(10,20,30),prop=c(0.15,0.65,0.2)),
  data.frame(age=c(10,20,30),prop=c(0.15,0.6,0.25)),
  data.frame(age=c(10,20,30),prop=c(0.1,0.6,0.3))
)
birth_decades<-c(1940,1950,1960,1970,1980)

# If data told us birth date of children:
par(mar=c(5,4,0,2),mfrow=c(1,1),xpd=FALSE)
plot(c(seq(1950,2020,10)),c(rep(beta1,3),rep(beta2,5)),ylim=c(beta2-0.01,beta1+0.01),xlab="Decade of child's birth",ylab="Effect estimate")
abline(h=0,col="lightgrey")
abline(v=year_of_change,col="red",lty=2)

# Plot with a line through the year of change

par(oma=c(0,0,3,0),mar=c(2,4,0,2),mfrow=c(5,1),xpd=FALSE)
bp<-barplot(c(0,distributions[[1]]$prop,0,0,0,0,0),names.arg=seq(1940,2020,10),ylab="1940",xlab=NA)
abline(v=5,col="red",lty=2)
par(mar=c(2,4,0,2))
bp<-barplot(c(0,0,distributions[[2]]$prop,0,0,0,0),names.arg=seq(1940,2020,10),ylab="1950",xlab=NA)
abline(v=5,col="red",lty=2)
bp<-barplot(c(0,0,0,distributions[[3]]$prop,0,0,0),names.arg=seq(1940,2020,10),ylab="1960",xlab=NA)
abline(v=5,col="red",lty=2)
bp<-barplot(c(0,0,0,0,distributions[[4]]$prop,0,0),names.arg=seq(1940,2020,10),ylab="1970",xlab=NA)
abline(v=5,col="red",lty=2)
bp<-barplot(c(0,0,0,0,0,distributions[[5]]$prop,0),names.arg=seq(1940,2020,10),ylab="1980",xlab=NA)
abline(v=5,col="red",lty=2)
mtext("Proportion of children born in each decade",line=1, outer = TRUE, cex = 1)

# known decade effects
decade_effects<-data.frame(decade=birth_decades,effects=NA)
for(d in 1:length(birth_decades)){
  decade_effects$effects[d]<-ifelse(birth_decades[d]+10<year_of_change,beta1,beta2)*distributions[[d]]$prop[1] +
                              ifelse(birth_decades[d]+20<year_of_change,beta1,beta2)*distributions[[d]]$prop[2] +
                              ifelse(birth_decades[d]+30<year_of_change,beta1,beta2)*distributions[[d]]$prop[3]
}

decade_effects

min_beta_function<-function(betas,year,distributions,decade_effects,birth_decades){
  est_decade_effects<-rep(0,length(birth_decades))
  for(d in 1:length(birth_decades)){
    est_decade_effects[d]<-ifelse(birth_decades[d]+10<year,betas[1],betas[2])*distributions[[d]]$prop[1] +
      ifelse(birth_decades[d]+20<year,betas[1],betas[2])*distributions[[d]]$prop[2] +
      ifelse(birth_decades[d]+30<year,betas[1],betas[2])*distributions[[d]]$prop[3]
  }
  return(sum((est_decade_effects - decade_effects$effects)^2))
}


min_year_function<-function(year_range,distributions,decade_effects,birth_decades){
  year_optim<-vector("list",length(year_range))
  for (i in 1:length(year_optim)){
    year_optim[[i]]<-optim(c(0,0),min_beta_function,gr=NULL,year_range[i],distributions,decade_effects,birth_decades)[c('par','value')]
    year_optim[[i]]$year<-year_range[i]
  }
  values<-sapply(year_optim,function(x) x$value)
  return(year_optim[which(values==min(values))])
}
result<-min_year_function(seq(1950,2010,10),distributions,decade_effects,birth_decades)


# Simulate the number of children each person has in their lifetime
# Assume the 20,000 people born in each decade have 20,000 children

set.seed(2981847)
nid_per_decade<-20000
dataset<-data.frame(id=1:(nid_per_decade*decades),genotype=c(0,1),trait=c(0,1),birth_decade=rep(birth_decades,each=nid_per_decade),number_of_children=0)


for (d in 1:length(birth_decades)){
  parents<-dataset[dataset$birth_decade==birth_decades[d],]
  age_of_parent<-sample(c(10,20,30),nid_per_decade,T,distributions[[d]]$prop)
  p_10<-sample(parents$id,sum(age_of_parent==10),TRUE,exp(1+parents$trait*ifelse(birth_decades[d]+10<year_of_change,beta1,beta2)))
  p_20<-sample(parents$id,sum(age_of_parent==20),TRUE,exp(1+parents$trait*ifelse(birth_decades[d]+20<year_of_change,beta1,beta2)))
  p_30<-sample(parents$id,sum(age_of_parent==30),TRUE,exp(1+parents$trait*ifelse(birth_decades[d]+30<year_of_change,beta1,beta2)))
  tabParents <- table(c(p_10,p_20,p_30))
  dataset$number_of_children[as.numeric(unlist(dimnames(tabParents)))] <- tabParents
}

## Overall estimate
## genotype to trait
mod_gx<-summary(lm(dataset$trait~dataset$genotype))
b_gx<-mod_gx$coefficients[2,1]
se_gx<-mod_gx$coefficients[2,2]

## genotype to number of children (glm method)
mod_glm_gy<-summary(glm(dataset$number_of_children~dataset$genotype,family=poisson("log")))
b_glm_gy<-mod_glm_gy$coefficients[2,1]
se_glm_gy<-mod_glm_gy$coefficients[2,2]
## Wald Ratio
wr<-mr_wald_ratio(b_gx,b_glm_gy,se_gx,se_glm_gy)
wr

## Estimate by decade of birth
decade_effect_ests<-data.frame(decade=birth_decades,effects=NA,se=NA,pval=NA,nsnp=NA)
for (d in 1:length(birth_decades)){
  d_dataset<-dataset[dataset$birth_decade==birth_decades[d],]
  ## genotype to trait
  mod_gx<-summary(lm(d_dataset$trait~d_dataset$genotype))
  b_gx<-mod_gx$coefficients[2,1]
  se_gx<-mod_gx$coefficients[2,2]
  
  ## genotype to number of children (glm method)
  mod_glm_gy<-summary(glm(d_dataset$number_of_children~d_dataset$genotype,family=poisson("log")))
  b_glm_gy<-mod_glm_gy$coefficients[2,1]
  se_glm_gy<-mod_glm_gy$coefficients[2,2]
  ## Wald Ratio
  wr<-mr_wald_ratio(b_gx,b_glm_gy,se_gx,se_glm_gy)
  decade_effect_ests$effects[d]<-wr$b
  decade_effect_ests$se[d]<-wr$se
  decade_effect_ests$pval[d]<-wr$pval
  decade_effect_ests$nsnp[d]<-wr$nsnp
}

sim_result<-min_year_function(seq(1950,2010,10),distributions,decade_effect_ests,birth_decades)



##############
# Simulate with more noise

set.seed(2904723)
nid_per_decade<-20000
dataset<-data.frame(id=1:(nid_per_decade*decades),birth_decade=rep(birth_decades,each=nid_per_decade),number_of_children=0)
dataset$genotype<-rbinom(nrow(dataset),2,0.5)
dataset$trait<-5*dataset$genotype+rnorm(nrow(dataset))

for (d in 1:length(birth_decades)){
  parents<-dataset[dataset$birth_decade==birth_decades[d],]
  age_of_parent<-sample(c(10,20,30),nid_per_decade,T,distributions[[d]]$prop)
  p_10<-sample(parents$id,sum(age_of_parent==10),TRUE,exp(1+parents$trait*ifelse(birth_decades[d]+10<year_of_change,beta1,beta2)+rnorm(nrow(parents),0,0.1)))
  p_20<-sample(parents$id,sum(age_of_parent==20),TRUE,exp(1+parents$trait*ifelse(birth_decades[d]+20<year_of_change,beta1,beta2)+rnorm(nrow(parents),0,0.1)))
  p_30<-sample(parents$id,sum(age_of_parent==30),TRUE,exp(1+parents$trait*ifelse(birth_decades[d]+30<year_of_change,beta1,beta2)+rnorm(nrow(parents),0,0.1)))
  tabParents <- table(c(p_10,p_20,p_30))
  dataset$number_of_children[as.numeric(unlist(dimnames(tabParents)))] <- tabParents
}

## Overall estimate
## genotype to trait
mod_gx<-summary(lm(dataset$trait~dataset$genotype))
b_gx<-mod_gx$coefficients[2,1]
se_gx<-mod_gx$coefficients[2,2]

## genotype to number of children (glm method)
mod_glm_gy<-summary(glm(dataset$number_of_children~dataset$genotype,family=poisson("log")))
b_glm_gy<-mod_glm_gy$coefficients[2,1]
se_glm_gy<-mod_glm_gy$coefficients[2,2]
## Wald Ratio
wr<-mr_wald_ratio(b_gx,b_glm_gy,se_gx,se_glm_gy)
wr

## Estimate by decade of birth
decade_effect_ests<-data.frame(decade=birth_decades,effects=NA,se=NA,pval=NA,nsnp=NA)
for (d in 1:length(birth_decades)){
  d_dataset<-dataset[dataset$birth_decade==birth_decades[d],]
  ## genotype to trait
  mod_gx<-summary(lm(d_dataset$trait~d_dataset$genotype))
  b_gx<-mod_gx$coefficients[2,1]
  se_gx<-mod_gx$coefficients[2,2]
  
  ## genotype to number of children (glm method)
  mod_glm_gy<-summary(glm(d_dataset$number_of_children~d_dataset$genotype,family=poisson("log")))
  b_glm_gy<-mod_glm_gy$coefficients[2,1]
  se_glm_gy<-mod_glm_gy$coefficients[2,2]
  ## Wald Ratio
  wr<-mr_wald_ratio(b_gx,b_glm_gy,se_gx,se_glm_gy)
  decade_effect_ests$effects[d]<-wr$b
  decade_effect_ests$se[d]<-wr$se
  decade_effect_ests$pval[d]<-wr$pval
  decade_effect_ests$nsnp[d]<-wr$nsnp
}

sim_result<-min_year_function(seq(1950,2010,10),distributions,decade_effect_ests,birth_decades)




##############
# Simulate with a gradually decreasing beta
beta<-data.frame(years=c(1950,1960,1970,1980,1990,2000,2010),beta=c(0.1,0.075,0.05,0.025,0,-0.025,-0.05))

# If data told us birth date of children:
par(mar=c(5,4,0,2),mfrow=c(1,1),xpd=FALSE)
plot(beta$years,beta$beta,xlab="Decade of child's birth",ylab="Effect estimate")
abline(h=0,col="lightgrey")

set.seed(765775)
nid_per_decade<-20000
dataset<-data.frame(id=1:(nid_per_decade*decades),birth_decade=rep(birth_decades,each=nid_per_decade),number_of_children=0)
dataset$genotype<-rbinom(nrow(dataset),2,0.5)
dataset$trait<-5*dataset$genotype+rnorm(nrow(dataset))

for (d in 1:length(birth_decades)){
  parents<-dataset[dataset$birth_decade==birth_decades[d],]
  age_of_parent<-sample(c(10,20,30),nid_per_decade,T,distributions[[d]]$prop)
  p_10<-sample(parents$id,sum(age_of_parent==10),TRUE,exp(1+parents$trait*beta$beta[d]+rnorm(nrow(parents),0,0.1)))
  p_20<-sample(parents$id,sum(age_of_parent==20),TRUE,exp(1+parents$trait*beta$beta[d+1]+rnorm(nrow(parents),0,0.1)))
  p_30<-sample(parents$id,sum(age_of_parent==30),TRUE,exp(1+parents$trait*beta$beta[d+2]+rnorm(nrow(parents),0,0.1)))
  tabParents <- table(c(p_10,p_20,p_30))
  dataset$number_of_children[as.numeric(unlist(dimnames(tabParents)))] <- tabParents
}

## Overall estimate
## genotype to trait
mod_gx<-summary(lm(dataset$trait~dataset$genotype))
b_gx<-mod_gx$coefficients[2,1]
se_gx<-mod_gx$coefficients[2,2]

## genotype to number of children (glm method)
mod_glm_gy<-summary(glm(dataset$number_of_children~dataset$genotype,family=poisson("log")))
b_glm_gy<-mod_glm_gy$coefficients[2,1]
se_glm_gy<-mod_glm_gy$coefficients[2,2]
## Wald Ratio
wr<-mr_wald_ratio(b_gx,b_glm_gy,se_gx,se_glm_gy)
wr

## Estimate by decade of birth
decade_effect_ests<-data.frame(decade=birth_decades,effects=NA,se=NA,pval=NA,nsnp=NA)
for (d in 1:length(birth_decades)){
  d_dataset<-dataset[dataset$birth_decade==birth_decades[d],]
  ## genotype to trait
  mod_gx<-summary(lm(d_dataset$trait~d_dataset$genotype))
  b_gx<-mod_gx$coefficients[2,1]
  se_gx<-mod_gx$coefficients[2,2]
  
  ## genotype to number of children (glm method)
  mod_glm_gy<-summary(glm(d_dataset$number_of_children~d_dataset$genotype,family=poisson("log")))
  b_glm_gy<-mod_glm_gy$coefficients[2,1]
  se_glm_gy<-mod_glm_gy$coefficients[2,2]
  ## Wald Ratio
  wr<-mr_wald_ratio(b_gx,b_glm_gy,se_gx,se_glm_gy)
  decade_effect_ests$effects[d]<-wr$b
  decade_effect_ests$se[d]<-wr$se
  decade_effect_ests$pval[d]<-wr$pval
  decade_effect_ests$nsnp[d]<-wr$nsnp
}

sim_result<-min_year_function(seq(1950,2010,10),distributions,decade_effect_ests,birth_decades)

# beta1 is beta in first decade, beta2 is the difference between decades
min_linear_beta_function<-function(betas,distributions,decade_effects,birth_decades){
  est_decade_effects<-rep(0,length(birth_decades))
  for(d in 1:length(birth_decades)){
    est_decade_effects[d]<-(betas[1]+betas[2]*(d-1))*distributions[[d]]$prop[1] +
      (betas[1]+betas[2]*d)*distributions[[d]]$prop[2] +
      (betas[1]+betas[2]*(d+1))*distributions[[d]]$prop[3]
  }
  return(sum((est_decade_effects - decade_effects$effects)^2))
}


sim_linear_result<-optim(c(0,0),min_linear_beta_function,gr=NULL,distributions,decade_effect_ests,birth_decades)[c('par','value')]

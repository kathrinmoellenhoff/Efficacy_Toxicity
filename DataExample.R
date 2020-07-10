######################################################################
#########Code for evaluating the data example from the paper##########
######"Equivalence tests for binary efficacy-toxicity responses"######
#############method implemented can also be used for #################
##########reproducing the simulation results##########################
####################Author: Kathrin Möllenhoff########################

library(MultiOrd) #package for generation of binary data
library(alabama) #package for constrained optimization

#(bivariate) gumbel cdf
gumbel <- function(nu){
  cdf <- function(u,v){1/(1+exp(-u))*1/(1+exp(-v))*(1+(nu*exp(-u-v))/((1+exp(-u))*(1+exp(-v))))
  }
  return(cdf)
}

#standardized doses, see definition in the paper 
dose <- function(beta,gamma){
  dose <- function(d){beta+gamma*d}
  return(dose)
}

#marginal prob. curve
marg_prob <- function(beta,gamma){return(function(d){1/(1+exp(-dose(beta,gamma)(d)))})} 

#cell prob. curves
p11 <- function(d){function(beta1,gamma1,beta2,gamma2,nu){gumbel(nu)(dose(beta1,gamma1)(d),dose(beta2,gamma2)(d))}}
p10 <- function(d){function(beta1,gamma1,beta2,gamma2,nu){marg_prob(beta1,gamma1)(d)-gumbel(nu)(dose(beta1,gamma1)(d),dose(beta2,gamma2)(d))}}
p01 <- function(d){function(beta1,gamma1,beta2,gamma2,nu){marg_prob(beta2,gamma2)(d)-gumbel(nu)(dose(beta1,gamma1)(d),dose(beta2,gamma2)(d))}}
p00 <- function(d){function(beta1,gamma1,beta2,gamma2,nu){1-marg_prob(beta1,gamma1)(d)-marg_prob(beta2,gamma2)(d)+gumbel(nu)(dose(beta1,gamma1)(d),dose(beta2,gamma2)(d))}}

#data-generating function
n <- 30 #patients per dose level
gen_data <- function(doses,beta1,gamma1,beta2,gamma2,nu){
  for(i in 1:length(doses)){
    
    do <- doses[i]  
    d1 <- dose(beta1,gamma1)(do)
    d2 <- dose(beta2,gamma2)(do)
    
    e1 <- 1/(1+exp(-d1))
    e2 <- 1/(1+exp(-d2))
    a <- (nu*exp(-d1-d2))/((1+exp(-d1))^2*(1+exp(-d2))^2)
    
    #four cell probabilities in dependence of dose 
    p11s <- e1*e2+a
    p10s <- e1-e1*e2-a
    p01s <- e2-e1*e2-a
    p00s <- 1-e1-e2+e1*e2+a
    
    #saving cell probs
    cellp <- matrix(c(p00s,p01s,p10s,p11s),nrow=2,ncol=2,byrow=TRUE) 
    cellp_list1[[i]] <- cellp
    
    margp <- c(e1,e2) #marginal probabilities
    oddsr <- (p00s*p01s)/(p10s*p11s) #odds ratio
    
    corr <- nu/((exp(d1/2)+exp(-d1/2))*(exp(d2/2)+exp(-d2/2))) #correlation (in dependence of dose)
    cr <- matrix(c(1,corr,corr,1),nrow=2,ncol=2) #correlation matrix
    
    sim_data <- generate.binary(nObs=n,prop.vec.bin=margp,corr.mat=cr)
    
    data1[[i]] <- sim_data
    cellprob1[i,] <- c(p00s,p01s,p10s,p11s)
  }
  return(data1)
}

#dose levels under consideration
doses <- c(0,0.05,0.2,0.5,1)
doses2 <- c(0,0.1,0.3,0.6,1)

#create a grid on (0,1) to search for the maximum 
y <- seq(min(doses),max(doses),0.01)

#for storing results
cellprob1 <- matrix(nrow=length(doses),ncol=4)
cellp_list1 <- list()
cellprob2 <- matrix(nrow=length(doses),ncol=4)
cellp_list2 <- list()

#number of bootstrap repetitions and vectors for storing results
B <- 1000
boot <- vector()
boot2 <- vector()
critval <- vector()
pval <- vector()

#functions yielding inequality constraints for the estimation of the Gumbel model in order to guarantee that a joint distribution can exist
hin1 <- function(x) {
  p1 <- 1/(1+exp(-dose(x[1],x[2])(doses)))
  p2 <- 1/(1+exp(-dose(x[3],x[4])(doses)))
  corr_c <- x[5]/(exp(dose(x[1],x[2])(doses)/2)+exp(-dose(x[1],x[2])(doses)/2)*(exp(dose(x[3],x[4])(doses)/2)+exp(-dose(x[3],x[4])(doses)/2)))
  h <- rep(NA, 1)
  h[1] <- x[5]+((1+exp(-dose(x[1],x[2])(min(doses))))*(1+exp(-dose(x[3],x[4])(min(doses)))))/(exp(-dose(x[1],x[2])(min(doses))-dose(x[3],x[4])(min(doses))))
  h[2] <- min(corr_c-max(-sqrt(p1*p2/((1-p1)*(1-p2))),-sqrt((1-p1)*(1-p2)/(p1*p2))))
  h[3] <- min(min(sqrt((p1*(1-p2))/((1-p1)*p2)),sqrt((1-p1)*p2/((p1*(1-p2)))))-corr_c)
  h
} 

#joint version of hin1 for a simultaneous estimation of both models (for the constrained optimization)
hin <- function(x) {
  p1 <- 1/(1+exp(-dose(x[1],x[2])(doses)))
  p2 <- 1/(1+exp(-dose(x[3],x[4])(doses)))
  corr_c <- x[5]/(exp(dose(x[1],x[2])(doses)/2)+exp(-dose(x[1],x[2])(doses)/2)*(exp(dose(x[3],x[4])(doses)/2)+exp(-dose(x[3],x[4])(doses)/2)))
  p1b <- 1/(1+exp(-dose(x[6],x[7])(doses)))
  p2b <- 1/(1+exp(-dose(x[8],x[9])(doses)))
  corr_cb <- x[10]/(exp(dose(x[6],x[7])(doses)/2)+exp(-dose(x[6],x[7])(doses)/2)*(exp(dose(x[8],x[9])(doses)/2)+exp(-dose(x[8],x[9])(doses)/2)))
  h <- rep(NA, 1)
  h[1] <- x[5]+((1+exp(-dose(x[1],x[2])(min(doses))))*(1+exp(-dose(x[3],x[4])(min(doses)))))/(exp(-dose(x[1],x[2])(min(doses))-dose(x[3],x[4])(min(doses))))
  h[2] <- x[10]+((1+exp(-dose(x[6],x[7])(min(doses))))*(1+exp(-dose(x[8],x[9])(min(doses)))))/(exp(-dose(x[6],x[7])(min(doses))-dose(x[8],x[9])(min(doses))))
  h[3] <- min(corr_c-max(-sqrt(p1*p2/((1-p1)*(1-p2))),-sqrt((1-p1)*(1-p2)/(p1*p2))))
  h[4] <- min(corr_cb-max(-sqrt(p1b*p2b/((1-p1b)*(1-p2b))),-sqrt((1-p1b)*(1-p2b)/(p1b*p2b))))
  h[5] <- min(min(sqrt((p1b*(1-p2b))/((1-p1b)*p2b)),sqrt((1-p1b)*p2b/((p1b*(1-p2b)))))-corr_cb)
  h[6] <- min(min(sqrt((p1*(1-p2))/((1-p1)*p2)),sqrt((1-p1)*p2/((p1*(1-p2)))))-corr_c)
  h
} 

#function for constrained optimization producing models with the same beta parameters (placebos)
shared_placebos <- function(v){
  beta1=v[1]
  gamma1=v[2]
  beta2=v[3]
  gamma2=v[4]
  nu=v[5]
  beta1b=v[6]
  gamma1b=v[7]
  beta2b=v[8]
  gamma2b=v[9]
  nub=v[10]
  h <- rep(NA, 1)
  h[1]=beta1-beta1b
  h[2]=beta2-beta2b
  h
}

#load datasets from the Data Example
dataset <- readRDS("CaseStudyData.rds")
data1 <- dataset$d1
data2 <- dataset$d2

######################################################
######################################################
#Option 1: two seperate models
#fit a Gumbel model to each treatment group

#group 1
#MLE
sum <- 0
likelihood1 <- function(v){
  beta1=v[1]
  gamma1=v[2]
  beta2=v[3]
  gamma2=v[4]
  nu=v[5]
  for(i in 1:length(doses)){
    d <- doses[i]  
    sum <- sum + log(p11(d)(beta1,gamma1,beta2,gamma2,nu)^sum(data1[[i]][,1]==1 & data1[[i]][,2]==1)*p01(d)(beta1,gamma1,beta2,gamma2,nu)^sum(data1[[i]][,1]==0 & data1[[i]][,2]==1)*p10(d)(beta1,gamma1,beta2,gamma2,nu)^sum(data1[[i]][,1]==1 & data1[[i]][,2]==0)*p00(d)(beta1,gamma1,beta2,gamma2,nu)^sum(data1[[i]][,1]==0 & data1[[i]][,2]==0))
  }
  return(-sum) #for minimizing
}
par_gumbel1 <- auglag(c(0,1,0,1,0.5),likelihood1,hin=hin1,control.outer=list(method="nlminb",trace=FALSE))$par  

#group2
#MLE
sum <- 0
likelihood2 <- function(v){
  beta1=v[1]
  gamma1=v[2]
  beta2=v[3]
  gamma2=v[4]
  nu=v[5]
  for(i in 1:length(doses2)){
    d <- doses2[i]  
    sum <- sum + log(p11(d)(beta1,gamma1,beta2,gamma2,nu)^sum(data2[[i]][,1]==1 & data2[[i]][,2]==1)*p01(d)(beta1,gamma1,beta2,gamma2,nu)^sum(data2[[i]][,1]==0 & data2[[i]][,2]==1)*p10(d)(beta1,gamma1,beta2,gamma2,nu)^sum(data2[[i]][,1]==1 & data2[[i]][,2]==0)*p00(d)(beta1,gamma1,beta2,gamma2,nu)^sum(data2[[i]][,1]==0 & data2[[i]][,2]==0))
  }
  return(-sum) #for minimizing
}
par_gumbel2 <- auglag(c(0,1,0,1,0.5),likelihood2,hin=hin1,control.outer=list(method="nlminb",trace=FALSE))$par  

#calculate test statistics 
d_e_hat <- max(abs(marg_prob(par_gumbel1[1],par_gumbel1[2])(y)-marg_prob(par_gumbel2[1],par_gumbel2[2])(y))) #max.abs.distance efficacy
d_t_hat <- max(abs(marg_prob(par_gumbel1[3],par_gumbel1[4])(y)-marg_prob(par_gumbel2[3],par_gumbel2[4])(y))) #max. abs distance toxicity
t.stat <- d_e_hat
t.stat2 <- d_t_hat
######################################################
#bootstrap test for three different equivalence margins
######################################################

for(epsilon in c(0.1,0.15,0.2)){
#constrained optimization for the efficacy curves
if (t.stat>=epsilon){minimum <- c(par_gumbel1,par_gumbel2)} else {
  sum <- 0
  joint_likelihood=function(v){
    beta1=v[1]
    gamma1=v[2]
    beta2=v[3]
    gamma2=v[4]
    nu=v[5]
    beta1b=v[6]
    gamma1b=v[7]
    beta2b=v[8]
    gamma2b=v[9]
    nub=v[10]
    for(i in 1:5){
      d <- c(doses,doses2)[i]
      d2 <- c(doses,doses2)[5+i]
      sum <- sum +
        log(p11(d)(beta1,gamma1,beta2,gamma2,nu)^sum(data1[[i]][,1]==1 & data1[[i]][,2]==1)*p01(d)(beta1,gamma1,beta2,gamma2,nu)^sum(data1[[i]][,1]==0 & data1[[i]][,2]==1)*p10(d)(beta1,gamma1,beta2,gamma2,nu)^sum(data1[[i]][,1]==1 & data1[[i]][,2]==0)*p00(d)(beta1,gamma1,beta2,gamma2,nu)^sum(data1[[i]][,1]==0 & data1[[i]][,2]==0))+
        log(p11(d2)(beta1b,gamma1b,beta2b,gamma2b,nub)^sum(data2[[i]][,1]==1 & data2[[i]][,2]==1)*p01(d2)(beta1b,gamma1b,beta2b,gamma2b,nub)^sum(data2[[i]][,1]==0 & data2[[i]][,2]==1)*p10(d2)(beta1b,gamma1b,beta2b,gamma2b,nub)^sum(data2[[i]][,1]==1 & data2[[i]][,2]==0)*p00(d2)(beta1b,gamma1b,beta2b,gamma2b,nub)^sum(data2[[i]][,1]==0 & data2[[i]][,2]==0))
    }
    return(-sum)
  }
  softmax <- function(v){
    beta1=v[1]
    gamma1=v[2]
    beta2=v[3]
    gamma2=v[4]
    nu=v[5]
    beta1b=v[6]
    gamma1b=v[7]
    beta2b=v[8]
    gamma2b=v[9]
    nub=v[10]
    diff_e=function(x){abs(marg_prob(beta1,gamma1)(x)-marg_prob(beta1b,gamma1b)(x))}
    (sum(diff_e(y)*exp(100*diff_e(y))))/(sum(exp(100*diff_e(y))))-epsilon}
  minimum <- auglag(par=c(par_gumbel1,par_gumbel2),joint_likelihood,heq=softmax,hin=hin,control.outer=list(method="nlminb",trace=FALSE))$par
}

#bootstrap
#generate bootstrap data using the parameters obtained by the constrained optimization
data1b <- list()
data2b <- list()

for (k in 1:B)
{
  data1b <- gen_data(doses,minimum[1],minimum[2],minimum[3],minimum[4],minimum[5])
  
  #MLE
  sum <- 0
  likelihood1b <- function(v){
    beta1=v[1]
    gamma1=v[2]
    beta2=v[3]
    gamma2=v[4]
    nu=v[5]
    for(i in 1:length(doses)){
      d <- doses[i]  
      sum <- sum + log(p11(d)(beta1,gamma1,beta2,gamma2,nu)^sum(data1b[[i]][,1]==1 & data1b[[i]][,2]==1)*p01(d)(beta1,gamma1,beta2,gamma2,nu)^sum(data1b[[i]][,1]==0 & data1b[[i]][,2]==1)*p10(d)(beta1,gamma1,beta2,gamma2,nu)^sum(data1b[[i]][,1]==1 & data1b[[i]][,2]==0)*p00(d)(beta1,gamma1,beta2,gamma2,nu)^sum(data1b[[i]][,1]==0 & data1b[[i]][,2]==0))
    }
    return(-sum) #for minimizing
  }
  par_gumbel1b <- optim(minimum[1:5],likelihood1b)$par  

  #group2
  data2b <- gen_data(doses2,minimum[6],minimum[7],minimum[8],minimum[9],minimum[10])
  
  #MLE
  sum <- 0
  likelihood2b <- function(v){
    beta1=v[1]
    gamma1=v[2]
    beta2=v[3]
    gamma2=v[4]
    nu=v[5]
    for(i in 1:length(doses2)){
      d <- doses2[i]  
      sum <- sum + log(p11(d)(beta1,gamma1,beta2,gamma2,nu)^sum(data2b[[i]][,1]==1 & data2b[[i]][,2]==1)*p01(d)(beta1,gamma1,beta2,gamma2,nu)^sum(data2b[[i]][,1]==0 & data2b[[i]][,2]==1)*p10(d)(beta1,gamma1,beta2,gamma2,nu)^sum(data2b[[i]][,1]==1 & data2b[[i]][,2]==0)*p00(d)(beta1,gamma1,beta2,gamma2,nu)^sum(data2b[[i]][,1]==0 & data2b[[i]][,2]==0))
    }
    return(-sum) #for minimizing
  }
  par_gumbel2b <- optim(minimum[6:10],likelihood2b)$par

  d_e_star <- max(abs(marg_prob(par_gumbel1b[1],par_gumbel1b[2])(y)-marg_prob(par_gumbel2b[1],par_gumbel2b[2])(y))) #max.abs.distance efficacy
  d_t_star <- max(abs(marg_prob(par_gumbel1b[3],par_gumbel1b[4])(y)-marg_prob(par_gumbel2b[3],par_gumbel2b[4])(y))) #max. abs distance toxicity
  t.star <- d_e_star
  
  boot[k] <- t.star
  #save time if a seperate bootstrap for toxicity is not needed
  #conditions imply that fot both tests the constrained estimates are the same as the unconstrained
  if (t.stat2>=epsilon & t.stat>=epsilon){boot2[k] <- d_t_star}
}

crit_val<-quantile(boot,0.05)
pval1 <- ecdf(boot)(t.stat)

#second bootstrap (toxicity curves)
if (t.stat2<epsilon | t.stat<epsilon){
  if (t.stat2>=epsilon){minimum <- c(par_gumbel1,par_gumbel2)} else {
    sum <- 0
    joint_likelihood=function(v){
      beta1=v[1]
      gamma1=v[2]
      beta2=v[3]
      gamma2=v[4]
      nu=v[5]
      beta1b=v[6]
      gamma1b=v[7]
      beta2b=v[8]
      gamma2b=v[9]
      nub=v[10]
      for(i in 1:5){
        d <- c(doses,doses2)[i]
        d2 <- c(doses,doses2)[5+i]
        sum <- sum +
          log(p11(d)(beta1,gamma1,beta2,gamma2,nu)^sum(data1[[i]][,1]==1 & data1[[i]][,2]==1)*p01(d)(beta1,gamma1,beta2,gamma2,nu)^sum(data1[[i]][,1]==0 & data1[[i]][,2]==1)*p10(d)(beta1,gamma1,beta2,gamma2,nu)^sum(data1[[i]][,1]==1 & data1[[i]][,2]==0)*p00(d)(beta1,gamma1,beta2,gamma2,nu)^sum(data1[[i]][,1]==0 & data1[[i]][,2]==0))+
          log(p11(d2)(beta1b,gamma1b,beta2b,gamma2b,nub)^sum(data2[[i]][,1]==1 & data2[[i]][,2]==1)*p01(d2)(beta1b,gamma1b,beta2b,gamma2b,nub)^sum(data2[[i]][,1]==0 & data2[[i]][,2]==1)*p10(d2)(beta1b,gamma1b,beta2b,gamma2b,nub)^sum(data2[[i]][,1]==1 & data2[[i]][,2]==0)*p00(d2)(beta1b,gamma1b,beta2b,gamma2b,nub)^sum(data2[[i]][,1]==0 & data2[[i]][,2]==0))
      }
      return(-sum)
    }
    softmax <- function(v){
      beta1=v[1]
      gamma1=v[2]
      beta2=v[3]
      gamma2=v[4]
      nu=v[5]
      beta1b=v[6]
      gamma1b=v[7]
      beta2b=v[8]
      gamma2b=v[9]
      nub=v[10]
      diff_t=function(x){abs(marg_prob(beta2,gamma2)(x)-marg_prob(beta2b,gamma2b)(x))}
      (sum(diff_t(y)*exp(100*diff_t(y))))/(sum(exp(100*diff_t(y))))-epsilon}
    minimum <- auglag(par=c(par_gumbel1,par_gumbel2),joint_likelihood,heq=softmax,hin=hin,control.outer=list(method="nlminb",trace=FALSE))$par
  }
  
  data1b <- list()
  data2b <- list()
  
  for (k in 1:B)
  {
    data1b <- gen_data(doses,minimum[1],minimum[2],minimum[3],minimum[4],minimum[5])
    
    #MLE
    sum <- 0
    likelihood1b <- function(v){
      beta1=v[1]
      gamma1=v[2]
      beta2=v[3]
      gamma2=v[4]
      nu=v[5]
      for(i in 1:length(doses)){
        d <- doses[i]  
        sum <- sum + log(p11(d)(beta1,gamma1,beta2,gamma2,nu)^sum(data1b[[i]][,1]==1 & data1b[[i]][,2]==1)*p01(d)(beta1,gamma1,beta2,gamma2,nu)^sum(data1b[[i]][,1]==0 & data1b[[i]][,2]==1)*p10(d)(beta1,gamma1,beta2,gamma2,nu)^sum(data1b[[i]][,1]==1 & data1b[[i]][,2]==0)*p00(d)(beta1,gamma1,beta2,gamma2,nu)^sum(data1b[[i]][,1]==0 & data1b[[i]][,2]==0))
      }
      return(-sum) #for minimizing
    }
    par_gumbel1b <- optim(minimum[1:5],likelihood1b)$par  

    #group2
    data2b <- gen_data(doses2,minimum[6],minimum[7],minimum[8],minimum[9],minimum[10])
    
    #MLE
    sum <- 0
    likelihood2b <- function(v){
      beta1=v[1]
      gamma1=v[2]
      beta2=v[3]
      gamma2=v[4]
      nu=v[5]
      for(i in 1:length(doses2)){
        d <- doses2[i]  
        sum <- sum + log(p11(d)(beta1,gamma1,beta2,gamma2,nu)^sum(data2b[[i]][,1]==1 & data2b[[i]][,2]==1)*p01(d)(beta1,gamma1,beta2,gamma2,nu)^sum(data2b[[i]][,1]==0 & data2b[[i]][,2]==1)*p10(d)(beta1,gamma1,beta2,gamma2,nu)^sum(data2b[[i]][,1]==1 & data2b[[i]][,2]==0)*p00(d)(beta1,gamma1,beta2,gamma2,nu)^sum(data2b[[i]][,1]==0 & data2b[[i]][,2]==0))
      }
      return(-sum) #for minimizing
    }
    par_gumbel2b <- optim(minimum[6:10],likelihood2b)$par   

    d_e_star <- max(abs(marg_prob(par_gumbel1b[1],par_gumbel1b[2])(y)-marg_prob(par_gumbel2b[1],par_gumbel2b[2])(y))) #max.abs.distance efficacy
    d_t_star <- max(abs(marg_prob(par_gumbel1b[3],par_gumbel1b[4])(y)-marg_prob(par_gumbel2b[3],par_gumbel2b[4])(y))) #max. abs distance toxicity
    t.star <- d_t_star
    
    boot2[k] <- t.star
  }
}
crit_val2 <- quantile(boot2,0.05)
pval2 <- ecdf(boot2)(t.stat2)

critval <- c(critval,crit_val,crit_val2)
pval <- c(pval,pval1,pval2)
}

######################################################
######################################################
########### Option 2 (shared placebos) ###############
######################################################
#bootstrap test for three different equivalence margins
######################################################
######################################################

sum <- 0
joint_likelihood=function(v){
  beta1=v[1]
  gamma1=v[2]
  beta2=v[3]
  gamma2=v[4]
  nu=v[5]
  beta1b=v[6]
  gamma1b=v[7]
  beta2b=v[8]
  gamma2b=v[9]
  nub=v[10]
  for(i in 1:5){
    d <- c(doses,doses2)[i]
    d2 <- c(doses,doses2)[5+i]
    sum <- sum +
      log(p11(d)(beta1,gamma1,beta2,gamma2,nu)^sum(data1[[i]][,1]==1 & data1[[i]][,2]==1)*p01(d)(beta1,gamma1,beta2,gamma2,nu)^sum(data1[[i]][,1]==0 & data1[[i]][,2]==1)*p10(d)(beta1,gamma1,beta2,gamma2,nu)^sum(data1[[i]][,1]==1 & data1[[i]][,2]==0)*p00(d)(beta1,gamma1,beta2,gamma2,nu)^sum(data1[[i]][,1]==0 & data1[[i]][,2]==0))+
      log(p11(d2)(beta1b,gamma1b,beta2b,gamma2b,nub)^sum(data2[[i]][,1]==1 & data2[[i]][,2]==1)*p01(d2)(beta1b,gamma1b,beta2b,gamma2b,nub)^sum(data2[[i]][,1]==0 & data2[[i]][,2]==1)*p10(d2)(beta1b,gamma1b,beta2b,gamma2b,nub)^sum(data2[[i]][,1]==1 & data2[[i]][,2]==0)*p00(d2)(beta1b,gamma1b,beta2b,gamma2b,nub)^sum(data2[[i]][,1]==0 & data2[[i]][,2]==0))
  }
  return(-sum)
}

par_gumbel <- auglag(par=c(0,1,0,1,0.5,0,1,0,1,0.5),joint_likelihood,heq=shared_placebos,hin=hin,control.outer=list(method="nlminb",trace=FALSE))$par
par_gumbel1 <- par_gumbel[1:5]
par_gumbel2 <- par_gumbel[6:10]

######################################################
######################################################

#calculate test statistics 
d_e_hat <- max(abs(marg_prob(par_gumbel1[1],par_gumbel1[2])(y)-marg_prob(par_gumbel2[1],par_gumbel2[2])(y))) #max.abs.distance efficacy
d_t_hat <- max(abs(marg_prob(par_gumbel1[3],par_gumbel1[4])(y)-marg_prob(par_gumbel2[3],par_gumbel2[4])(y))) #max. abs distance toxicity
t.stat <- d_e_hat
t.stat2 <- d_t_hat

for(epsilon in c(0.1,0.15,0.2)){
  #constrained optimization for the efficacy curves
  if (t.stat>=epsilon){minimum <- c(par_gumbel1,par_gumbel2)} else {
    sum <- 0
    joint_likelihood=function(v){
      beta1=v[1]
      gamma1=v[2]
      beta2=v[3]
      gamma2=v[4]
      nu=v[5]
      beta1b=v[6]
      gamma1b=v[7]
      beta2b=v[8]
      gamma2b=v[9]
      nub=v[10]
      for(i in 1:5){
        d <- c(doses,doses2)[i]
        d2 <- c(doses,doses2)[5+i]
        sum <- sum +
          log(p11(d)(beta1,gamma1,beta2,gamma2,nu)^sum(data1[[i]][,1]==1 & data1[[i]][,2]==1)*p01(d)(beta1,gamma1,beta2,gamma2,nu)^sum(data1[[i]][,1]==0 & data1[[i]][,2]==1)*p10(d)(beta1,gamma1,beta2,gamma2,nu)^sum(data1[[i]][,1]==1 & data1[[i]][,2]==0)*p00(d)(beta1,gamma1,beta2,gamma2,nu)^sum(data1[[i]][,1]==0 & data1[[i]][,2]==0))+
          log(p11(d2)(beta1b,gamma1b,beta2b,gamma2b,nub)^sum(data2[[i]][,1]==1 & data2[[i]][,2]==1)*p01(d2)(beta1b,gamma1b,beta2b,gamma2b,nub)^sum(data2[[i]][,1]==0 & data2[[i]][,2]==1)*p10(d2)(beta1b,gamma1b,beta2b,gamma2b,nub)^sum(data2[[i]][,1]==1 & data2[[i]][,2]==0)*p00(d2)(beta1b,gamma1b,beta2b,gamma2b,nub)^sum(data2[[i]][,1]==0 & data2[[i]][,2]==0))
      }
      return(-sum)
    }
    softmax <- function(v){
      beta1=v[1]
      gamma1=v[2]
      beta2=v[3]
      gamma2=v[4]
      nu=v[5]
      beta1b=v[6]
      gamma1b=v[7]
      beta2b=v[8]
      gamma2b=v[9]
      nub=v[10]
      diff_e=function(x){abs(marg_prob(beta1,gamma1)(x)-marg_prob(beta1b,gamma1b)(x))}
      h <- rep(NA, 1)
      h[1]=beta1-beta1b
      h[2]=beta2-beta2b
      h[3]=(sum(diff_e(y)*exp(100*diff_e(y))))/(sum(exp(100*diff_e(y))))-epsilon
      h
      }
    minimum <- auglag(par=c(par_gumbel1,par_gumbel2),joint_likelihood,heq=softmax,hin=hin,control.outer=list(method="nlminb",trace=FALSE))$par
  }
  
  #bootstrap
  #generate bootstrap data using the parameters obtained by the constrained optimization
  data1b <- list()
  data2b <- list()
  
  for (k in 1:B)
  {
    data1b <- gen_data(doses,minimum[1],minimum[2],minimum[3],minimum[4],minimum[5])
    data2b <- gen_data(doses2,minimum[6],minimum[7],minimum[8],minimum[9],minimum[10])
    
    #MLE
    #only 8 parameters due to the shared placebos
    sum <- 0
    likelihoodb <- function(v){
      beta1=v[1]
      gamma1=v[2]
      beta2=v[3]
      gamma2=v[4]
      nu=v[5]
      beta1b=v[6]
      gamma1b=v[7]
      beta2b=v[8]
      gamma2b=v[9]
      nub=v[10]
      for(i in 1:length(doses)){
        d <- doses[i]  
        sum <- sum + log(p11(d)(beta1,gamma1,beta2,gamma2,nu)^sum(data1b[[i]][,1]==1 & data1b[[i]][,2]==1)*p01(d)(beta1,gamma1,beta2,gamma2,nu)^sum(data1b[[i]][,1]==0 & data1b[[i]][,2]==1)*p10(d)(beta1,gamma1,beta2,gamma2,nu)^sum(data1b[[i]][,1]==1 & data1b[[i]][,2]==0)*p00(d)(beta1,gamma1,beta2,gamma2,nu)^sum(data1b[[i]][,1]==0 & data1b[[i]][,2]==0))
      + log(p11(d)(beta1,gamma1b,beta2,gamma2b,nub)^sum(data2b[[i]][,1]==1 & data2b[[i]][,2]==1)*p01(d)(beta1,gamma1b,beta2,gamma2b,nub)^sum(data2b[[i]][,1]==0 & data2b[[i]][,2]==1)*p10(d)(beta1,gamma1b,beta2,gamma2b,nub)^sum(data2b[[i]][,1]==1 & data2b[[i]][,2]==0)*p00(d)(beta1,gamma1b,beta2,gamma2b,nub)^sum(data2b[[i]][,1]==0 & data2b[[i]][,2]==0))
      }
      return(-sum) #for minimizing
    }
    par_gumbelb <- auglag(par=minimum,likelihoodb,heq=shared_placebos,hin=hin,control.outer=list(method="nlminb",trace=FALSE))$par
    par_gumbel1b <- par_gumbelb[1:5]
    par_gumbel2b <- par_gumbelb[6:10]
    
    
    d_e_star <- max(abs(marg_prob(par_gumbel1b[1],par_gumbel1b[2])(y)-marg_prob(par_gumbel2b[1],par_gumbel2b[2])(y))) #max.abs.distance efficacy
    d_t_star <- max(abs(marg_prob(par_gumbel1b[3],par_gumbel1b[4])(y)-marg_prob(par_gumbel2b[3],par_gumbel2b[4])(y))) #max. abs distance toxicity
    t.star <- d_e_star
    
    boot[k] <- t.star

    #save time if a seperate bootstrap for toxicity is not needed
    #conditions imply that fot both tests the constrained estimates are the same as the unconstrained
    if (t.stat2>=epsilon & t.stat>=epsilon){boot2[k] <- d_t_star}
  }
  
  crit_val<-quantile(boot,0.05)
  pval1 <- ecdf(boot)(t.stat)
  
  #second bootstrap (toxicity curves)
  if (t.stat2<epsilon | t.stat<epsilon){
    if (t.stat2>=epsilon){minimum <- c(par_gumbel1,par_gumbel2)} else {
      sum <- 0
      joint_likelihood=function(v){
        beta1=v[1]
        gamma1=v[2]
        beta2=v[3]
        gamma2=v[4]
        nu=v[5]
        beta1b=v[6]
        gamma1b=v[7]
        beta2b=v[8]
        gamma2b=v[9]
        nub=v[10]
        for(i in 1:5){
          d <- c(doses,doses2)[i]
          d2 <- c(doses,doses2)[5+i]
          sum <- sum +
            log(p11(d)(beta1,gamma1,beta2,gamma2,nu)^sum(data1[[i]][,1]==1 & data1[[i]][,2]==1)*p01(d)(beta1,gamma1,beta2,gamma2,nu)^sum(data1[[i]][,1]==0 & data1[[i]][,2]==1)*p10(d)(beta1,gamma1,beta2,gamma2,nu)^sum(data1[[i]][,1]==1 & data1[[i]][,2]==0)*p00(d)(beta1,gamma1,beta2,gamma2,nu)^sum(data1[[i]][,1]==0 & data1[[i]][,2]==0))+
            log(p11(d2)(beta1b,gamma1b,beta2b,gamma2b,nub)^sum(data2[[i]][,1]==1 & data2[[i]][,2]==1)*p01(d2)(beta1b,gamma1b,beta2b,gamma2b,nub)^sum(data2[[i]][,1]==0 & data2[[i]][,2]==1)*p10(d2)(beta1b,gamma1b,beta2b,gamma2b,nub)^sum(data2[[i]][,1]==1 & data2[[i]][,2]==0)*p00(d2)(beta1b,gamma1b,beta2b,gamma2b,nub)^sum(data2[[i]][,1]==0 & data2[[i]][,2]==0))
        }
        return(-sum)
      }
      softmax <- function(v){
        beta1=v[1]
        gamma1=v[2]
        beta2=v[3]
        gamma2=v[4]
        nu=v[5]
        beta1b=v[6]
        gamma1b=v[7]
        beta2b=v[8]
        gamma2b=v[9]
        nub=v[10]
        diff_e=function(x){abs(marg_prob(beta2,gamma2)(x)-marg_prob(beta2b,gamma2b)(x))}
        h <- rep(NA, 1)
        h[1]=beta1-beta1b
        h[2]=beta2-beta2b
        h[3]=(sum(diff_e(y)*exp(100*diff_e(y))))/(sum(exp(100*diff_e(y))))-epsilon
        h
      }
      minimum <- auglag(par=c(par_gumbel1,par_gumbel2),joint_likelihood,heq=softmax,hin=hin,control.outer=list(method="nlminb",trace=FALSE))$par
    }
    
    data1b <- list()
    data2b <- list()
    
      for (k in 1:B)
      {
        data1b <- gen_data(doses,minimum[1],minimum[2],minimum[3],minimum[4],minimum[5])
        data2b <- gen_data(doses2,minimum[6],minimum[7],minimum[8],minimum[9],minimum[10])
        
        #MLE
        sum <- 0
        likelihoodb <- function(v){
          beta1=v[1]
          gamma1=v[2]
          beta2=v[3]
          gamma2=v[4]
          nu=v[5]
          beta1b=v[6]
          gamma1b=v[7]
          beta2b=v[8]
          gamma2b=v[9]
          nub=v[10]
          for(i in 1:length(doses)){
            d <- doses[i]  
            sum <- sum + log(p11(d)(beta1,gamma1,beta2,gamma2,nu)^sum(data1b[[i]][,1]==1 & data1b[[i]][,2]==1)*p01(d)(beta1,gamma1,beta2,gamma2,nu)^sum(data1b[[i]][,1]==0 & data1b[[i]][,2]==1)*p10(d)(beta1,gamma1,beta2,gamma2,nu)^sum(data1b[[i]][,1]==1 & data1b[[i]][,2]==0)*p00(d)(beta1,gamma1,beta2,gamma2,nu)^sum(data1b[[i]][,1]==0 & data1b[[i]][,2]==0))
            + log(p11(d)(beta1,gamma1b,beta2,gamma2b,nub)^sum(data2b[[i]][,1]==1 & data2b[[i]][,2]==1)*p01(d)(beta1,gamma1b,beta2,gamma2b,nub)^sum(data2b[[i]][,1]==0 & data2b[[i]][,2]==1)*p10(d)(beta1,gamma1b,beta2,gamma2b,nub)^sum(data2b[[i]][,1]==1 & data2b[[i]][,2]==0)*p00(d)(beta1,gamma1b,beta2,gamma2b,nub)^sum(data2b[[i]][,1]==0 & data2b[[i]][,2]==0))
          }
          return(-sum) #for minimizing
        }
        par_gumbelb <- auglag(par=minimum,likelihoodb,heq=shared_placebos,hin=hin,control.outer=list(method="nlminb",trace=FALSE))$par
        par_gumbel1b <- par_gumbelb[1:5]
        par_gumbel2b <- par_gumbelb[6:10]
        
        
        d_e_star <- max(abs(marg_prob(par_gumbel1b[1],par_gumbel1b[2])(y)-marg_prob(par_gumbel2b[1],par_gumbel2b[2])(y))) #max.abs.distance efficacy
        d_t_star <- max(abs(marg_prob(par_gumbel1b[3],par_gumbel1b[4])(y)-marg_prob(par_gumbel2b[3],par_gumbel2b[4])(y))) #max. abs distance toxicity
        t.star <- d_t_star
    
      boot2[k] <- t.star
    }
  }
  crit_val2 <- quantile(boot2,0.05)
  pval2 <- ecdf(boot2)(t.stat2)
  
  critval <- c(critval,crit_val,crit_val2)
  pval <- c(pval,pval1,pval2)
}

rm(list=ls())

library(MASS)
library(glmnet)
library(tidyverse)

###################
# 1. read in data #
###################
load("rprotblds.RData")

#####################################################
# 2. simulation based power analysis
####################################################

# input parameters
set.seed(435) #seed for reproducibility
nsamp_true=704 #sample size
nrep=100 #number of replicates

# function to generate random draws from quasi-poisson
rqpois <- function(n, mu, theta) {
  rnbinom(n = n, mu = mu, size = mu/(theta-1))
}

# initialize data containers
simdata=list()
lassomod=list()
lassofit=list()
pow=list()
fdp=list()

# generate multivariate normal to be used to simulate 3,072 protein expression levels 
mus=rgamma(3072,shape=6.3,rate=1.2)
vars=rlnorm(3072,meanlog=-0.87,sdlog=0.76)
covs=rnorm((3072^2-3072)/2,mean=0.095,sd=0.10)

# find nearest positive definite covariance matrix
sigma=diag(vars)
sigma[lower.tri(sigma)]=covs
sigma[upper.tri(sigma)]=covs
sigma = (sigma + t(sigma))/2
sigma.evd=eigen(sigma)
sigma.psd = sigma.evd$vectors %*% diag(pmax(sigma.evd$values,0)) %*% t(sigma.evd$vectors)

# simulate data and fit LASSO for each replicate
for (r in 1:nrep){

  # bootstrap draw nsamp patients from pilot data with covariates
  nsamp=nsamp_true*2 #generating more samples since some draws are discarded due to having too high event rates
  df=rprotblds[sample(1:nrow(rprotblds),size=nsamp,replace=TRUE),c(1,2,96:99,115)]
  
  # simulate protein expression levels from multivariate normal
  df.prot=mvrnorm(n=nsamp,mu=mus,Sigma=sigma.psd)
  names(df.prot) = paste("protein",1:3072,sep="_")
  df.all=cbind(df,df.prot)
  colnames(df.all)[8:3079]= paste("protein",1:3072,sep="_")
  
  # simulate pseudo-responses for 6-month number of rehospitalizations or death
  respmus=numeric(length=nsamp)
  pseudoresp=numeric(length=nsamp)
  for (i in 1:nsamp){
    logmu=0.002*df.all$age[i]-0.006*(df.all$sex[i]==1)-0.41*(df.all$race___4[i]==1)+
      0.16*(df.all$hf_cat[i]==1)+ # covariates
      # treatment effect
      (df.all$intervention_1_control_0[i]==1)*(1-0.3*df.all$protein_1[i])+
      (df.all$intervention_1_control_0[i]==1)*(1-0.3*df.all$protein_2[i])+
      (df.all$intervention_1_control_0[i]==1)*(1-0.3*df.all$protein_3[i])+
      # protein effect under control only
      (df.all$intervention_1_control_0[i]==0)*(-3+0.62*df.all$protein_1[i])+
      (df.all$intervention_1_control_0[i]==0)*(-3+0.62*df.all$protein_2[i])+
      (df.all$intervention_1_control_0[i]==0)*(-3+0.62*df.all$protein_3[i])
  
    respmus[i]=exp(logmu)
    
    # generate pseudoresponses from quasi-poisson with 1.7 dispersion (matching pilot data)
    pseudoresp[i]=rqpois(1,respmus[i],1.7)
  }
  
  # drop observations that have an event rate above 10 (too high for 6 months)
  if(length(which(respmus>10))>0){
    df.final=df.all[-which(respmus>10),] 
    pseudoresp=pseudoresp[-which(respmus>10)]
    respmus=respmus[-which(respmus>10)]
  } else{
    df.final=df.all
  }
  
  # drop uneeded extra draws to get down to original sample size
  idx=sample(1:nrow(df.final),nsamp_true)
  df.final=df.final[idx,]
  respmus=respmus[idx]
  pseudoresp=pseudoresp[idx]
  
  # combine everything into a single dataframe
  pseudo_data <- cbind(pseudoresp,df.final[,c("intervention_1_control_0","age","sex","race___4","hf_cat")],
                         df.final[,8:3079])
  simdata[[r]] = pseudo_data
  
  # use LASSO to estimate model parameters
  f <- as.formula(pseudoresp ~ .+intervention_1_control_0*.)
  pseudo_y <- pseudo_data$pseudoresp
  pseudo_x <- model.matrix(f, pseudo_data)[, -1]

  fit <- cv.glmnet(pseudo_x, pseudo_y,family=quasipoisson())
  plot(fit)
  lassomod[[r]]=fit
  lassofit[[r]]=coef(fit,lamba=fit$lambda.min)[which(coef(fit,lambda=fit$lambda.min)!=0),]
  
  # get power and fdp for the run 
  # (correct identification of emm if both main effect and interaction terms are significant)
  pow[[r]] = ((sum(coef(fit,lamba=fit$lambda.min)[c(7,3083)]!=0)==2) +
    (sum(coef(fit,lamba=fit$lambda.min)[c(8,3084)]!=0)==2)+
      (sum(coef(fit,lamba=fit$lambda.min)[c(9,3085)]!=0)==2))/3
  
  tmpfdp=0
  for (i in 4:3072){
    tmpfdp = tmpfdp+(sum(coef(fit,lamba=fit$lambda.min)[c(6+i,3082+i)]!=0)==2)/3069
  }
  fdp[[r]]=tmpfdp
  
  print(paste("Replicate:",r,", Power:",pow[[r]],", FDP:",fdp[[r]]))
  
  save.image("PowerAnalysis_092823.RData")
  
  # plot single run pseudoresponse distribution and protein EMM behavior
  if (r==1){
    
    pdf("pseudoresponse_distribution.pdf", onefile = TRUE)
    
    # visualize distribution of pseudoresponse by treatment
    ggplot(data=data.frame(pseudoresp=pseudoresp,trt=df.final$intervention_1_control_0),
           mapping=aes(x=pseudoresp,fill=trt))+geom_histogram(aes(y = after_stat(count / sum(count))), 
                                                              position='identity',alpha=0.4)+
      labs(x="Pseudo Response (6-month # hospitalizations or death)",
           y="Proportion",title="Distribution of Pseudo Response by Treatment")
    
    # visualize treatment effect modification for EMM proteins (1,2,3) and non-EMM proteins (4,5)
    p=list()
    for (i in 1:5){
      prot=names(df.final)[1:15+7][i]
      p[[i]]=ggplot(data=df.final,mapping=aes(x=!!ensym(prot),y=pseudoresp,
                                              color=intervention_1_control_0))+geom_point()+
        geom_smooth(method='glm',method.args=list(family='quasipoisson'))+
        labs(title="Treatment Effect Modification")
      print(p[[i]])
    }
    dev.off()
  }

}

# output power and FDR
mean(unlist(pow))
mean(unlist(fdp))

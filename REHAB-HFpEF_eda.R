rm(list=ls())

library(tidyverse)
library(AER)
library(fitdistrplus)

###################
# 1. read in data #
###################
load("rprotblds.RData")

#####################################################
# 2. EDA for primary outcome (number of rehospitalizations and death events in 6 months)
####################################################

# check for overdispersion
mean(rprotblds$nrehosp_death) #1.300412
var(rprotblds$nrehosp_death) #2.268884
ggplot(data=rprotblds,mapping=aes(x=nrehosp_death))+geom_bar()
# variance is larger than mean, so potentially overdisperse

# test for overdispersion with covariates only (no proteins)
mod1=glm(data=rprotblds,formula=nrehosp_death ~ age + sex + race___4 + hf_cat
         + intervention_1_control_0,family='poisson')
dispersiontest(mod1)
# overdisperse conditionally with disperson of 1.682113

# reference for dispersion test:
#Cameron, A.C. and Trivedi, P.K. (1990). Regression-based Tests for Overdispersion
# in the Poisson Model. Journal of Econometrics, 46, 347–364.

# first build independent regression models for each protein individually
# and make volcano plot
intxn_df = matrix(NA,nrow=92,ncol=2)
colnames(intxn_df) = c("effect","p-value")
for (i in 3:94){
  prot=names(rprotblds)[i]
  formtmp=paste("nrehosp_death ~ age + sex + race___4 + hf_cat + `",prot,
             "`*intervention_1_control_0",sep="")
  tmpmod=glm(data=rprotblds,formula=formtmp,family='quasipoisson')
  intxn_df[i-2,]=coef(summary(tmpmod))[length(tmpmod$coefficients),c(1,4)]
}

intxn_df_all = data.frame(Protein=names(rprotblds)[3:94],intxn_df,
                      Significance=ifelse(intxn_df[,2]<0.05,"Significant","Not Significant"))

pdf("nrehospdeath_volcano_effectmod.pdf", onefile = TRUE)
ggplot(data=intxn_df_all,mapping=aes(x=effect,y=-log10(p.value),color=Significance))+geom_point()+
  geom_hline(yintercept = -log10(0.05), linetype="dotted")+
  labs(x = "Effect modification", y = "-log10(p-value)") +
  OlinkAnalyze::set_plot_theme() + ggtitle('NRD Effect Modification (Unadjusted P-values)')+
  ggrepel::geom_label_repel(data = intxn_df_all[intxn_df_all$Significance=="Significant",],
                            ggplot2::aes(label = Protein), box.padding = 1, show.legend = FALSE,max.overlaps=Inf) 
dev.off()

#PI3 only one that is significant (unadjusted for multiple testing), so visualize treatment effect modification
pdf("nrehospdeath_PI3_effectmod.pdf", onefile = TRUE)
ggplot(data=rprotblds,mapping=aes(x=`PI3`,y=nrehosp_death,color=intervention_1_control_0))+
  geom_point()+geom_smooth(method='glm',se=TRUE,method.args=list(family='quasipoisson'))+
  ylab("Number of Rehospitalizations and Death")+xlab("Baseline PI3")+
  scale_color_discrete(labels=c('Control', 'Intervention'))+
  theme(legend.title= element_blank())
dev.off()

#####################################################
# 3. Investigate distribution of protein expression levels
####################################################

pdf("proteinexpressiondistributions.pdf", onefile = TRUE)
# Use a multivariate normal distribution to simulate protein expression levels for 3,072 proteins.
# We will first approximate the distribution of the mean protein expression levels, as well as the variances
# and covariances of the protein expression levels using pilot data on 92 proteins.  

# Distribution of mean protein expression levels across proteins
summary(fitdist(as.numeric(apply(rprotblds[,3:94],2,function(x) mean(x,na.rm=TRUE))),distr='gamma',method='mme'))
plot(fitdist(as.numeric(apply(rprotblds[,3:94],2,function(x) mean(x,na.rm=TRUE))),distr='gamma',method='mme'))#draw means for protein expressions from normal(mean=5.2,sd=2.1)
#draw means from gamma with shape 6.3 and rate 1.2

#Distribution of variances of protein expression levels across proteins
summary(fitdist(as.numeric(apply(rprotblds[,3:94],2,function(x) var(x,na.rm=TRUE))),distr='lnorm',method='mme'))
plot(fitdist(as.numeric(apply(rprotblds[,3:94],2,function(x) var(x,na.rm=TRUE))),distr='lnorm',method='mme'))#draw means for protein expressions from normal(mean=5.2,sd=2.1)
#draw variances from log normal with meanlog -0.87 and sdlog=0.76

#Distribution of covariances of protein expression levels across proteins
covs=cov(rprotblds[,3:94])[lower.tri(cov(rprotblds[,3:94]))]
covs=covs[!is.na(covs)]
summary(fitdist(covs,distr='norm',method='mme'))
plot(fitdist(covs,distr='norm',method='mme'))#draw means for protein expressions from normal(mean=5.2,sd=2.1)
#draw covariances from normal mean=0.095, sd=0.10

dev.off()

#####################################################
# 4. Explore treatment and covariate effects in pilot data
####################################################

#Estimate treatment and covariate effects across protein-specific models
trtmt_df = matrix(NA,nrow=92,ncol=2)
age_df = matrix(NA,nrow=92,ncol=2)
sex_df = matrix(NA,nrow=92,ncol=2)
race_df = matrix(NA,nrow=92,ncol=2)
hfcat_df = matrix(NA,nrow=92,ncol=2)

for (i in 3:94){
  prot=names(rprotblds)[i]
  formtmp=paste("nrehosp_death ~ age + sex + race___4 + hf_cat + `",prot,
                "`*intervention_1_control_0",sep="")
  tmpmod=glm(data=rprotblds,formula=formtmp,family='quasipoisson')
  trtmt_df[i-2,]=coef(summary(tmpmod))[length(tmpmod$coefficients)-1,c(1,4)]
  age_df[i-2,]=coef(summary(tmpmod))[2,c(1,4)]
  sex_df[i-2,]=coef(summary(tmpmod))[3,c(1,4)]
  race_df[i-2,]=coef(summary(tmpmod))[4,c(1,4)]
  hfcat_df[i-2,]=coef(summary(tmpmod))[5,c(1,4)]
}

pdf("covariateeffects.pdf", onefile = TRUE)

ggplot(data=as.data.frame(trtmt_df),mapping=aes(x=V1))+geom_histogram()+
  labs(x="Treatment Effect")+geom_vline(xintercept=median(trtmt_df[,1]),color='red',linetype='dashed')
summary(trtmt_df[,1]) #median is 0.46

ggplot(data=as.data.frame(age_df),mapping=aes(x=V1))+geom_histogram()+
  labs(x="Age Effect")+geom_vline(xintercept=median(age_df[,1]),color='red',linetype='dashed')
summary(age_df[,1]) #median is 0.0002

ggplot(data=as.data.frame(sex_df),mapping=aes(x=V1))+geom_histogram()+
  labs(x="Sex Effect")+geom_vline(xintercept=median(sex_df[,1]),color='red',linetype='dashed')
summary(sex_df[,1]) #median is -0.06 (Sex=1 is female)

ggplot(data=as.data.frame(race_df),mapping=aes(x=V1))+geom_histogram()+
  labs(x="Race Effect")+geom_vline(xintercept=median(race_df[,1]),color='red',linetype='dashed')
summary(race_df[,1]) #median is -0.41 (Race=1 is white)

ggplot(data=as.data.frame(hfcat_df),mapping=aes(x=V1))+geom_histogram()+
  labs(x="Heart Failure Category Effect")+geom_vline(xintercept=median(hfcat_df[,1]),color='red',linetype='dashed')
summary(hfcat_df[,1]) #median is 0.16 (hfcat==1 is HFpEF) 

dev.off()
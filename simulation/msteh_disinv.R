################################################################################
# msteh_disinv.R
# Misspecification of Treatment Effect Heterogeneity - Disinvestment example
# Rhys Bowden
################################################################################
library("lme4")
library("readxl")
pld <- as.data.frame(read_excel("~/monash work/R code/teh/disinvestment/pmed.1002412.s002.xlsx")) # load disinvestment data
pld = pld[!duplicated(pld$episodenumber),] # remove duplicate entries
# important columns:
# outcome: acute_los
# treatment indicator: no_we_exposure
# time: sw_step (1:7 is disinvestment, 8:14 is reinvestment)
# cluster: index_ward
# hospital: site (1 has all wards and time periods, 2 does not have ward 5 in the second trial, or time period 14)
pld$log_acute_los = log(pld$acute_los) # outcome
trial = pld[pld$sw_step>7&pld$site==1,c('log_acute_los','no_we_exposure','sw_step','index_ward')] # select relevant columns
trial[,c(2,3,4,5)] = lapply(trial[,c(2,3,4)],as.factor)
names(trial) = c("Y","treatf","timef","clusterf")
trial[,"clustimef"] = interaction(trial[,"timef"],trial[,"clusterf"])

options(warn = 0)
set.seed(1)
design <- "sw" # sw is stepped-wedge, cxo is cluster crossover
if(design=="sw"){
  source("fgls5.R")
  source("tevpie.R")
}

# simulation parameters
individual_level = TRUE # whether to do the analysis at the individual level
include_nu = TRUE # whether to include a cluster-period random effect

converge_fail = 1
converge_fail2 = 1
var_theta_hat2 = 0 # variance of ideal, GLS estimator of theta (with known variance components)
est_var_theta1 = 0 # theoretical approximation to expected value of estimate of variance effect, estimated using misspecified model (model 1) 
hat_sigma2_a = 0 
hat_sigma2_eps = 0
hat_sigma2_a2 = 0
hat_sigma2_u2 = 0
hat_sigma_au2 = 0
hat_sigma2_eps2 = 0
hat_thetaMS = 0
hat_theta2 = 0
hat_var_theta_hat1 = 0 # lme4 estimate of the variance of the treatment effect, assuming Model 1 ("misspecified"). 
hat_var_theta_hat2 = 0 # lme4 estimate of the variance of the treatment effect, assuming Model 2 ("correctly specified"). 
if(include_nu){
  hat_sigma2_nu = 0
  hat_sigma2_nu2 = 0
}
TT = 7 # number of time periods
SS = TT-1 # number of treatment sequences
cc = 1 # number of clusters per sequence
CC = cc*SS # total number of clusters
ni = nrow(trial)/(TT*CC) # average number of observations per cluster period

# ---- inference ----
# model 1 (misspecified)
output = tryCatch({
  if(!include_nu){
    fitY = lmer(Y ~ timef + treatf + (1 | clusterf),data=trial)
  }else{
    fitY = lmer(Y ~ timef + treatf + (1 | clusterf) + (1 | clustimef),data=trial)
  }
  list(fitY,2)
},
warning=function(war){
  print(paste("MY_WARNING:  ",war))
  return(list(fitY,1))
})
fitY = output[[1]]
converge_fail=output[[2]]
# model 2 (correctly specified)
output = tryCatch({
  if(!include_nu){
    fitY2 = lmer(Y ~ timef + treatf + (treatf | clusterf),data=trial)
  }else{
    fitY2 = lmer(Y ~ timef + treatf + (treatf | clusterf) + (1 | clustimef),data=trial)
  }
  list(fitY2,2)
},
warning=function(war){
  #converge_fail2=1
  print(paste("MY_WARNING:  ",war))
  if(exists("fitY2")){
    return(list(fitY2,1))
  }else{
    return(list(-1,1))
  }
},
error=function(err){
  stop(paste("MY_ERROR: ",err))
  return(list(-1,1))
})
fitY2 = output[[1]]
converge_fail2=output[[2]]

# ---- fixed effects estimates ----
hat_thetaMS = fixef(fitY)["treatf1"]
hat_theta2 = fixef(fitY2)["treatf1"]

# ---- variance component estimates (hat_sigma,hat_var_theta_hat) ----
hat_sigmas = as.data.frame(VarCorr(fitY))
hat_sigma2_a = hat_sigmas[hat_sigmas$grp=='clusterf','vcov']
hat_sigma2_eps = hat_sigmas[hat_sigmas$grp=='Residual','vcov']
hat_sigmas2 = as.data.frame(VarCorr(fitY2))
hat_sigma2_a2 = hat_sigmas2[hat_sigmas2$grp=='clusterf'&hat_sigmas2$var1=='(Intercept)'&is.na(hat_sigmas2$var2),'vcov']
hat_sigma2_u2 = hat_sigmas2[hat_sigmas2$grp=='clusterf'&hat_sigmas2$var1=='treatf1'&is.na(hat_sigmas2$var2),'vcov']
hat_sigma_au2 = hat_sigmas2[hat_sigmas2$grp=='clusterf'&hat_sigmas2$var1=='(Intercept)'&!is.na(hat_sigmas2$var2)&hat_sigmas2$var2=='treatf1','vcov']
hat_sigma2_eps2 = hat_sigmas2[hat_sigmas2$grp=='Residual','vcov']
if(include_nu){
  hat_sigma2_nu = hat_sigmas[hat_sigmas$grp=='clustimef','vcov']
  hat_sigma2_nu2 = hat_sigmas2[hat_sigmas2$grp=='clustimef','vcov']
}

hat_var_theta_hat1 = (summary(fitY)$coef[TT+1, 2])^2 # square of the std error of the estimate for treatment effect
if(typeof(fitY2)=="double"){
  my_error_log2=1
} else {
  hat_var_theta_hat2 = (summary(fitY2)$coef[TT+1, 2])^2 # square of the std error of the estimate for treatment effect
}

# calculation of var ratio
# everything is done at cluster mean level
cluster_labels = rep(1:CC,each=TT)
period_labels = rep(1:TT,CC)
Ntot = CC*TT
B = 0+(matrix(rep(cluster_labels,each=Ntot),ncol=Ntot)==matrix(rep(cluster_labels,Ntot),ncol=Ntot))
treat_grid_s = matrix(0,nrow=SS,ncol=TT)
treat_grid_s[upper.tri(treat_grid_s,diag=FALSE)]=1
treat_grid_c = kronecker(treat_grid_s,matrix(1,nrow=cc,ncol=1))
treat_v = matrix(t(treat_grid_c),ncol=1)
Bu = B*matrix(rep(treat_v,each=Ntot),ncol=Ntot)*matrix(rep(treat_v,Ntot),ncol=Ntot) # Bu(i,j) = 1 if same cluster and both treated
Buc = B*(matrix(rep(treat_v,each=Ntot),ncol=Ntot)|matrix(rep(treat_v,Ntot),ncol=Ntot)) # Buc(i,j) = 1 if same cluster and either treated
XT = kronecker(matrix(1,nrow=CC,ncol=1),diag(TT))
X = cbind(treat_v,XT)

sigma2A=hat_sigma2_a2
sigma2U=hat_sigma2_u2
sigmaAU=hat_sigma_au2
sigma2E=hat_sigma2_eps2
sigma2NU=hat_sigma2_nu2
sigma2G= sigma2E/ni+sigma2NU # cluster mean level variance
# long versions of parameters:
Rkls = fgls5(TT,cc,ni)
Rkls = Rkls[,,1,1]
sigma2_hats = Rkls%*%matrix(c(sigma2G,sigma2A,sigma2U,sigmaAU),ncol=1)
est_var_theta1 = tevpie(sigma2_hats[1,1],sigma2_hats[2,1],TT,cc)

true_covar = diag(CC*TT)*sigma2G+B*sigma2A+Bu*(sigma2U+sigmaAU)+Buc*sigmaAU# equivalent variance of estimator if using cell means,vars vs full expanded model
true_covar_eta = solve(t(X)%*%solve(true_covar)%*%X) # variance matrix for fixed effects
var_theta_hat2 = true_covar_eta[1,1] # variance of estimate of theta under model 2 with known correlations

lme_var_ratio = hat_var_theta_hat1/hat_var_theta_hat2 # misspecification error ratio for lme4
var_ratio = est_var_theta1/var_theta_hat2 # misspecification error ratio using theoretical approximation and ideal variance of GLS estimator

# estimates of cluster level correlations
rhoc  = hat_sigma2_a/(hat_sigma2_a+hat_sigma2_nu+hat_sigma2_eps)
rhocc = hat_sigma2_a2/(hat_sigma2_a2+hat_sigma2_nu2+hat_sigma2_eps2)
rhott = (hat_sigma2_a2+2*hat_sigma_au2+hat_sigma2_u2)/(hat_sigma2_a2+2*hat_sigma_au2+hat_sigma2_u2+hat_sigma2_nu2+hat_sigma2_eps2)
rhoct = (hat_sigma2_a2+hat_sigma_au2)/sqrt((hat_sigma2_a2+hat_sigma2_nu2+hat_sigma2_eps2)*(hat_sigma2_a2+2*hat_sigma_au2+hat_sigma2_u2+hat_sigma2_nu2+hat_sigma2_eps2))
# estimates of cluster-period level correlations
rhocp = (hat_sigma2_nu+hat_sigma2_a)/(hat_sigma2_a+hat_sigma2_nu+hat_sigma2_eps)
rhocpc = (hat_sigma2_nu+hat_sigma2_a)/(hat_sigma2_a+hat_sigma2_nu+hat_sigma2_eps)
rhocpt = (hat_sigma2_a2+2*hat_sigma_au2+hat_sigma2_u2+hat_sigma2_nu2)/(hat_sigma2_a2+2*hat_sigma_au2+hat_sigma2_u2+hat_sigma2_nu2+hat_sigma2_eps2)

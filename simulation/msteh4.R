#####################################################
# Misspecification of Treatment Effect Heterogeneity
# Rhys Bowden
#####################################################
# with thanks to Ren et al. 2019 A Simulation Study of Statistical Approaches to Data for the form of the calls to lme4.
# v4: paper version
# Analysis of the Stepped Wedge Design
library("lme4")
library("tictoc")
library("abind")
options(warn = 0)
source("rho_to_sigma_ratio2.R")
set.seed(1)
design <- "sw" # sw is stepped-wedge, cxo is cluster crossover
if(design=="sw"){
  source("fgls5.R")
  source("tevpie.R")
}

# simulation parameters
num_repeats = 20000
individual_level = TRUE
include_nu = TRUE
store_full_data = FALSE
output_text_file = FALSE # requires store_full_data to be TRUE
plot_results = FALSE

# design and model parameters
T_list = 5 # number of periods
c_list = 5 # number of clusters per period
N_list = 20 # number of individual observations per cluster-period
rhocc_list = 0.01
rhott_list = seq(0.01,0.05,by=0.002)
rhoct_list = 0.01 
rhoac_list = 0.8 # cluster auto-correlation
sigma_eps2 = 1 # error variance. The other variance parameters are determined by this and the correlation values above

if(design=="cxo"){
  T_list = T_list-T_list%%2 # ensure T is even
}

# generate each combination of parameters
param_grid = expand.grid(T_list,c_list,N_list,rhocc_list,rhott_list,rhoct_list,rhoac_list)
names(param_grid) = c('T','c','ni','rhocc_long','rhott_long','rhoct_long','rhoac_long')
num_param_tuples = dim(param_grid)[1]

# determine variances from correlation parameters
sigrats = rho_to_sigma_ratio2(param_grid$rhocc_long,param_grid$rhott_long,param_grid$rhoct_long,param_grid$rhoac_long)
param_grid$sigmaE2 = rep_len(sigma_eps2,num_param_tuples)
param_grid$sigmaA2 = sigrats[[1]]*param_grid$sigmaE2
param_grid$sigmaU2 = sigrats[[2]]*param_grid$sigmaE2
param_grid$sigmaAU = sigrats[[3]]*param_grid$sigmaE2
param_grid$sigmaNU2 = sigrats[[4]]*param_grid$sigmaE2

if(design == "cxo"){
  Cvec = 2*param_grid[['c']]
  Tvec = param_grid[['T']]
  Nvec = param_grid[['ni']]
  sE2vec = param_grid[['sigmaE2']] # sigma_epsilon^2
  sNU2vec = param_grid[['sigmaNU2']] # sigma_nu^2
  sG2vec = sE2vec/Nvec+sigmaNU2list # sigma_gamma^2
  sU2vec = param_grid[['sigmaU2']] # sigma_U^2
  th_correct_var = (4*sG2vec/Tvec+sU2vec)/Cvec
  th_misspec_var = (4*sG2vec/Tvec+((Cvec-1)/(Cvec*Tvec-Cvec-Tvec))*sU2vec)/Cvec
}else if(design=="sw"){
  th_misspec_var = vector(mode = "list", length = num_param_tuples)
  for(param_ind in 1:num_param_tuples){
    # determining theoretical plug-in estimated value of treatment effect variance
    TT = param_grid[param_ind,"T"]
    cc = param_grid[param_ind,"c"]
    ni = param_grid[param_ind,"ni"]
    Rkls = fgls5(TT,cc,ni)
    Rkls = Rkls[,,1,1]
    sG2 = param_grid[param_ind,"sigmaE2"]/ni+param_grid[param_ind,"sigmaNU2"]
    sigma2_hats = Rkls%*%matrix(c(sG2,param_grid[param_ind,"sigmaA2"],param_grid[param_ind,"sigmaU2"],param_grid[param_ind,"sigmaAU"]),ncol=1)
    est_var = tevpie(sigma2_hats[1,1],sigma2_hats[2,1],TT,cc)
    th_misspec_var[param_ind] = est_var
  }
}

# results variables----
fitYs = array(0,dim=c(num_param_tuples,num_repeats))
fitY2s = array(0,dim=c(num_param_tuples,num_repeats))
bar_sigma2_a1 = array(0,dim=c(num_param_tuples))
bar_sigma2_eps1 = array(0,dim=c(num_param_tuples))
bar_sigma2_a = array(0,dim=c(num_param_tuples))
bar_sigma2_eps = array(0,dim=c(num_param_tuples))
bar_sigma2_a2 = array(0,dim=c(num_param_tuples))
bar_sigma2_eps2 = array(0,dim=c(num_param_tuples))
bar_sigma2_u2 = array(0,dim=c(num_param_tuples))
bar_sigma_au2 = array(0,dim=c(num_param_tuples))
if(include_nu){
  bar_sigma2_nu1 = array(0,dim=c(num_param_tuples))
  bar_sigma2_nu = array(0,dim=c(num_param_tuples))
  bar_sigma2_nu2 = array(0,dim=c(num_param_tuples))
}
bar_hat_var_theta_hat1 = array(0,dim=c(num_param_tuples))
bar_hat_var_theta_hatMS = array(0,dim=c(num_param_tuples))
bar_hat_var_theta_hat2 = array(0,dim=c(num_param_tuples))
bar_theta_hat1 = array(0,dim=c(num_param_tuples))
bar_theta_hatMS = array(0,dim=c(num_param_tuples))
bar_theta_hat2 = array(0,dim=c(num_param_tuples))
emp_var_theta_hat1 = array(0,dim=c(num_param_tuples))
emp_var_theta_hatMS = array(0,dim=c(num_param_tuples))
emp_var_theta_hat2 = array(0,dim=c(num_param_tuples))
bar_converge = array(0,dim=c(num_param_tuples))
bar_converge2 = array(0,dim=c(num_param_tuples))
bar_error2 = array(0,dim=c(num_param_tuples))
hat_theta1 = array(0,dim=c(num_param_tuples,num_repeats))
hat_thetaMS = array(0,dim=c(num_param_tuples,num_repeats))
hat_theta2 = array(0,dim=c(num_param_tuples,num_repeats))
hat_sigma2_a1 = array(0,dim=c(num_param_tuples,num_repeats))
hat_sigma2_eps1 = array(0,dim=c(num_param_tuples,num_repeats))
hat_sigma2_a = array(0,dim=c(num_param_tuples,num_repeats))
hat_sigma2_eps = array(0,dim=c(num_param_tuples,num_repeats))
hat_sigma2_a2 = array(0,dim=c(num_param_tuples,num_repeats))
hat_sigma2_eps2 = array(0,dim=c(num_param_tuples,num_repeats))
hat_sigma2_u2 = array(0,dim=c(num_param_tuples,num_repeats))
hat_sigma_au2 = array(0,dim=c(num_param_tuples,num_repeats))
hat_sigma2_nu1 = array(0,dim=c(num_param_tuples,num_repeats))
hat_sigma2_nu = array(0,dim=c(num_param_tuples,num_repeats))
hat_sigma2_nu2 = array(0,dim=c(num_param_tuples,num_repeats))
converge_fail1 = array(0,dim=c(num_param_tuples,num_repeats))
converge_fail = array(0,dim=c(num_param_tuples,num_repeats))
converge_fail2 = array(0,dim=c(num_param_tuples,num_repeats))
my_error_log = array(0,dim=c(num_param_tuples,num_repeats))
my_error_log2 = array(0,dim=c(num_param_tuples,num_repeats))
hat_var_theta_hat1 = array(0,dim=c(num_param_tuples,num_repeats))
hat_var_theta_hatMS = array(0,dim=c(num_param_tuples,num_repeats))
hat_var_theta_hat2 = array(0,dim=c(num_param_tuples,num_repeats))
if(design=="cxo"){
  bar_hat_theta_mean_diff = array(0,dim=c(num_param_tuples))
  var_hat_theta_mean_diff = array(0,dim=c(num_param_tuples))
  hat_theta_mean_diff = array(0,dim=c(num_param_tuples,num_repeats))
}
# ----
if(store_full_data){
  YY1ifull = vector("list",length=num_param_tuples)
  YY1full = vector("list",length=num_param_tuples)
  YYifull = vector("list",length=num_param_tuples)
  YYfull = vector("list",length=num_param_tuples)
}
tic("outer loop")
for (param_ind in 1:num_param_tuples) {
  # design parameters ----
  TT = param_grid[param_ind,'T'] # number of time periods
  cc = param_grid[param_ind,'c'] # number of clusters per treatment sequence
  ni = param_grid[param_ind,'ni'] # number of individuals per cluster
  if(design=="sw"){
    S = TT-1 # number of treatment sequences
    treat_mat_s <- matrix(0,S+1,TT) # treatment matrix: treatment sequences are rows, times are columns
    treat_mat_s[upper.tri(treat_mat_s)] <- 1
    treat_mat_s <- treat_mat_s[1:S,1:TT]
    treatM <- treat_mat_s[rep(seq_len(S),each=cc),]
  }
  if(design=="cxo"){
    S = 2 # number of treatment sequences
    treat_mat_s = matrix(c(1,0,0,1),2,2)
    treatM = treat_mat_s %x% matrix(1,nrow=cc,ncol=(TT/2))
  }
  C <- S*cc # total number of clusters
  if(store_full_data){
    YY1ifull[[param_ind]] = matrix(nrow=num_repeats,ncol=ni*C*TT)
    YY1full[[param_ind]] = matrix(nrow=num_repeats,ncol=ni*C*TT)
    YYifull[[param_ind]] = matrix(nrow=num_repeats,ncol=ni*C*TT)
    YYfull[[param_ind]] = matrix(nrow=num_repeats,ncol=ni*C*TT)
  }
  
  # model parameters ----
  sigmaE = sqrt(param_grid[param_ind,'sigmaE2']) # sd of epsilon (errors)
  sigmaA = sqrt(param_grid[param_ind,'sigmaA2']) # sd of a (cluster effects)
  sigmaU = sqrt(param_grid[param_ind,'sigmaU2']) # sd of U (TEH in a cluster)
  sigmaNU = sqrt(param_grid[param_ind,'sigmaNU2']) # sd of nu (cluster-period random effect)
  sigmaAU = param_grid[param_ind,'sigmaAU']#covariance of a and U
  beta <- numeric(TT) # vector of fixed time effects
  beta <- c(1,-1,1,0,-1)
  theta <- 1
  sigma1 <- sigmaAU/sigmaA
  sigma2 <- sqrt(sigmaU^2-sigma1^2)
  
  # factors ----
  timeM <- t(matrix(1:TT,TT,C))
  clusterM <- matrix(1:C,C,TT)
  clustimeM <- matrix(1:(C*TT),C,TT)
  timeV <- timeM
  dim(timeV)<- c(C*TT,1)
  clusterV <- clusterM 
  dim(clusterV) <- c(C*TT,1)
  clustimeV <- clustimeM
  dim(clustimeV) <- c(C*TT,1)
  treatV <- treatM
  dim(treatV)<-c(C*TT,1)
  
  treatVi <- rep(treatV,each=ni)
  clusterVi<- rep(clusterV,each=ni)
  timeVi <- rep(timeV,each=ni)
  clustimeVi <- rep(clustimeV,each=ni)
  if(!individual_level){
    treatf <- factor(treatV)
    timef <- factor(timeV)
    clusterf <- factor(clusterV)
    clustimef <- interaction(clusterf,timef)
  }else{
    treatf <- factor(treatVi)
    timef <- factor(timeVi)
    clusterf <- factor(clusterVi)
    clustimef <- interaction(clusterf,timef)
  }
  
  # iteration ---- 
  for(i in 1:num_repeats){
    # ---- sample ----
    epsi <- rnorm(C*TT*ni,mean=0,sd=sigmaE) # epsilon (error term)
    eps <- epsi
    dim(eps) <- c(ni,C*TT)
    eps <- colSums(eps)/ni
    za <- rnorm(C)
    zu <- rnorm(C)
    a <- za*sigmaA # a (cluster random effects)
    aV <- rep(a,TT)
    u <- (za*sigma1+zu*sigma2) # u (treatment effect heterogeneity random effects)
    uV <- rep(u,TT)*treatV
    nu <- sigmaNU*rnorm(C*TT)
    nuV <- nu
    
    betaM <- t(matrix(beta,TT,C))
    betaV <- betaM
    dim(betaV) <- c(TT*C,1)
    thetaV <- theta*treatV
    thetaVi <- theta*treatVi
    
    aVi <- rep(aV,each=ni)
    uVi <- rep(uV,each=ni)
    nuVi <- rep(nuV,each=ni)
    betaVi <- rep(betaV,each=ni)
    
    YY1 =  eps +aV     +betaV +thetaV +nuV  # model 1: no TEH
    YY1i = epsi+aVi    +betaVi+thetaVi+nuVi # individual level version of  model 1
    YY =   eps +aV +uV +betaV +thetaV +nuV  # model 2: TEH
    YYi =  epsi+aVi+uVi+betaVi+thetaVi+nuVi # individuals version of model 2
    # ----- data allocation and factors -----
    if(individual_level){
      Y1 = YY1i
      Y = YYi
      Ydf1 = data.frame(YY1i,treatf,timef,clusterf,clustimef)
      Ydf = data.frame(YYi,treatf,timef,clusterf,clustimef)
    }else{
      Y1 = YY1
      Y = YY
      Ydf1 = data.frame(YY1,treatf,timef,clusterf,clustimef)
      Ydf = data.frame(YY,treatf,timef,clusterf,clustimef)
    }
    # make the data column name the same
    names(Ydf1)[1] = "Y1"
    names(Ydf)[1] = "Y"
    # ----- store data -----
    if(store_full_data){
      YY1ifull[[param_ind]][i,] = YY1i
      YY1full[[param_ind]][i,] = YY1
      YYifull[[param_ind]][i,] = YYi
      YYfull[[param_ind]][i,] = YY
    }
    
    # ---- inference ----
    # model 1 (correctly specified)
    output = tryCatch({
      if(!include_nu){
        fitY1 = lmer(Y1 ~ timef + treatf + (1 | clusterf),data=Ydf1)
      }else{
        fitY1 = lmer(Y1 ~ timef + treatf + (1 | clusterf) + (1 | clustimef),data=Ydf1)
      }
      list(fitY1,2)
    },
    warning=function(war){
      print(paste("MY_WARNING:  ",war))
      return(list(fitY1,1))
    })
    fitY1 = output[[1]]
    converge_fail1[param_ind,i]=output[[2]]
    
    # model 1 (misspecified)
    output = tryCatch({
      if(!include_nu){
        fitY = lmer(Y ~ timef + treatf + (1 | clusterf),data=Ydf)
      }else{
        fitY = lmer(Y ~ timef + treatf + (1 | clusterf) + (1 | clustimef),data=Ydf)
      }
      list(fitY,2)
    },
    warning=function(war){
      #converge_fail[i]=1
      print(paste("MY_WARNING:  ",war))
      return(list(fitY,1))
    })
    fitY = output[[1]]
    converge_fail[param_ind,i]=output[[2]]
    # model 2 (correctly specified)
    output = tryCatch({
      if(!include_nu){
        fitY2 = lmer(Y ~ timef + treatf + (treatf | clusterf),data=Ydf)
      }else{
        fitY2 = lmer(Y ~ timef + treatf + (treatf | clusterf) + (1 | clustimef),data=Ydf)
      }
      list(fitY2,2)
    },
    warning=function(war){
      #converge_fail2[i]=1
      #print(paste("MY_WARNING:  ",war))
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
    converge_fail2[param_ind,i]=output[[2]]
    
    # ---- fixed effects estimates ----
    hat_theta1[param_ind,i] = fixef(fitY1)["treatf1"]
    hat_thetaMS[param_ind,i] = fixef(fitY)["treatf1"]
    hat_theta2[param_ind,i] = fixef(fitY2)["treatf1"]
    if(design=="cxo"){
      hat_theta_mean_diff[param_ind,i] = (2/(C*TT))*((2*t(treatV)-1)%*%YY)
    }
    
    # ---- variance component estimates (hat_sigma,hat_var_theta_hat) ----
    hat_sigmas1 = as.data.frame(VarCorr(fitY1))
    hat_sigma2_a1[param_ind,i] = hat_sigmas1[hat_sigmas1$grp=='clusterf','vcov']
    hat_sigma2_eps1[param_ind,i] = hat_sigmas1[hat_sigmas1$grp=='Residual','vcov']
    hat_sigmas = as.data.frame(VarCorr(fitY))
    hat_sigma2_a[param_ind,i] = hat_sigmas[hat_sigmas$grp=='clusterf','vcov']
    hat_sigma2_eps[param_ind,i] = hat_sigmas[hat_sigmas$grp=='Residual','vcov']
    hat_sigmas2 = as.data.frame(VarCorr(fitY2))
    hat_sigma2_a2[param_ind,i] = hat_sigmas2[hat_sigmas2$grp=='clusterf'&hat_sigmas2$var1=='(Intercept)'&is.na(hat_sigmas2$var2),'vcov']
    hat_sigma2_u2[param_ind,i] = hat_sigmas2[hat_sigmas2$grp=='clusterf'&hat_sigmas2$var1=='treatf1'&is.na(hat_sigmas2$var2),'vcov']
    hat_sigma_au2[param_ind,i] = hat_sigmas2[hat_sigmas2$grp=='clusterf'&hat_sigmas2$var1=='(Intercept)'&!is.na(hat_sigmas2$var2)&hat_sigmas2$var2=='treatf1','vcov']
    hat_sigma2_eps2[param_ind,i] = hat_sigmas2[hat_sigmas2$grp=='Residual','vcov']
    if(include_nu){
      hat_sigma2_nu1[param_ind,i] = hat_sigmas1[hat_sigmas1$grp=='clustimef','vcov']
      hat_sigma2_nu[param_ind,i] = hat_sigmas[hat_sigmas$grp=='clustimef','vcov']
      hat_sigma2_nu2[param_ind,i] = hat_sigmas2[hat_sigmas2$grp=='clustimef','vcov']
    }
    
    hat_var_theta_hat1[param_ind,i] = (summary(fitY1)$coef[TT+1, 2])^2
    hat_var_theta_hatMS[param_ind,i] = (summary(fitY)$coef[TT+1, 2])^2 # square of the std error of the estimate for treatment effect
    if(typeof(fitY2)=="double"){
      my_error_log2[param_ind,i]=1
    } else {
      hat_var_theta_hat2[param_ind,i] = (summary(fitY2)$coef[TT+1, 2])^2 # square of the std error of the estimate for treatment effect
    }
    if(param_ind==1 & output_text_file){
      cat(sprintf("Data %d\n",i),file='teh/fit data/fitY2_param11.txt',append=TRUE)
      sink('teh/fit data/fitY2_param11.txt',append=TRUE)
      print(summary(fitY2))
      print(VarCorr(fitY2))
      write("\n--------------------------------------\n",file='teh/fit data/fitY2_param11.txt',append=TRUE)
      sink()
    }
  } # end of reps
  # ---- means for each parameter set ----
  bar_theta_hat1[param_ind] = mean(hat_theta1[param_ind,])
  bar_theta_hatMS[param_ind] = mean(hat_thetaMS[param_ind,])
  bar_theta_hat2[param_ind] = mean(hat_theta2[param_ind,])
  emp_var_theta_hat1[param_ind] = var(hat_theta1[param_ind,])
  emp_var_theta_hatMS[param_ind] = var(hat_thetaMS[param_ind,])
  emp_var_theta_hat2[param_ind] = var(hat_theta2[param_ind,])
  if(design=="cxo"){
    bar_hat_theta_mean_diff[param_ind] = mean(hat_theta_mean_diff[param_ind,])
    var_hat_theta_mean_diff[param_ind] = var(hat_theta_mean_diff[param_ind,])
  }
  bar_sigma2_a1[param_ind] = mean(hat_sigma2_a1[param_ind,])
  bar_sigma2_eps1[param_ind] = mean(hat_sigma2_eps1[param_ind,])
  bar_sigma2_a[param_ind] = mean(hat_sigma2_a[param_ind,])
  bar_sigma2_eps[param_ind] = mean(hat_sigma2_eps[param_ind,])
  bar_sigma2_a2[param_ind] = mean(hat_sigma2_a2[param_ind,])
  bar_sigma2_eps2[param_ind] = mean(hat_sigma2_eps2[param_ind,])
  bar_sigma2_u2[param_ind] = mean(hat_sigma2_u2[param_ind,])
  bar_sigma_au2[param_ind] = mean(hat_sigma_au2[param_ind,])
  if(include_nu){
    bar_sigma2_nu1[param_ind] = mean(hat_sigma2_nu1[param_ind,])
    bar_sigma2_nu[param_ind] = mean(hat_sigma2_nu[param_ind,])
    bar_sigma2_nu2[param_ind] = mean(hat_sigma2_nu2[param_ind,])
  }
  bar_hat_var_theta_hat1[param_ind] = mean(hat_var_theta_hat1[param_ind,])
  bar_hat_var_theta_hatMS[param_ind] = mean(hat_var_theta_hatMS[param_ind,])
  bar_hat_var_theta_hat2[param_ind] = mean(hat_var_theta_hat2[param_ind,my_error_log2[param_ind,]==0])
  bar_converge[param_ind] = mean(2-converge_fail[param_ind,])
  bar_converge2[param_ind] = mean(2-converge_fail2[param_ind,])
  bar_error2[param_ind] = mean(my_error_log2[param_ind,])
}   # end of parameter choice loop
total_time = toc()
if(include_nu){
  bar_sigma2_df = data.frame('sigA1'=bar_sigma2_a1,'sigE1'=bar_sigma2_eps1,'sigNU1'=bar_sigma2_nu1,'sigA'=bar_sigma2_a,'sigE'=bar_sigma2_eps,'sigNU'=bar_sigma2_nu,'sigA2'=bar_sigma2_a2,'sigU2'=bar_sigma2_u2,'sigAU'=bar_sigma_au2,'sigE2'=bar_sigma2_eps2,'sigNU1'=bar_sigma2_nu2)
}else{
  bar_sigma2_df = data.frame('sigA1'=bar_sigma2_a1,'sigE1'=bar_sigma2_eps1,'sigA'=bar_sigma2_a,'sigE'=bar_sigma2_eps,'sigA2'=bar_sigma2_a2,'sigU2'=bar_sigma2_u2,'sigAU'=bar_sigma_au2,'sigE2'=bar_sigma2_eps2)
}
# put all the variance parameters together in one array
sig_names = c('sigA1','sigE1','sigA','sigE','sigA2','sigU2','sigAU','sigE2')
hat_sigma2_array = abind(hat_sigma2_a1,hat_sigma2_eps1,hat_sigma2_a,hat_sigma2_eps,hat_sigma2_a2,hat_sigma2_u2,hat_sigma_au2,hat_sigma2_eps2,new.names=sig_names,along=3) # combine all estimates of variance components into one array
sig_zeroes = apply(hat_sigma2_array==0,c(1,3),mean) # proportion of estimates of each variance parameter that are equal to 0, for each parameter set

if(plot_results){
  rhott_long = param_grid$rhott_long
  # comparison shown in the paper - plug-in estimate of variance of treatment effect estimator compared with estimated variance of treatment effect estimator given by simulation, vs rho_TT 
  plot(rhott_long,bar_hat_var_theta_hatMS,ylim=range(c(th_misspec_var,bar_hat_var_theta_hatMS)),xlab=expression(paste(rho[TT])),ylab=expression(paste('var(',theta[1],')')))#ylim=c(min(th_misspec_rho0),1)
  lines(rhott_long,th_misspec_var,col="red",pch=4)
  # for cxo trials: comparison of the empirical variance of the treatment effect estimator, the theoretical variance of the estimator given in the paper, and the empirical variance of the "mean difference between treated and untreated cells" estimator
  plot(rhott_long,emp_var_theta_hatMS,ylim=range(c(th_misspec_var,emp_var_theta_hatMS)),xlab=expression(paste(rho[TT])),ylab=expression(paste('var2(',theta[1],')')))#ylim=c(min(th_misspec_rho0),1)
  lines(rhott_long,th_misspec_var,col="red",pch=4)
  lines(rhott_long,var_hat_theta_mean_diff,col="blue",pch=4)
}
if(output_text_file){
  Y = YYifull[[1]]
  factors = data.frame(treat=treatVi,time=timeVi,cluster=clusterVi,clustime = clustimeVi)
  write.table(cbind(factors,as.data.frame(t(Y))),"Y_11a.txt",row.names = FALSE)
  write.table(factors,"factors_cols.txt")
  write.table(as.integer(hat_sigma2_a2[1,]==0),"hat_sigma2_a2_is0_param11a.txt")
  write.table(param_grid,"parameters.txt")
}

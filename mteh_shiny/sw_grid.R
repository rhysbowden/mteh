#library('matlib')
source('rho_to_sigma_ratio2.R')
source('fgls5.R')
source('tevpie.R')
sw_grid = function(T_list,rhott_list,panel_var_list,fixed_var,rho_cc=0.01,rho_ct=0.01,rho_ac=0.8,panel_var=1){
  # hard-coded parameters and options
  small_sample_correction = 1
  # calculate sigmas from rhos
  sigma2E = 1 # individual error rate
  rhocc_list = rep(rho_cc,length(rhott_list))
  rhoct_list = rep(rho_ct,length(rhott_list))
  rhoac_list = rep(rho_ac,length(rhott_list))
  sigma_scenarios_list = rho_to_sigma_ratio2(rhocc_list,rhott_list,rhoct_list,rhoac_list)
  sigma2Alist = sigma2E*sigma_scenarios_list[[1]]
  sigma2Ulist = sigma2E*sigma_scenarios_list[[2]]
  sigmaAUlist = sigma2E*sigma_scenarios_list[[3]]
  sigma2NUlist = sigma2E*sigma_scenarios_list[[4]]
  # initialise
  long_ind = 1
  num_iters = length(rhott_list)*length(T_list)*length(panel_var_list)
  rhott_long = matrix(0,ncol=1,nrow=num_iters)
  TT_long = matrix(0,ncol=1,nrow=num_iters)
  cc_long = matrix(0,ncol=1,nrow=num_iters)
  ni_long = matrix(0,ncol=1,nrow=num_iters)
  est_var_long = matrix(0,ncol=1,nrow=num_iters)
  true_var_long = matrix(0,ncol=1,nrow=num_iters)
  sigma2U_long = matrix(0,ncol=1,nrow=num_iters)
  sigma2G_long = matrix(0,ncol=1,nrow=num_iters)
  for(T_ind in 1:length(T_list)){
    TT = T_list[T_ind]
    SS = TT-1
    for(panel_ind in 1:length(panel_var_list)){
      if(panel_var==1){
        cc = panel_var_list[panel_ind]
        ni = fixed_var
      }else{
        ni = panel_var_list[panel_ind]
        cc = fixed_var
      }
      CC = cc*SS
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
      
      for(scen_ind in 1:length(rhott_list)){
        cat('(',T_ind,panel_ind,scen_ind,')')
        sigma2A=sigma2Alist[scen_ind]
        sigma2U=sigma2Ulist[scen_ind]
        sigmaAU=sigmaAUlist[scen_ind]
        sigma2NU=sigma2NUlist[scen_ind]
        sigma2G= sigma2E/ni+sigma2NU # cluster mean level variance
        # long versions of parameters:
        rhott_long[long_ind] = rhott_list[scen_ind]
        TT_long[long_ind] = TT
        cc_long[long_ind] = cc
        ni_long[long_ind] = ni
        sigma2U_long[long_ind]=sigma2U
        sigma2G_long[long_ind]=sigma2G
        Rkls = fgls5(TT,cc,ni)
        Rkls = Rkls[,,1,1]
        if(small_sample_correction){
          Rkls = Rkls/matrix(rep(c(Rkls[1,1],Rkls[2,2]),4),nrow=2,ncol=4)
        }
        sigma2_hats = Rkls%*%matrix(c(sigma2G,sigma2A,sigma2U,sigmaAU),ncol=1)
        est_var = tevpie(sigma2_hats[1,1],sigma2_hats[2,1],TT,cc)
        est_var_long[long_ind] = est_var
        
        true_covar = diag(CC*TT)*sigma2G+B*sigma2A+Bu*(sigma2U+sigmaAU)+Buc*sigmaAU# equivalent variance of estimator if using cell means,vars vs full expanded model
        true_covar_eta = solve(t(X)%*%solve(true_covar)%*%X)
        true_var = true_covar_eta[1,1]
        true_var_long[long_ind] = true_var
        long_ind = long_ind+1
      }
    }
  }
  var_ratio = est_var_long/true_var_long
  var_ratio_frame = data.frame(var_ratio,real_var=true_var_long,est_var=est_var_long,TT=TT_long,cc=cc_long,ni=ni_long,rhott=rhott_long,sigma2U=sigma2U_long,sigma2G=sigma2G_long)
  return(var_ratio_frame)
}
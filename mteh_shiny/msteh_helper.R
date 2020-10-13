# Helper functions for msteh shiny app
library(gridExtra)
library(ggplot2)
library(reshape2)
library(scales)
library(grid)
library(dplyr)
library(latex2exp)
source('rho_to_sigma_ratio2.R',local=TRUE)
source('fgls5.R',local=TRUE)
source('sw_grid.R')

T_len = 21

# The following two functions are supplied by Jessica Kasza
SWdesmat <- function(T) {
  Xsw <- matrix(data=0, ncol = T, nrow = (T-1))
  for(i in 1:(T-1)) {
    Xsw[i,(i+1):T] <- 1
  }
  return(Xsw)
}
CRXOdesmat<- function(T) {
  if((T-1)%%2 == 0) {
    
    Xcrxo <- matrix(data=0, ncol = T, nrow = T-1)
    Xcrxo[1:(T-1)/2, seq(1,T,2)] <- 1
    Xcrxo[((T-1)/2 + 1):(T-1), seq(2,T,2)] <- 1
    return(Xcrxo)    
  }
  if((T-1)%%2 == 1) {
    
    Xcrxo <- matrix(data=0, ncol = T, nrow = T)
    Xcrxo[1:(T)/2, seq(1,T,2)] <- 1
    Xcrxo[((T)/2 + 1):(T), seq(2,T,2)] <- 1
    return(Xcrxo)    
  }
}

# gives the highest allowable value of rhoct
rhoct_max = function(rhocc,rhott_min,rhott_max){
  max_rhoct = function(rhottm,rhoccm){return((rhoccm/(1-rhoccm)+rhottm/(1-rhottm))*0.5*sqrt((1-rhoccm)*(1-rhottm)))}
  optimize(max_rhoct,rhoccm=rhocc,lower=rhott_min,upper=rhott_max)
}

# return the appropriate graph based on which design is being used.
my_outer_graph = function(design, panel_var, cc, ni, rhoCC, rhoCT, rhoCAC, rhoTTrange){
  # design: "Stepped wedge" = 1, "CRXO" = 2, "SW/CXO ratio"=3
  # panel_var: "c"=1,"N"=2
  if(design==1){
    if(panel_var==1){
      return(msteh_results_grid_swR(ni=ni, cc=c(1,3,5), rhoCC=rhoCC, rhoCT=rhoCT, rhoCAC=rhoCAC, panel_var=panel_var, rhoTTrange=rhoTTrange))
    }
    if(panel_var==2){
      return(msteh_results_grid_swR(ni=c(5,20,100,500), cc=cc, rhoCC=rhoCC, rhoCT=rhoCT, rhoCAC=rhoCAC, panel_var=panel_var, rhoTTrange=rhoTTrange))
    }
  }
  if(design==2){
    # panel_var: "c" = 1, "N" = 2
    if(panel_var==1){
      return(msteh_results_grid_cxoC(ni=ni, rhoCC=rhoCC, rhoCT=rhoCT, rhoCAC=rhoCAC, rhoTTrange=rhoTTrange))
    }
    if(panel_var==2){
      return(msteh_results_grid_cxoN(cc=cc, rhoCC=rhoCC, rhoCT=rhoCT, rhoCAC=rhoCAC, rhoTTrange=rhoTTrange))
    }
  }
  if(design==3){
    #panel_var is irrelevant in this case: c is fixed to 1
    return(msteh_results_grid_ratio(ni=ni, rhoCC=rhoCC, rhoCT=rhoCT, rhoCAC=rhoCAC, rhoTTrange=rhoTTrange))
  }
}

# make a data frame of sw variance ratios
msteh_results_grid_swR = function(Tlist=seq(4,15,1),cc,sigmaE2=1,ni,rhoCC,rhoCT,rhoCAC,panel_var=1,rhoTTrange=seq(0.01,0.02,0.005)){
  #Tlist = seq(4, 20 ,1) 
  #Clist = c(2,5,20)
  #sigmaE2 = 1
  #n=20
  cat("Calculating ratios\n")
  rhott_list = seq(rhoTTrange[1],rhoTTrange[2],length.out=T_len)
  if(panel_var==1){
    ratio_data = sw_grid(T_list=Tlist,rhott_list=rhott_list,panel_var_list=cc,fixed_var=ni,rho_cc=rhoCC,rho_ct=rhoCT,rho_ac=rhoCAC,panel_var=panel_var)
  }else{
    ratio_data = sw_grid(T_list=Tlist,rhott_list=rhott_list,panel_var_list=ni,fixed_var=cc,rho_cc=rhoCC,rho_ct=rhoCT,rho_ac=rhoCAC,panel_var=panel_var)
  }
  cat("Plotting ratios\n")
  myplot = plot_panels(ratio_data,panel_var)
  return(myplot)
}

# make a data frame of cxo variance ratios with multiple values of c
msteh_results_grid_cxoC = function(Tlist=seq(4,20,1),clist=c(2,5,20),sigmaE2=1,ni=20,rhoCC=0.01,rhoCT=0.01,rhoCAC=0.8,rhoTTrange=seq(0.01,0.02,0.005)){
  #Tlist = seq(4, 20 ,1) 
  #Clist = c(2,5,20)
  #sigmaE2 = 1
  #n=20
  rhott_list = seq(rhoTTrange[1],rhoTTrange[2],length.out=T_len)
  ratio_data = expand.grid(TT=Tlist,cc=clist,rhott=rhott_list)
  rhott_long = ratio_data$rhott
  scenarios = t(matrix(rep(c(rhoCC,0.01,rhoCT,rhoCAC),length(rhott_long)),nrow=4))
  scenarios[,2] = rhott_long
  sig2rats = rho_to_sigma_ratio2(scenarios[,1],scenarios[,2],scenarios[,3],scenarios[,4])
  sigmaU2 = sig2rats[[2]]*sigmaE2
  sigmaNU2 = sig2rats[[4]]*sigmaE2
  sigEbar2 = sigmaNU2+sigmaE2/ni
  Tlong = ratio_data$TT
  Clong = 2*ratio_data$cc
  
  est_var = (4/(Clong*Tlong))*sigEbar2+((Clong-2)/(Clong*(Clong*Tlong-Clong-Tlong)))*sigmaU2
  real_var = (4/(Clong*Tlong))*sigEbar2+sigmaU2/Clong
  ratio_data$est_var = est_var
  ratio_data$real_var = real_var
  ratio_data$var_ratio = est_var/real_var
  
  myplot = plot_panels(ratio_data,1)
  return(myplot)
}

# make a data frame of cxo variance ratios with multiple values of N
msteh_results_grid_cxoN = function(Tlist=seq(4,20,1),cc=5,sigmaE2=1,Nlist=c(5,20,100,1000),rhoCC=0.01,rhoCT=0.01,rhoCAC=0.8,rhoTTrange=seq(0.01,0.02,0.005)){
  # Tlist = seq(4, 20 ,1) 
  # Nlist = c(5,20,100,1000)
  # sigmaE2 = 1
  # C=5
  rhott_list = seq(rhoTTrange[1],rhoTTrange[2],length.out=T_len)
  ratio_data = expand.grid(TT=Tlist,N=Nlist,rhott=rhott_list)
  rhott_long = ratio_data$rhott
  N          = ratio_data$N
  scenarios = t(matrix(rep(c(rhoCC,0.01,rhoCT,rhoCAC),length(rhott_long)),nrow=4))
  scenarios[,2] = rhott_long
  sig2rats = rho_to_sigma_ratio2(scenarios[,1],scenarios[,2],scenarios[,3],scenarios[,4])
  sigmaU2  = sig2rats[[2]]*sigmaE2
  sigmaNU2 = sig2rats[[4]]*sigmaE2
  sigEbar2 = sigmaNU2+sigmaE2/N
  Tlong = ratio_data$TT
  Clong = 2*cc
  
  est_var = (4/(Clong*Tlong))*sigEbar2+((Clong-1)/(Clong*(Clong*Tlong-Clong-Tlong)))*sigmaU2
  real_var = (4/(Clong*Tlong))*sigEbar2+sigmaU2/Clong
  ratio_data$est_var = est_var
  ratio_data$real_var = real_var
  ratio_data$var_ratio = est_var/real_var
  
  myplot = plot_panels(ratio_data,2)
  return(myplot)
}

# find the ratio of CXO misspecification variance ratio to SW misspecification variance ratio 
msteh_results_grid_ratio = function(Tlist=seq(4,15,1),sigmaE2=1,ni,rhoCC,rhoCT,rhoCAC,rhoTTrange=seq(0.01,0.02,0.005)){
  rhott_list = seq(rhoTTrange[1],rhoTTrange[2],length.out=T_len)
  ratio_dataSW = sw_grid(T_list=Tlist,rhott_list=rhott_list,panel_var_list=2,fixed_var=ni,rho_cc=rhoCC,rho_ct=rhoCT,rho_ac=rhoCAC,panel_var=1)
  # (var_ratio,real_var=true_var_long,est_var=est_var_long,TT=TT_long,cc=cc_long,ni=ni_long,rhott=rhott_long,sigma2U=sigma2U_long,sigma2G=sigma2G_long)
  var_ratioSW = ratio_dataSW$var_ratio
  rhocc = 0.01
  rhoct = 0.01
  TlongSW = ratio_dataSW$TT
  TlongX = ratio_dataSW$TT
  ClongX = (TlongX-1)*2
  rhottlongX = ratio_dataSW$rhott
  sigma2UlongX = ratio_dataSW$sigma2U # just for dimension, will reallocate below
  sigma2GlongX = ratio_dataSW$sigma2G
  est_varX = (4/(ClongX*TlongX))*sigma2GlongX+((ClongX-2)/(ClongX*(ClongX*TlongX-ClongX-TlongX)))*sigma2UlongX
  real_varX = (4/(ClongX*TlongX))*sigma2GlongX+sigma2UlongX/ClongX
  var_ratioX = est_varX/real_varX
  
  TlongSWX = rep(TlongX,3)
  rhottlongSWX = rep(rhottlongX,3)
  var_ratioSWX = c(var_ratioX,var_ratioSW,var_ratioX/var_ratioSW)
  design = c(rep(1,length(TlongX)),rep(2,length(TlongSW)),rep(3,length(TlongSW)))
  ratio_dataSWX = data.frame(TT=TlongSWX,rhott=rhottlongSWX,var_ratio=var_ratioSWX,Design=design)
  
  ratio_dataSWX$Design <- as.factor(ratio_dataSWX$Design)
  levels(ratio_dataSWX$Design)[levels(ratio_dataSWX$Design)=="1"]   <-"CXO"
  levels(ratio_dataSWX$Design)[levels(ratio_dataSWX$Design)=="2"]   <-"SW"
  levels(ratio_dataSWX$Design)[levels(ratio_dataSWX$Design)=="3"]   <-"CXO/SW"
  return(plot_panels(ratio_dataSWX,3))
}

# turn a long form data frame into a graph
# panel_var is 1 for c panels and 2 for N panels
plot_panels = function(ratio_data,panel_var){
  if(panel_var==1){# c panels
    names(ratio_data)[names(ratio_data)=='cc']='c' # to make the labelling work
    myplot <- ggplot(ratio_data, aes(x=rhott, y=TT)) + 
      facet_grid( c ~ .,labeller=labeller(c=label_both)) + 
      geom_tile(aes(fill = var_ratio)) + 
      scale_fill_gradientn(colours=c("red","white"),
                           values  = c(0,1)) +
      scale_x_continuous(expand=c(0,0)) +
      scale_y_continuous(expand=c(0,0)) +
      theme(aspect.ratio=3/8, legend.position ="none", legend.key.size = unit(1, "cm"), 
            legend.text=element_text(size=12), 
            legend.background = element_rect(fill="grey95")) +
      coord_fixed() + xlab(TeX("$\\rho_{TT}$")) +  ylab("Number of periods, T") +
      guides(fill=guide_legend(nrow=1, keywidth=2, unit="cm")) +
      geom_text(aes(rhott, TT, label = round(var_ratio,2)), color = "black", size = 1.5) 
    return(myplot)
  }
  if(panel_var==2){# N panels
    names(ratio_data)[names(ratio_data)=='ni']='N' # to make the labelling work
    myplot <- ggplot(ratio_data, aes(x=rhott, y=TT)) + 
      facet_grid( N ~ .,labeller=labeller(N=label_both)) + 
      geom_tile(aes(fill = var_ratio)) + 
      scale_fill_gradientn(colours=c("red","white"),
                           values  = c(0,1)) +
      scale_x_continuous(expand=c(0,0)) +
      scale_y_continuous(expand=c(0,0)) +
      theme(aspect.ratio=3/8, legend.position ="none", legend.key.size = unit(1, "cm"), 
            legend.text=element_text(size=12), 
            legend.background = element_rect(fill="grey95")) +
      coord_fixed() + xlab(TeX("$\\rho_{TT}$")) +  ylab("Number of periods, T") +
      guides(fill=guide_legend(nrow=1, keywidth=2, unit="cm")) +
      geom_text(aes(rhott, TT, label = round(var_ratio,2)), color = "black", size = 1)
    return(myplot)
  }
  if(panel_var==3){
    myplot <- ggplot(ratio_data, aes(x=rhott, y=TT)) + 
      facet_grid( Design ~ .) + 
      geom_tile(aes(fill = var_ratio)) + 
      scale_fill_gradientn(colours=c("red","white","blue"),
                           values  = rescale(c(min(ratio_data$var_ratio), 1, max(ratio_data$var_ratio)))) +
      scale_x_continuous(expand=c(0,0)) +
      scale_y_continuous(expand=c(0,0)) +
      theme(aspect.ratio=3/8, legend.position ="none", legend.key.size = unit(1, "cm"),
            legend.text=element_text(size=12),
            legend.background = element_rect(fill="grey95")) +
      coord_fixed() + xlab(TeX("$\\rho_{TT}$")) +  ylab("Number of periods, T") +
      guides(fill=guide_legend(nrow=1, keywidth=2, unit="cm")) +
      geom_text(aes(rhott, TT, label = round(var_ratio,2)), color = "black", size = 1.5) 
    return(myplot)
  }
}
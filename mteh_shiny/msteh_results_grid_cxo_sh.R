# graph results of msteh in a grid
# stepped wedge, T vs sigU vs rho
library(gridExtra)
library(ggplot2)
library(reshape2)
library(scales)
library(grid)
library(dplyr)
library(latex2exp)

specify_ICCs = 1

source('./TEH/rho_to_sigma_ratio2.R',local=TRUE)
ht=10
wd=17
Tlist = seq(4, 20 ,1) 
Clist = c(2,5,20)
sigmaE2 = 1
n=20
if(specify_ICCs){
  rhott_list = seq(0.01,0.05,0.002)
  ratio_data = expand.grid(TT=Tlist,C=Clist,rhott=rhott_list)
  rhott_long = ratio_data$rhott
  scenarios = t(matrix(rep(c(0.01,0.01,0.01,0.8),length(rhott_long)),nrow=4))
  scenarios[,2] = rhott_long
  sig2rats = rho_to_sigma_ratio2(scenarios[,1],scenarios[,2],scenarios[,3],scenarios[,4])
  sigmaU2 = sig2rats[[2]]*sigmaE2
  sigmaNU2 = sig2rats[[4]]*sigmaE2
  sigEbar2 = sigmaNU2+sigmaE2/n
}else{
  sigU2list = seq(0,4,0.2) # var
  ratio_data = expand.grid(TT=Tlist,C=Clist,sigmaU2=sigU2list)
  sigmaU2 = ratio_data$sigmaU2
  sigEbar2 = sigmaE2/n
}
Tlong = ratio_data$TT
Clong = ratio_data$C

est_var = (4/(Clong*Tlong))*sigEbar2+((Clong-2)/(Clong*(Clong*Tlong-Clong-Tlong)))*sigmaU2
real_var = (4/(Clong*Tlong))*sigEbar2+sigmaU2/Clong
ratio_data$est_var = est_var
ratio_data$real_var = real_var
ratio_data$var_ratio = est_var/real_var

myplot <- ggplot(ratio_data, aes(x=rhott, y=TT)) + 
  facet_grid( C ~ .,labeller=labeller(c=label_both)) + 
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

ggsave("rplots/VarRatioCXO_C.pdf", plot=myplot, width= wd, height=ht, units="cm", dpi=600)
#ggsave("VarRatioHH.eps", plot=myplot, width= wd, height=ht, units="cm", dpi=600)


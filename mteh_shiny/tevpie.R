tevpie = function(sehat,sahat,TT,cc){
SS = TT-1
II = SS*cc # total number of clusters
sigma2 = sehat# sigma^2 from Hussey and Hughes
tau2 = sahat # tau^2 from Hussey and Hughes
treat_grid_s = matrix(0,nrow=SS,ncol=TT)
treat_grid_s[upper.tri(treat_grid_s,diag=FALSE)]=1
treat_grid_c = kronecker(treat_grid_s,matrix(1,nrow=cc,ncol=1))
UU = sum(colSums(treat_grid_c))
WW = sum(colSums(treat_grid_c)^2)
VV = sum(rowSums(treat_grid_c)^2)
vari = II*sigma2*(sigma2+TT*tau2)/((II*UU-WW)*sigma2+(UU^2+II*TT*UU-TT*WW-II*VV)*tau2)
return(vari)
}
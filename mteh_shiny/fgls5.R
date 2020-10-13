# FGLS for stepped wedge designs
# returns array of coefficients of sigma terms
# Rkls[1,] is sigma_gamma^2 as a linear combination of sigma_eps^2,sigma_a^2,sigma_u^2,sigma_{au}
# Rkls[2,] is sigma_a^2 as a linear combination of the same.
#function Rkls = fgls5(Tlist,Clist,ni)
library('matlib')
fgls5=function(Tlist,Clist,ni){
    #Tlist = 5#4:20#[4 10 20 30 40 50];%[5 10 20 30 40 50 60 70];
    #Clist = 6#c(1,2,3)#[1 2 3 4];
    #ni=3# ni = 20;
    Rkls = array(data=0,dim=c(2,4,length(Tlist),length(Clist)))#-3*ones(2,5,length(Tlist),length(Clist));
    # scaled versions
    for(T_ind in 1:length(Tlist)){
        for(C_ind in 1:length(Clist)){
            TT = Tlist[T_ind]
            S = TT-1# number of sequences
            cc=Clist[C_ind]
            C=S
            ##fprintf(1,'T: %2d, S: %2d, C:%2d, n:%d\n',T,S,c,ni);
            seq_grid = matrix(0,nrow=S,ncol=TT)
            seq_grid[upper.tri(seq_grid,diag=FALSE)]=1 # sequence period design grid (rows are sequences, columns are time periods, cell is 1 iff treatment exists in that cell)
            cluster_gridb = seq_grid # cluster period design grid (rows are clusters). Base version (c=1) is same as seq_grid
            cluster_vecb = matrix(t(cluster_gridb),ncol=1)
            T_labelsb = t(matrix(1:TT,TT,C))#T_labelsb = repmat(1:T,C,1);
            Tl_vecb = matrix(t(T_labelsb),ncol=1) #Tl_vecb = reshape(T_labelsb',[],1); % Tl_vec(i) is the time period to which the ith data point belongs
            C_labelsb = matrix(1:C,nrow=C,ncol=TT)#C_labelsb = repelem((1:(C))',1,T);
            Cl_vecb = matrix(t(C_labelsb),ncol=1)#Cl_vecb = reshape(C_labelsb',[],1); % Cl_vec(i) is the cluster to which the ith data point belongs
            TEH_vecb = Cl_vecb*cluster_vecb;
            X_betab = 0+(kronecker(matrix(1,nrow=1,ncol=TT),Tl_vecb)==matrix(rep(1:TT,each=C*TT),C*TT,TT))# part of the design matrix corresponding to time period effects, it is C*S*T by T. Rows are data points, columns are time effects. X_beta(i,j)=1 iff ith data point is in jth time period.
            X_thetab = cluster_vecb # treatment indicators
            Xb = cbind(X_thetab,X_betab);
            
            # estimation
            # assume y is data vector (column of length CT)
            olsmb = inv(t(Xb)%*%Xb)%*%t(Xb)# % matrix for getting ordinary least squares estimates of fixed effects, olsm*y
            
            # newer estimation
            Pallb = diag(C*TT);
            P0b = (1/(C*TT))*matrix(1,nrow=C*TT,ncol=C*TT)# projection matrix for mean subspace
            Pcb = (1/(TT))*kronecker(diag(C),matrix(1,TT,TT))# projection matrix for cluster subspace
            Ptb = (1/C)*(X_betab%*%t(X_betab))# projection matrix for time subspace
            Wcb = Pcb-P0b# projection matrix for cluster subspace sans mean
            Wtb = Ptb-P0b# projection matrix for time subspace sans mean
            Wrb = Pallb-Wcb-Wtb-P0b
            gRb = sum((Wrb%*%X_thetab)^2) # bias in R due to TE
            gCb = sum((Wcb%*%X_thetab)^2) # bias in C due to TE
            dR = cc*S*TT*ni-TT-cc*S+1 # degrees of freedom (dim) of residual subspace
            dC = cc*S-1 # degrees of freedom (dim) of cluster subspace, wrt the mean
            
            e1 = matrix(0,nrow=TT+1,ncol=1);
            e1[1] = 1 # elementary unit vector
            Hthetab = t(e1)%*%olsmb# hat matrix for theta
            
            Xtehb = matrix(0,C*TT,C);
            for(k in 1:C){
                Xtehb[,k] = 0+(TEH_vecb ==k);
            }
            # aliases
            Zcb = kronecker(diag(C),matrix(1,TT,1))
            Zeb = diag(C*TT)
            Zub = Xtehb
            
            # matrix K of base components
            # Tr(Z'*P*Z), where 
            # rows correspond to Z=(Ze,Zc,Zu,Zu x Zc)
            # columns correspond to P=(P0,Pc,Pt,Pall,Htheta)
            ZZ1 = list(Zeb,Zcb,Zub,Zcb)
            ZZ2 = list(Zeb,Zcb,Zub,2*Zub)
            PP = list(P0b,Pcb,Ptb,Pallb,t(Hthetab)%*%Hthetab)
            numZ= length(ZZ1)
            numP = length(PP)
            Kbase = matrix(0,numZ,numP)
            for(ZZind in 1:numZ){
                for(PPind in 1:numP){
                    Kbase[ZZind,PPind] = tr(t(ZZ1[[ZZind]])%*%PP[[PPind]]%*%ZZ2[[ZZind]])
                }
            }
            Kscaling_factors = t(matrix(c(1,cc,1,cc*ni,1, ni,cc*ni,ni,cc*ni,ni, ni,cc*ni,ni,cc*ni,ni, ni,cc*ni,ni,cc*ni,ni),5,4))
            KK = Kbase*Kscaling_factors
            Rkls[1,1:4,T_ind,C_ind] = t((1/dR)*KK%*%matrix(c(1,-1,-1,1,-gRb),5,1))
            Rkls[2,1:4,T_ind,C_ind] = (1/(dC*ni*TT))*(matrix(c(-1,1,0,0,-gCb),1,5)%*%t(KK)-dC*Rkls[1,1:4,T_ind,C_ind]);
        }
    }
    return(Rkls)
}

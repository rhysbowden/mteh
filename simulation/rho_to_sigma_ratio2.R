# v2 includes cluster-period term
# Helper function for converting correlation terms to variance components (actually the ratio of each variance component to sigma_epsilon^2, the individual error term)
rho_to_sigma_ratio2 <- function(rhocc,rhott,rhoct,rhoac){
  sigArat = rhocc/(1-rhocc);
  sigAUrat = rhoct/sqrt((1-rhocc)*(1-rhott)) - sigArat;
  sigUrat = rhocc/(1-rhocc)+rhott/(1-rhott)-2*rhoct/sqrt((1-rhocc)*(1-rhott));
  sigNUcorrection_factor = rhoac*(1-rhocc)/(rhoac-rhocc);
  sigArat = sigArat*sigNUcorrection_factor;
  sigUrat = sigUrat*sigNUcorrection_factor;
  sigAUrat = sigAUrat*sigNUcorrection_factor;
  sigNUrat = (1/rhoac-1)*sigArat;
  return(list(sigArat,sigUrat,sigAUrat,sigNUrat))
}
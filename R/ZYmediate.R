ZYmediate <- function(x, y, med, nboot = 2000, alpha = 0.05, kappa = 0.05, ...){
  #
  # Robust mediation analysis using M-estimator as
  # described in Zu and Yuan, 2010, MBR, 45, 1--44.
  #
  # x[,1] is predictor
  # x[,2] is mediator variable (m)
  #  y is outcome variable.
  cl <- match.call()
  SEED=FALSE
  xout=FALSE
  ep=0.00000001  # convergence criteria
  B=nboot         # the number of bootstrap replications
  #kappa    # the percent of cases to be controlled when robust method is used
  # Zu and Yuan used .05, so this is the default value used here.
  level=alpha    # alpha level

  Z=elimna(cbind(x,med,y))

  p=3
  n=nrow(Z)
  HT=HuberTun(kappa,p)
  r=HT$r
  tau=HT$tau
  H=robEst(Z,r,tau,ep)
  R.v=H$u2*tau
  oH=order(R.v)
  oCaseH=(1:n)[oH]        # case number with its Ri increases
  oR.v=R.v[oH]

  thetaH=H$theta
  aH=thetaH[1]
  bH=thetaH[2]
  abH=aH*bH

  muH=H$mu
  SigmaH=H$Sigma
  dH=H$d


  ### Use robust method
  # point estimate
  thetaH=H$theta
  aH=thetaH[1]
  bH=thetaH[2]
  abH=aH*bH

  muH=H$mu
  SigmaH=H$Sigma
  dH=H$d

#   #Standard errors
#   RH=SErob(Z,muH,SigmaH,thetaH,dH,r,tau)
#
#   Zr=RH$Zr
#   SEHI=RH$inf
#   SEHS=RH$sand
#
#   #Standard errors
#   RH=SErob(Z,muH,SigmaH,thetaH,dH,r,tau)
#
#   Zr=RH$Zr
#   SEHI=RH$inf
#   SEHS=RH$sand
#
  #Standard errors
  RH=SErob(Z,muH,SigmaH,thetaH,dH,r,tau)

  Zr=RH$Zr
  SEHI=RH$inf
  SEHS=RH$sand
  ParEstH<-round(cbind(thetaH,SEHI[1:6],SEHS[1:6]),3)
  rnames<-c("a","b","c","vx","vem","vey")
  ParEstH<-cbind(rnames,ParEstH)
  res=t(ParEstH)
  #
  Res=BCI(Z,Zr,ab=3,abH,B,level)
  result <- list(CI.ab=Res$CI,p.value=Res$pv,a.est=aH,b.est=bH,ab.est=abH, alpha = alpha, call = cl)
  class(result) <- "robmed"
  result
}

yuen.effect.ci<-function(formula, data, tr = 0.2, nboot = 400, alpha = 0.05, ...){
  #
  # Compute a 1-alpha  confidence interval
  # for a robust, heteroscedastic  measure of effect size
  #  The absolute value of the measure of effect size is used.
  #
  if (missing(data)) {
    mf <- model.frame(formula)
  } else {
    mf <- model.frame(formula, data)
  }
  cl <- match.call()

  xy <- split(model.extract(mf, "response"), mf[,2])
  faclevels <- names(xy)
  x <- xy[[1]]
  y <- xy[[2]]

  x=elimna(x)
  y=elimna(y)
  bvec=0
  datax<-matrix(sample(x,size=length(x)*nboot,replace=TRUE),nrow=nboot)
  datay<-matrix(sample(y,size=length(x)*nboot,replace=TRUE),nrow=nboot)
  for(i in 1:nboot){
    bvec[i]=yuenv2(datax[i,],datay[i,],tr=tr,SEED=FALSE)$Effect.Size
  }
  bvec<-sort(abs(bvec))
  crit<-alpha/2
  icl<-round(crit*nboot)+1
  icu<-nboot-icl
  ci<-NA
  ci[1]<-bvec[icl]
  pchk=yuen(formula=formula, data = mf,tr=tr)$p.value
  if(pchk>alpha)ci[1]=0
  ci[2]<-bvec[icu]
  if(ci[1]<0)ci[1]=0
  es=abs(yuenv2(x,y,tr=tr)$Effect.Size)
  list(effsize = es, alpha = alpha, CI=ci)
}

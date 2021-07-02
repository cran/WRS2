akp.effect <- function(formula, data, EQVAR = TRUE, tr = 0.2, nboot = 200, alpha = 0.05, ...){
  #
  # Computes the robust effect size suggested by
  #Algina, Keselman, Penfield Psych Methods, 2005, 317-328

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

  x<-elimna(x)
  y<-elimna(y)
  n1<-length(x)
  n2<-length(y)

  ## effect size computation
  s1sq=winvar(x,tr=tr)
  s2sq=winvar(y,tr=tr)
  spsq<-(n1-1)*s1sq+(n2-1)*s2sq
  sp<-sqrt(spsq/(n1+n2-2))
  cterm=1
  if(tr>0)cterm=area(dnormvar,qnorm(tr),qnorm(1-tr))+2*(qnorm(tr)^2)*tr
  cterm=sqrt(cterm)
  if(EQVAR)dval<-cterm*(tmean(x,tr)-tmean(y,tr))/sp
  if(!EQVAR) dval<-cterm*(tmean(x,tr)-tmean(y,tr))/sqrt(s1sq)

  ## bootstrap CI
  be.f=NA
  for(i in 1:nboot){
    X=sample(x,n1,replace=TRUE)
    Y=sample(y,n2,replace=TRUE)
    s1sq=winvar(X,tr=tr)
    s2sq=winvar(Y,tr=tr)
    spsq<-(n1-1)*s1sq+(n2-1)*s2sq
    sp<-sqrt(spsq/(n1+n2-2))
    cterm=1
    if(tr>0)cterm=area(dnormvar,qnorm(tr),qnorm(1-tr))+2*(qnorm(tr)^2)*tr
    cterm=sqrt(cterm)
    if(EQVAR)dval_b<-cterm*(tmean(X,tr)-tmean(Y,tr))/sp
    if(!EQVAR)dval_b<-cterm*(tmean(X,tr)-tmean(Y,tr))/sqrt(s1sq)
    be.f[i] <- dval_b
  }
  L=alpha*nboot/2
  U=nboot-L
  be.f=sort(be.f)
  ci=be.f[L+1]
  ci[2]=be.f[U]

  ## output
  result <- list(AKPeffect = dval, AKPci = ci, call = cl)
  class(result) = "AKP"
  result
}


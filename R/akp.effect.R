akp.effect <- function(formula, data, EQVAR = TRUE, tr = 0.2){
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
  s1sq=winvar(x,tr=tr)
  s2sq=winvar(y,tr=tr)
  spsq<-(n1-1)*s1sq+(n2-1)*s2sq
  sp<-sqrt(spsq/(n1+n2-2))
  cterm=1
  if(tr>0)cterm=area(dnormvar,qnorm(tr),qnorm(1-tr))+2*(qnorm(tr)^2)*tr
  cterm=sqrt(cterm)
  if(EQVAR)dval<-cterm*(tmean(x,tr)-tmean(y,tr))/sp
  if(!EQVAR){
    dval<-cterm*(tmean(x,tr)-tmean(y,tr))/sqrt(s1sq)
    dval[2]=cterm*(tmean(x,tr)-tmean(y,tr))/sqrt(s2sq)
  }
  dval
}
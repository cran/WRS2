dep.effect <-function(x, y, tr = 0.2, nboot = 1000, ...){
  #
  #
  # For two dependent groups,
  # compute confidence intervals for four measures of effect size based on difference scores:
  #
  #  AKP: robust standardized difference similar to  Cohen's d
  #  QS:  Quantile shift based on the median of the distribution of difference scores,
  #  QStr: Quantile shift based on the trimmed mean of the distribution of X-Y
  #  SIGN:  P(X<Y), probability that for a random pair, the first is less than the second.
  #
  #  OPT=TRUE: No effect, difference scores are symmetric about zero.
  #  Under normality, suppose a shift of .2, .5 and .8 standard deviation of the difference score
  #  is considered small, medium and large. The corresponding values for QS and SIGN are printed.
  #
  #
  #  REL.MAG: suppose it is decided that AKP values .1, .3 and .5 are viewed as small, medium and large under normality Then
  #
  #  if OPT=T and
  # REL.MAG=c(.1,.3,.5)  the default,
  # means that the function will compute the corresponding values for the measures of effect size used here.

  REL.MAG = NULL
  SEED <- FALSE
  ecom=c(0.10,  0.54,  0.54,  0.46, 0.30,  0.62, 0.62,  0.38, 0.50,  0.69,  0.69,  0.31)
  REL.EF=matrix(NA,4,3)

  if(!is.null(y))x=x-y
  x=elimna(x)
  n=length(x)
  output=matrix(NA,ncol=7,nrow=4)
  dimnames(output)=list(c('AKP','QS (median)','QStr','SIGN'),c('NULL','Est','S','M','L','ci.low','ci.up'))
  output[1,1:2]=c(0,D.akp.effect(x, tr=tr))
  output[2,1:2]=c(0.5,depQS(x)$Q.effect)
  output[3,1:2]=c(0.5,depQS(x,locfun=mean,tr=tr)$Q.effect)
  output[4,1:2]=c(0.5,mean(x[x!=0]<0))
  if(is.null(REL.MAG)){
    REL.MAG=c(.1,.3,.5)
    REL.EF=matrix(ecom,4,3)
  }

  if(output[1,2]<0)REL.EF[1,]=0-REL.EF[1,]
  if(output[2,2]<0.5)REL.EF[2,]=.5-(REL.EF[2,]-.5)
  if(output[3,2]<0.5)REL.EF[3,]=.5-(REL.EF[3,]-.5)
  if(output[4,2]>0.5)REL.EF[4,]=.5-(REL.EF[4,]-.5)
  output[,3:5]=REL.EF
  output[1,3:5]=REL.MAG
  #}
  output[1,6:7]=D.akp.effect.ci(x,SEED=SEED,tr=tr,nboot=nboot)$ci
  output[2,6:7]=depQSci(x,SEED=SEED,nboot=nboot)$ci
  output[3,6:7]=depQSci(x,locfun=tmean,SEED=SEED,tr=tr,nboot=nboot)$ci
  Z=sum(x<0)
  output[4,6:7]=binom.conf(Z,n)$ci

  class(output)=c("matrix", "array", "dep.effect")

  output
}


bi2KMSv2<-function(r1=sum(x),n1=length(x),r2=sum(y),n2=length(y),
                   x=NA,y=NA,nullval=0){
  #
  # Test the hypothesis that two independent binomials have equal
  # probability of success using method KMS.
  #
  #  Unlike the function bi2KMS, a p-value is returned
  #
  # r1=number of successes in group 1
  # n1=number of observations in group 1
  #
  # Uses Kulinskaya et al. method American Statistician, 2010, 64, 350-
  #
  #  null value is the hypothesized value for p1-p2
  #
  alph<-c(1:99)/100
  for(i in 1:99){
    irem<-i
    chkit<-bi2KMS(r1=r1,n1=n1,r2=r2,n2=n2,x=x,y=x,alpha=alph[i])
    if(chkit$ci[1]>nullval || chkit$ci[2]<nullval)break
  }
  p.value<-irem/100
  if(p.value<=.1){
    iup<-(irem+1)/100
    alph<-seq(.001,iup,.001)
    for(i in 1:length(alph)){
      p.value<-alph[i]
      chkit<-bi2KMS(r1=r1,n1=n1,r2=r2,n2=n2,x=x,y=x,alpha=alph[i])
      if(chkit$ci[1]>nullval || chkit$ci[2]<nullval)break
    }}
  est=bi2KMS(r1=r1,n1=n1,r2=r2,n2=n2,x=x,y=y)
  list(ci=est$ci,est.p1=est$p1,est.p2=est$p2,p.value=p.value)
}

bi2KMS<-function(r1=sum(x),n1=length(x),r2=sum(y),n2=length(y),
                 x=NA,y=NA,alpha=.05){
  #
  # Test the hypothesis that two independent binomials have equal
  # probability of success
  #
  # r1=number of successes in group 1
  # n1=number of observations in group 1
  #
  # Use Kulinskaya et al. method American Statistician, 2010, 64, 350-
  #
  N=n1+n2
  u=.5
  Dhat=(r1+.5)/(n1+1)-(r2+.5)/(n2+1)
  psihat=((r1+.5)/(n1+1)+(r2+.5)/(n2+1))/2
  nuhat=(1-2*psihat)*(.5-n2/N)
  what=sqrt(2*u*psihat*(1-psihat)+nuhat^2)
  se=qnorm(1-alpha/2)*sqrt(u/(2*n1*n2/N))
  val1=max(c(-1,(u*Dhat+nuhat)/what-se))
  ci=what*sin(asin(val1))/u-nuhat/u
  val2=min(c(1,(u*Dhat+nuhat)/what+se))
  ci[2]=what*sin(asin(val2))/u-nuhat/u
  list(ci=ci,p1=r1/n1,p2=r2/n2)
}


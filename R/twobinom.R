twobinom<-function(r1=sum(elimna(x)),n1=length(x),r2=sum(elimna(y)),n2=length(y),x=NA,y=NA,alpha=.05){
  #
  # Test the hypothesis that two independent binomials have equal
  # probability of success using the Storer--Kim method.
  #
  # r1=number of successes in group 1
  # n1=number of observations in group 1
  #
  n1p<-n1+1
  n2p<-n2+1
  n1m<-n1-1
  n2m<-n2-1
  chk<-abs(r1/n1-r2/n2)
  x<-c(0:n1)/n1
  y<-c(0:n2)/n2
  phat<-(r1+r2)/(n1+n2)
  m1<-outer(x,y,"-")
  m2<-matrix(1,n1p,n2p)
  flag<-(abs(m1)>=chk)
  m3<-m2*flag
  b1<-1
  b2<-1
  xv<-c(1:n1)
  yv<-c(1:n2)
  xv1<-n1-xv+1
  yv1<-n2-yv+1
  dis1<-c(1,pbeta(phat,xv,xv1))
  dis2<-c(1,pbeta(phat,yv,yv1))
  pd1<-NA
  pd2<-NA
  for(i in 1:n1)pd1[i]<-dis1[i]-dis1[i+1]
  for(i in 1:n2)pd2[i]<-dis2[i]-dis2[i+1]
  pd1[n1p]<-phat^n1
  pd2[n2p]<-phat^n2
  m4<-outer(pd1,pd2,"*")
  test<-sum(m3*m4)
  list(p.value=test,p1=r1/n1,p2=r2/n2,est.dif=r1/n1-r2/n2)
}

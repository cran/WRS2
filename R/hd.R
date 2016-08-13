hd<-function(x,q=.5,na.rm=TRUE,STAND=NULL){
  #
  #  Compute the Harrell-Davis estimate of the qth quantile
  #
  #  The vector x contains the data,
  #  and the desired quantile is q
  #  The default value for q is .5.
  #
  if(na.rm)x=elimna(x)
  n<-length(x)
  m1<-(n+1)*q
  m2<-(n+1)*(1-q)
  vec<-seq(along=x)
  w<-pbeta(vec/n,m1,m2)-pbeta((vec-1)/n,m1,m2)  # W sub i values
  y<-sort(x)
  hd<-sum(w*y)
  hd
}

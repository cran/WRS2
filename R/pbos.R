pbos<-function(x,beta=.2){
  #
  #    Compute the one-step percentage bend measure of location
  #
  #
  temp<-sort(abs(x-median(x)))
  omhatx<-temp[floor((1-beta)*length(x))]
  psi<-(x-median(x))/omhatx
  i1<-length(psi[psi<(-1)])
  i2<-length(psi[psi>1])
  sx<-ifelse(psi<(-1),0,x)
  sx<-ifelse(psi>1,0,sx)
  pbos<-(sum(sx)+omhatx*(i2-i1))/(length(x)-i1-i2)
  pbos
}
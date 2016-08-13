pbcor <- function(x, y = NULL, beta = 0.2){
  #   Compute the percentage bend correlation between x and y.
  #
  #   beta is the bending constant for omega sub N.
  #
  if(is.null(y[1])){
    y=x[,2]
    x=x[,1]
  }
  if(length(x)!=length(y))stop("The vectors do not have equal lengths!")
  cl <- match.call()
  
  m1<-cbind(x,y)
  m1<-na.omit(m1)
  nval<-nrow(m1)
  x<-m1[,1]
  y<-m1[,2]
  #  Have eliminated missing values
  temp<-sort(abs(x-median(x)))
  omhatx<-temp[floor((1-beta)*length(x))]
  temp<-sort(abs(y-median(y)))
  omhaty<-temp[floor((1-beta)*length(y))]
  a<-(x-pbos(x,beta))/omhatx
  b<-(y-pbos(y,beta))/omhaty
  a<-ifelse(a<=-1,-1,a)
  a<-ifelse(a>=1,1,a)
  b<-ifelse(b<=-1,-1,b)
  b<-ifelse(b>=1,1,b)
  pbcor<-sum(a*b)/sqrt(sum(a^2)*sum(b^2))
  test<-pbcor*sqrt((length(x) - 2)/(1 - pbcor^2))
  sig<-2*(1 - pt(abs(test),length(x)-2))
  result <- list(cor=pbcor,test=test,p.value=sig,n=nval, call = cl)
  class(result) <- "pbcor"
  result
}

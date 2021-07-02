pbcor <- function(x, y = NULL, beta = 0.2, ci = FALSE, nboot = 500, alpha = 0.05, ...){
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

  ## confidence interval
  if (ci) {
  data<-matrix(sample(length(y),size=length(y)*nboot,replace=TRUE),nrow=nboot)
  bvec <- apply(data, 1, function(i) {
    x1 <- x[i]
    y1 <- y[i]
    temp<-sort(abs(x1-median(x1)))
    omhatx<-temp[floor((1-beta)*length(x1))]
    temp<-sort(abs(y1-median(y1)))
    omhaty<-temp[floor((1-beta)*length(y1))]
    a<-(x1-pbos(x1,beta))/omhatx
    b<-(y1-pbos(y1,beta))/omhaty
    a<-ifelse(a<=-1,-1,a)
    a<-ifelse(a>=1,1,a)
    b<-ifelse(b<=-1,-1,b)
    b<-ifelse(b>=1,1,b)
    pbcor<-sum(a*b)/sqrt(sum(a^2)*sum(b^2))
    pbcor
  })
  ihi<-floor((1-alpha/2)*nboot+.5)
  ilow<-floor((alpha/2)*nboot+.5)
  bsort<-sort(bvec)
  corci<-1
  corci[1]<-bsort[ilow]
  corci[2]<-bsort[ihi]
  } else {
    corci <- NA
  }
  ## end CI

  result <- list(cor=pbcor,test=test,p.value=sig,n=nval, cor_ci = corci, call = cl)
  class(result) <- "pbcor"
  result
}

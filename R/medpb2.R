medpb2 <- function(formula, data, nboot = 2000, ...){
  #
  #   Compare 2 independent groups using medians.
  #
  #   A percentile bootstrap method is used, which performs well when
  #   there are tied values.
  #
  #   The data are assumed to be stored in x and y
  #
  #   Missing values are automatically removed.
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
  alpha = 0.05

  x=elimna(x)
  y=elimna(y)
  xx<-list()
  xx[[1]]<-x
  xx[[2]]<-y

  est1=median(xx[[1]])
  est2=median(xx[[2]])
  est.dif<-median(xx[[1]])-median(xx[[2]])
  crit<-alpha/2
  temp<-round(crit*nboot)
  icl<-temp+1
  icu<-nboot-temp
  bvec<-matrix(NA,nrow=2,ncol=nboot)

  for(j in 1:2){
    data<-matrix(sample(xx[[j]],size=length(xx[[j]])*nboot,replace=TRUE),nrow=nboot)
    bvec[j,]<-apply(data,1,median) # Bootstrapped medians for jth group
  }
  top<-bvec[1,]-bvec[2,]
  test<-sum(top<0)/nboot+.5*sum(top==0)/nboot
  if(test > .5)test<-1-test
  top<-sort(top)
  ci<-NA
  ci[1]<-top[icl]
  ci[2]<-top[icu]

  result <- list(test = est.dif, conf.int = ci, p.value = 2*test, call = cl)
  class(result) <- "pb2"
  result

}

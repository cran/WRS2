twopcor <- function(x1, y1, x2, y2, nboot = 599){
  #
  #   Compute a .95 confidence interval for
  #   the difference between two Pearson
  #   correlations corresponding to two independent
  #   goups.
  #
  #   This function uses an adjusted percentile bootstrap method that
  #   gives good results when the error term is heteroscedastic.
  #
  #   WARNING: If the number of boostrap samples is altered, it is
  #   unknown how to adjust the confidence interval when n1+n2 < 250.
  #
  cl <- match.call()
  if (nboot < 599) warning("It is unknown how to adjust the confidence interval when n1+n2 < 250.")  
  X<-elimna(cbind(x1,y1))
  x1<-X[,1]
  y1<-X[,2]
  X<-elimna(cbind(x2,y2))
  x2<-X[,1]
  y2<-X[,2]
  
  data1<-matrix(sample(length(y1),size=length(y1)*nboot,replace=TRUE),nrow=nboot)
  bvec1 <- apply(data1, 1, function(xx) cor(x1[xx], y1[xx]))
  
  data2<-matrix(sample(length(y2),size=length(y2)*nboot,replace=TRUE),nrow=nboot)
  bvec2<-apply(data2,1,function(xx) cor(x2[xx], y2[xx])) 
  
  bvec<-bvec1-bvec2
  ilow<-15
  ihi<-584
  if(length(y1)+length(y2) < 250){
    ilow<-14
    ihi<-585
  }
  if(length(y1)+length(y2) < 180){
    ilow<-11
    ihi<-588
  }
  if(length(y1)+length(y2) < 80){
    ilow<-8
    ihi<-592
  }
  if(length(y1)+length(y2) < 40){
    ilow<-7
    ihi<-593
  }
  bsort<-sort(bvec)
  r1<-cor(x1,y1)
  r2<-cor(x2,y2)
  ci<-c(bsort[ilow],bsort[ihi])
  result <- list(r1 = r1, r2 = r2, ci = ci, call = cl)
  class(result) <- "twocor"
  result
}

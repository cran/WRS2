wincor <- function(x, y = NULL, tr = 0.2, ci = FALSE, nboot = 500, alpha = 0.05, ...){
#   Compute the Winsorized correlation between x and y.
#
#   tr is the amount of Winsorization
#   This function also returns the Winsorized covariance
#
#    Pairwise deletion of missing values is performed.
#
  cl <- match.call()
  if(is.null(y[1])){
    y=x[,2]
    x=x[,1]
  }
  sig<-NA
  if(length(x)!=length(y))stop("Lengths of vectors are not equal")
  m1=cbind(x,y)
  m1<-elimna(m1)
  nval=nrow(m1)
  x<-m1[,1]
  y<-m1[,2]
  g<-floor(tr*length(x))
  xvec<-winval(x,tr)
  yvec<-winval(y,tr)
  wcor<-cor(xvec,yvec)
  wcov<-var(xvec,yvec)
  if(sum(x==y)!=length(x)){
    test<-wcor*sqrt((length(x)-2)/(1.-wcor^2))
    sig<-2*(1-pt(abs(test),length(x)-2*g-2))
  } else {
    test <- NA
    sig <- NA
  }

## confidence interval
  if (ci) {
    data<-matrix(sample(length(y),size=length(y)*nboot,replace=TRUE),nrow=nboot)
    bvec <- apply(data, 1, function(i) {
      x1 <- x[i]
      y1 <- y[i]
      g<-floor(tr*length(x1))
      xvec<-winval(x1,tr)
      yvec<-winval(y1,tr)
      wcor<-cor(xvec,yvec)
      wcor
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

result <- list(cor = wcor, cov = wcov, test = test, p.value = sig, n = nval, cor_ci = corci, call = cl)
class(result) <- "pbcor"
result
}

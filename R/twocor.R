twocor <- function(x1, y1, x2, y2, corfun = "pbcor", nboot = 599, tr = 0.2, beta = 0.2, ...){
  #
  #  Compute a .95 confidence interval for the
  #  difference between two correlation coefficients
  #  corresponding to two independent groups.
  #
  #   the function corfun is any R function that returns a
  #   correlation coefficient in corfun$cor. The functions pbcor and
  #   wincor follow this convention.
  #
  #   For Pearson's correlation, use
  #   the function twopcor instead.
  #
  #   The default number of bootstrap samples is nboot=599
  #
  cl <- match.call()
  alpha <- .05
  corfun <- match.arg(corfun, c("pbcor", "wincor"), several.ok = FALSE)

  data1<-matrix(sample(length(y1),size=length(y1)*nboot,replace=TRUE),nrow=nboot)
  if (corfun == "pbcor")  bvec1 <- apply(data1, 1, function(xx) pbcor(x1[xx], y1[xx], beta = beta, ci = FALSE)$cor)
  if (corfun == "wincor")  bvec1 <- apply(data1, 1, function(xx) wincor(x1[xx], y1[xx], tr = tr, ci = FALSE)$cor)

  data2<-matrix(sample(length(y2),size=length(y2)*nboot,replace=TRUE),nrow=nboot)
  if (corfun == "pbcor")  bvec2 <- apply(data2, 1, function(xx) pbcor(x2[xx], y2[xx], beta = beta,ci = FALSE)$cor)
  if (corfun == "wincor")  bvec2 <- apply(data2, 1, function(xx) wincor(x2[xx], y2[xx], tr = tr,ci = FALSE)$cor)

  bvec<-bvec1-bvec2
  bsort<-sort(bvec)

  nboot <- length(bsort)

  term<-alpha/2
  ilow<-round((alpha/2) * nboot)
  ihi<-nboot - ilow
  ilow<-ilow+1
  corci<-1
  corci[1]<-bsort[ilow]
  corci[2]<-bsort[ihi]
  #pv<-(sum(bvec<0)+.5*sum(bvec==0))/nboot
  pv<-(sum(bsort<0)+.5*sum(bsort==0))/nboot
  pv=2*min(c(pv,1-pv))

  if (corfun == "pbcor") {
    r1<-pbcor(x1,y1, beta,ci = FALSE)$cor
    r2<-pbcor(x2,y2, beta,ci = FALSE)$cor
  }
  if (corfun == "wincor") {
    r1<-wincor(x1,y1, tr,ci = FALSE)$cor
    r2<-wincor(x2,y2, tr,ci = FALSE)$cor
  }

  result <- list(r1=r1, r2 = r2, ci = corci, p.value = pv, call = cl)
  class(result) <- "twocor"
  result
}

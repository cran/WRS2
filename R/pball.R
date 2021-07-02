pball <- function(x, beta = 0.2, ...){
  #
  #    Compute the percentage bend correlation matrix for the
  #    data in the n by p matrix m.
  #
  #    This function also returns the two-sided significance level
  #    for all pairs of variables, plus a test of zero correlations
  #    among all pairs. (See chapter 6 for details.)
  #

  cl <- match.call()
  m <- x
  pbcorm<-matrix(0,ncol(m),ncol(m))
  temp<-matrix(1,ncol(m),ncol(m))
  siglevel<-matrix(NA,ncol(m),ncol(m))
  cmat<-matrix(0,ncol(m),ncol(m))
  for (i in 1:ncol(m)){
    ip1<-i
    for (j in ip1:ncol(m)){
      if(i<j){
        pbc<-pbcor(m[,i],m[,j],beta)
        pbcorm[i,j]<-pbc$cor
        temp[i,j]<-pbcorm[i,j]
        temp[j,i]<-pbcorm[i,j]
        siglevel[i,j]<-pbc$p.value
        siglevel[j,i]<-siglevel[i,j]
      }
    }
  }
  tstat<-pbcorm*sqrt((nrow(m)-2)/(1-pbcorm^2))
  cmat<-sqrt((nrow(m)-2.5)*log(1+tstat^2/(nrow(m)-2)))
  bv<-48*(nrow(m)-2.5)^2
  cmat<-cmat+(cmat^3+3*cmat)/bv-(4*cmat^7+33*cmat^5+240^cmat^3+855*cmat)/(10*bv^2+8*bv*cmat^4+1000*bv)
  H<-sum(cmat^2)
  df<-ncol(m)*(ncol(m)-1)/2
  h.siglevel<-1-pchisq(H,df)

  if (!is.null(colnames(x))) rownames(temp) <- colnames(temp) <- rownames(siglevel) <- colnames(siglevel) <- colnames(x)

  result <- list(pbcorm = temp, p.values = siglevel, H = H, H.p.value = h.siglevel, call = cl)
  class(result) <- "pball"
  result
}

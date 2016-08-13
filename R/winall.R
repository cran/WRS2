winall<-function(x, tr = 0.2){
#
#    Compute the Winsorized correlation and covariance matrix for the
#    data in the n by p matrix m.
#
#    This function also returns the two-sided significance level
#
m <- x  
cl <- match.call()
if(is.data.frame(m))m=as.matrix(m)
if(!is.matrix(m))stop("The data must be stored in a n by p matrix")
wcor<-matrix(1,ncol(m),ncol(m))
wcov<-matrix(0,ncol(m),ncol(m))
siglevel<-matrix(NA,ncol(m),ncol(m))
for (i in 1:ncol(m)){
  ip<-i
  for (j in ip:ncol(m)){
    val<-wincor(m[,i],m[,j],tr)
    wcor[i,j]<-val$cor
    wcor[j,i]<-wcor[i,j]
    if(i==j)wcor[i,j]<-1
    wcov[i,j]<-val$cov
    wcov[j,i]<-wcov[i,j]
    if(i!=j){
      siglevel[i,j]<-val$p.value
      siglevel[j,i]<-siglevel[i,j]
    }
  }
}

if (!is.null(colnames(x))) rownames(wcor) <- colnames(wcor) <- rownames(wcov) <- colnames(wcov) <- rownames(siglevel) <- colnames(siglevel) <- colnames(x)

result <- list(cor=wcor,cov=wcov,p.values=siglevel, call = cl)
class(result) <- "pball"
result
}


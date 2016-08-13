wincor <- function(x, y = NULL, tr = 0.2){
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
result <- list(cor = wcor, cov = wcov, test = test, p.value = sig, n = nval, call = cl)
class(result) <- "pbcor"
result
}

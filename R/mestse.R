mestse<-function(x,bend=1.28,...){
#
#   Estimate the standard error of M-estimator using Huber's Psi
#   using estimate of influence function
#
n<-length(x)
mestse<-sqrt(sum((ifmest(x,bend,op=2)^2))/(n*(n-1)))
mestse
}

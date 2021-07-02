
onesampb <- function(x, est = "onestep", nboot = 2000, nv = 0, alpha = 0.05, ...){
  #
  #   Compute a bootstrap, .95 confidence interval for the
  #   measure of location corresponding to the argument est.
  #   By default, a one-step
  #   M-estimator of location based on Huber's Psi is used.
  #   The default number of bootstrap samples is nboot=500
  #
  #    nv=null value when computing a p-value
  #
  cl <- match.call()
  est <- match.arg(est, c("mom", "onestep", "median"), several.ok = FALSE)
  est <- get(est)


  null.value <- NULL
  if(!is.null(null.value))nv=null.value
  #if(SEED)set.seed(2) # set seed of random number generator so that
  #             results can be duplicated.
  #print("Taking bootstrap samples. Please wait.")
  x=elimna(x)
  data<-matrix(sample(x,size=length(x)*nboot,replace=TRUE),nrow=nboot)
  bvec<-apply(data,1,est)
  bvec<-sort(bvec)
  low<-round((alpha/2)*nboot)
  up<-nboot-low
  low<-low+1
  pv=mean(bvec>nv)+.5*mean(bvec==nv)
  pv=2*min(c(pv,1-pv))
  estimate=est(x)
  result <- list(ci=c(bvec[low],bvec[up]), n=length(x), estimate=estimate, p.value=pv, alpha = alpha,call=cl)
  class(result) <- "onesampb"
  return(result)
}

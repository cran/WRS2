qcomhd <- function(formula, data, q = c(0.1, 0.25, 0.5, 0.75, 0.9), nboot = 2000, alpha = 0.05, ADJ.CI = TRUE){
  #
  # Compare quantiles using pb2gen
  # via hd estimator. Tied values are allowed.
  # When comparing lower or upper quartiles, both power and the probability of Type I error
  # compare well to other methods that have been derived.
  # q: can be used to specify the quantiles to be compared
  # q defaults to comparing the .1,.25,.5,.75, and .9 quantiles
  #   Function returns p-values and critical p-values based on Hochberg's method.
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
  
  pv=NULL
  output=matrix(NA,nrow=length(q),ncol=10)
  dimnames(output)<-list(NULL,c("q","n1","n2","est.1","est.2","est.1_minus_est.2","ci.low","ci.up", "p_crit","p-value"))
  for(i in 1:length(q)){
    output[i,1]=q[i]
    output[i,2]=length(elimna(x))
    output[i,3]=length(elimna(y))
    output[i,4]=hd(x,q=q[i])
    output[i,5]=hd(y,q=q[i])
    output[i,6]=output[i,4]-output[i,5]
    temp=pb2gen1(x,y,nboot=nboot,est=hd,q=q[i],SEED=FALSE,alpha=alpha,pr=FALSE)
    output[i,7]=temp$ci[1]
    output[i,8]=temp$ci[2]
    output[i,10]=temp$p.value
  }
  temp=order(output[,10],decreasing=TRUE)
  zvec=alpha/c(1:length(q))
  output[temp,9]=zvec
  if(ADJ.CI){
    for(i in 1:length(q)){
      temp=pb2gen1(x,y,nboot=nboot,est=hd,q=q[i],SEED=FALSE,alpha=output[i,9],pr=FALSE)
      output[i,7]=temp$ci[1]
      output[i,8]=temp$ci[2]
      output[i,10]=temp$p.value
    }
  }
  output <- data.frame(output)
  output$signif=rep("YES",nrow(output))
  temp=order(output[,10],decreasing=TRUE)
  for(i in 1:nrow(output)){
    if(output[temp[i],10]>output[temp[i],9])output$signif[temp[i]]="NO"
    if(output[temp[i],10]<=output[temp[i],9])break
  }
 
  output <- output[,-ncol(output)]
  colnames(output) <- c("q", "n1", "n2", "est1", "est2", "est1-est.2", "ci.low", "ci.up", "p.crit", "p.value")
  result <- list(partable = output, call = cl)
  class(result) <- "robtab"
  result
}

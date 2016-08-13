Dqcomhd<-function(x, y, q = c(1:9)/10, nboot = 1000, na.rm = TRUE){
  #
  # Compare the quantiles of the marginal distributions associated with  two dependent groups
  # via hd estimator. Tied values are allowed.
  # When comparing lower or upper quartiles, both power and the probability of Type I error
  # compare well to other methods have been derived.
  #
  #  x: data for group 1
  #  y: data for group 2
  #  q: the quantiles to be compared
  #  nboot: Number of bootstrap samples
  #
  #
  cl <- match.call()
  alpha <- 0.05
  if (na.rm) {
   xy=elimna(cbind(x,y))
   x=xy[,1]
   y=xy[,2]
  }
  pv=NULL
  output=matrix(NA,nrow=length(q),ncol=10)
  dimnames(output)<-list(NULL,c("q","n1","n2","est.1","est.2","est.1_minus_est.2","ci.low","ci.up","p_crit","p-value"))
  
  for(i in 1:length(q)){
    output[i,1]=q[i]
    output[i,2]=length(elimna(x))
    output[i,3]=length(elimna(y))
    output[i,4]=hd(x,q=q[i])
    output[i,5]=hd(y,q=q[i])
    output[i,6]=output[i,4]-output[i,5]
    if(na.rm){
      temp=bootdpci(x,y,est=hd,q=q[i],dif=FALSE,plotit=FALSE,pr=FALSE,nboot=nboot,alpha=alpha,SEED=FALSE)
      output[i,7]=temp$output[1,5]
      output[i,8]=temp$output[1,6]
      output[i,10]=temp$output[1,3]
    }
    if(!na.rm){
      temp=rmmismcp(x,y,est=hd,q=q[i],plotit=FALSE,pr=FALSE,nboot=nboot,alpha=alpha,SEED=FALSE)
      output[i,7]=temp$output[1,6]
      output[i,8]=temp$output[1,7]
      output[i,10]=temp$output[1,4]
    }
  }
  
  temp=order(output[,10],decreasing=TRUE)
  zvec=alpha/c(1:length(q))
  output[temp,9]=zvec
  output <- data.frame(output)
  output$signif=rep("YES",nrow(output))
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

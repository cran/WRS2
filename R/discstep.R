discstep<-function(formula, data, nboot = 500, alpha = 0.05){
  #
  #  Step-down multiple comparison procedure for comparing 
  #  J independent discrete random variables. 
  #  The method is based on a generalization of the Storer--Kim method
  #  comparing independent binomials; it can be sensitive to differences
  #  not detected by measures of location.
  #  
  #  x is a matrix with n rows and J columns
  #  or it can have list mode
  #
  
  if (missing(data)) {
    mf <- model.frame(formula)
  } else {
    mf <- model.frame(formula, data)
  }
  cl <- match.call()
  
  x <- split(model.extract(mf, "response"), mf[,2])   

  
  vals=lapply(x,unique)
  vals=sort(elimna(list2vec(vals)))
  K=length(unique(vals))
  n=lapply(x,length)
  n=list2vec(n)
  J=length(x)
  if(J==2)stop('For 2 groups use disc2com')
  if(J>5)stop('Designed for 5 groups or less')
  com <- com1 <- modgen(J)
  
  ## insert group labels
  lnames <- levels(mf[,2])
  com <- lapply(com, function(xx) lnames[xx])
  
  ntest=length(com)
  jp1=J+1
  com=com[jp1:length(com)]; com1 <- com1[jp1:length(com1)]
  ntest=length(com)
  
  mout= as.data.frame(matrix(NA,nrow=ntest,ncol=3))
  colnames(mout)=c('Groups','p-value','p.crit')
 
  test=NULL
  for(i in 1:ntest){
    test[i]=discANOVA.sub(x[com[[i]]])$test #$
    nmod=length(com[[i]])-1
    temp=c(nmod:0)
    #mout[i,1]=sum(com[[i]]*10^temp)
    mout[i,1] <- paste(com[[i]], collapse = "-")
  }
  mout[,3]=alpha
  xx=list()
  pv=NA
  jm2=J-2
  mout[,3]=alpha
  TB=matrix(NA,nrow=nboot,ncol=ntest)
  step1=discANOVA.sub(x)
  C1=step1$C1
  HT=NULL
  
  for(i in 1:K)HT[i]=mean(C1[i,])
  for(ib in 1:nboot){
    xx=list()
    for(j in 1:J){
      temp=rmultinomial(n[j],1,HT)
      xx[[j]]=which(temp[1,]==1)
      for(i in 2:n[j])xx[[j]][i]=which(temp[i,]==1)
    }
    for(k in 1:ntest)TB[ib,k]=discANOVA.sub(xx[com1[[k]]])$test #$
  }
  for(k in 1:ntest){
    mout[k,2]=1-mean(test[k]>TB[,k])-.5*mean(test[k]==TB[,k])
    pnum=length(com[[k]])
    pe=1-(1-alpha)^(pnum/J)
    if(length(com[[k]])<=jm2)mout[k,3]=pe
  }
  
  outtable <- mout[nrow(mout):1,]
  rownames(outtable) <- NULL
  result <- list(partable = outtable, call = cl)
  class(result) <- "robtab"
  result
}

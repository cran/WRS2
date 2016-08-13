discANOVA<-function(formula, data, nboot = 500){
  #
  #  Test the global hypothesis that for two or more independent groups,
  #  the corresponding discrete distributions are identical. 
  #  That is, test the hypothesis that independent groups have identical
  #  multinomial distributions. A generalization of the Storer--Kim method is used.
  #
  #  Could also use a chi-squared test via the function: disc2.chi.sq
  #
  #  The method is designed for situations where the cardinality of the
  #  sample space is relatively small. The method can be sensitive to
  #  differences that are missed using a measure of location.
  #
  #  Control over the Type I error probability is excellent, even when n=10
  #
  #  x is a matrix with n rows and J columns
  #  or it can have list mode.
  #
  
  if (missing(data)) {
    mf <- model.frame(formula)
  } else {
    mf <- model.frame(formula, data)
  }
  cl <- match.call()
  
  x <- split(model.extract(mf, "response"), mf[,2])   
  
  #if(is.matrix(x) || is.data.frame(x))x=listm(x)
  vals=lapply(x,unique)
  vals=sort(elimna(list2vec(vals)))
  K=length(unique(vals))
  n=lapply(x,length)
  n=list2vec(n)
  J=length(x)
  step1=discANOVA.sub(x)
  test=step1$test
  C1=step1$C1
  HT=NULL
  for(i in 1:K)HT[i]=mean(C1[i,])
  tv=NULL
  TB=NA
  VP=NA
  B1hat=NA
  xx=list()
  for(ib in 1:nboot){
    xx=list()
    for(j in 1:J){
      temp=rmultinomial(n[j],1,HT)
      xx[[j]]=which(temp[1,]==1)
      for(i in 2:n[j])xx[[j]][i]=which(temp[i,]==1)
    }
    TB[ib]=discANOVA.sub(xx)$test
  }
  pv=1-mean(test>TB)-.5*mean(test==TB)
  
  result <- list(test = test, crit.val = NA, p.value = pv, call = cl)
  class(result) <- "med1way"
  result
}


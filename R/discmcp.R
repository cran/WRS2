discmcp <- function(formula, data, alpha = 0.05, nboot = 500, ...){
  #
  #   Multiple comparisons for  J independent groups
  #   having discrete distributions.
  #   The method is based on a chi-squared test for each pair of groups to be compared
  #
  #   The data are assumed to be stored in x
  #   which either has list mode or is a matrix.  In the first case
  #   x[[1]] contains the data for the first group, x[[2]] the data
  #   for the second group, etc. Length(x)=the number of groups = J.
  #   If stored in a matrix, the columns of the matrix correspond
  #   to groups.
  #
  #   Missing values are allowed.
  #
  # Probability of one or more Type I errors controlled using Hochberg's method.
  #

  if (missing(data)) {
    mf <- model.frame(formula)
  } else {
    mf <- model.frame(formula, data)
  }
  cl <- match.call()

  x <- split(model.extract(mf, "response"), mf[,2])

  J<-length(x)
  ncon=(J^2-J)/2
  Jm<-J-1
  #
  # Determine critical values
  dvec=alpha/c(1:ncon)

  output<-as.data.frame(matrix(NA,nrow=ncon,ncol=4))
  colnames(output)<-c('Group 1','Group 2','p.value','p.crit')
  ic=0
  for(j in 1:J){
    for(k in 1:J){
      if(j<k){
        ic=ic+1
        #output[ic,1]=j
        output[ic,1] <- levels(mf[,2])[j]
        #output[ic,2]=k
        output[ic,2] <- levels(mf[,2])[k]
        output[ic,3]=disc2com(x[[j]],x[[k]],simulate.p.value = TRUE, B=nboot)$p.value
      }}}
  temp2<-order(0-output[,3])
  zvec<-dvec[1:ncon]
  #sigvec<-(test[temp2]>=zvec)
  output[temp2,4]<-zvec
  num.sig<-sum(output[,3]<=output[,4])

  outtable <- output
  result <- list(partable = outtable, call = cl)
  class(result) <- "robtab"
  result

  ## FIXME output
}

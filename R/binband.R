binband <- function(x, y, KMS = FALSE, alpha = 0.05, ADJ.P = FALSE){
  #
  #  Comparing two independent variables in terms of their probability function.
  #  For each value that occurs, say x, test P(X=x)=P(Y=x)
  #  So this method is useful when dealing with highly discrete data.
  #
  #  If KMS=T, use Kulinskaya, Morgenthaler and Staudte (2010)
  #   method for comparing binomials
  # Kulinskaya, E., Morgenthaler, S. and Staudte, R. (2010). 
  # Variance Stabilizing the Difference of two Binomial
  #  Proportions. {\em American Statistician, 64}, 
  #  350--356 DOI:10.1198/tast.2010.09096
  
  #  Otherwise use Storer and Kim.
  #
  #   ADJ.P=T means that critical p-value is adjusted to control FWE when the sample
  #   size is small (<50).
  #
  #
  #  Hochberg's method is used to determine critical p-values so that FWE=alpha
  #
  cl <- match.call()
  x=elimna(x)
  y=elimna(y)
  vals=sort(unique(c(x,y)))
  ncon=length(vals)
  n1=length(x)
  n2=length(y)
  p.values=NA
  adj=1
  cv=1
  if(!KMS){
    output=matrix(NA,ncol=6,nrow=length(vals))
    dimnames(output)=list(NULL,c("Value","p1.est","p2.est","p1-p2","p.value","p.crit"))
  }
  if(KMS){
    output=matrix(NA,ncol=8,nrow=length(vals))
    dimnames(output)=list(NULL,c("Value","p1.est","p2.est","p1-p2","ci.low","ci.up","p.value",
                                 "p.crit"))
  }
  for(i in 1:length(vals)){
    x1=sum(x==vals[i])
    y1=sum(y==vals[i])
    if(!KMS){
      output[i,5]=twobinom(x1,n1,y1,n2)$p.value
      output[i,2]=x1/n1
      output[i,3]=y1/n2
      output[i,1]=vals[i]
      output[i,4]=output[i,2]-output[i,3]
    }
    if(KMS){
      temp=bi2KMSv2(x1,n1,y1,n2)
      output[i,1]=vals[i]
      output[i,5]=temp$ci[1]
      output[i,6]=temp$ci[2]
      output[i,2]=x1/n1
      output[i,3]=y1/n2
      output[i,4]=output[i,2]-output[i,3]
      output[i,7]=temp$p.value
    }}
  # Determine adjusted  critical p-value using Hochberg method
  ncon=length(vals)
  dvec=alpha/c(1:ncon)
  if(ADJ.P){
    mn=max(c(n1,n2))
    cv=1
    if(ncon!=2){
      if(mn>50){
        cv=2-(mn-50)/50
        if(cv<1)cv=1
      }
      if(mn<=50)cv=2
    }
    if(KMS){
      flag=(output[,7]<=2*alpha)
      output[flag,8]=output[flag,8]/cv
    }
    if(!KMS){
      cv=1
      flag=(output[,5]<=2*alpha)
      if(min(c(n1,n2))<20 && n1!=n2 && ncon>=5)cv=2
      output[flag,5]=output[flag,5]/cv
    }}
  if(KMS){
    temp2=order(0-output[,7])
    output[temp2,8]=dvec
  }
  if(!KMS){
    temp2=order(0-output[,5])
    output[temp2,6]=dvec
  }
  output
  
  outtable <- data.frame(output)
  colnames(outtable)[4] <- "p1-p2"
  
  result <- list(partable = outtable, call = cl)
  class(result) <- "robtab"
  result
}

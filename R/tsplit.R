tsplit <- function(formula, random, data, tr = 0.2){
  
  if (missing(data)) {
    mf <- model.frame(formula)
  } else {
    mf <- model.frame(formula, data)
  }
  cl <- match.call()
  idRandom <- strsplit(strsplit(as.character(random)[2], "[>|]")[[1]][2], "[/]")[[1]]
  idvar <- gsub("\\s", "", idRandom[1])
  ranvar <- gsub("\\s", "", idRandom[2])
  depvar <- colnames(mf)[1]
  if(colnames(mf)[2] == ranvar) fixvar <- colnames(mf)[3] else fixvar <- colnames(mf)[2]
  
  K <- length(table(mf[,ranvar])) 
  J <- length(table(mf[,fixvar]))
  rowend <- table(mf[,fixvar])/K
  data$rowind <- c(1:rowend[1], 1:rowend[2])
  rowform <- as.formula(paste("rowind ~", fixvar, "+", ranvar))
  
  p <- J*K
  grp <- 1:p
  dataWide <- cast(data,  rowform, value = depvar)
  dataWide$rowind <- NULL
  
  data <- dataWide  
  
  if(is.data.frame(data))data=as.matrix(data)
  x<-data
  if(is.matrix(x)) {
    y <- list()
    for(j in 1:ncol(x))
      y[[j]] <- x[, j]
    data <- y
  }

  tmeans<-0
  h<-0
  v<-matrix(0,p,p)
  klow<-1-K
  kup<-0
  for (i in 1:p)tmeans[i]<-mean(data[[grp[i]]],tr,na.rm=TRUE)
  for (j in 1:J){
    h[j]<-length(data[[grp[j]]])-2*floor(tr*length(data[[grp[j]]]))
    #    h is the effective sample size for the jth level of factor A
    #   Use covmtrim to determine blocks of squared standard errors and
    #   covariances.
    klow<-klow+K
    kup<-kup+K
    sel<-c(klow:kup)
    v[sel,sel]<-covmtrim(data[grp[klow:kup]],tr)
  }
  ij<-matrix(c(rep(1,J)),1,J)
  ik<-matrix(c(rep(1,K)),1,K)
  jm1<-J-1
  cj<-diag(1,jm1,J)
  for (i in 1:jm1)cj[i,i+1]<-0-1
  km1<-K-1
  ck<-diag(1,km1,K)
  for (i in 1:km1)ck[i,i+1]<-0-1
  #  Do test for factor A
  cmat<-kron(cj,ik)  # Contrast matrix for factor A
  Qa<-johansp(cmat,tmeans,v,h,J,K)
  # Do test for factor B
  cmat<-kron(ij,ck)  # Contrast matrix for factor B
  Qb<-johansp(cmat,tmeans,v,h,J,K)
  # Do test for factor A by B interaction
  cmat<-kron(cj,ck)  # Contrast matrix for factor A by B
  Qab<-johansp(cmat,tmeans,v,h,J,K)
  
  result <- list(Qa=Qa$teststat, A.p.value=Qa$siglevel, Qb=Qb$teststat, B.p.value=Qb$siglevel, Qab=Qab$teststat, AB.p.value=Qab$siglevel, call = cl, varnames = c(depvar, fixvar, ranvar))
  class(result) <- c("t2way")
  result
}

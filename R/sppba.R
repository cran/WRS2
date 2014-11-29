sppba <- function(formula, random, data, est = "mom", avg = TRUE, nboot = 500, MDIS = FALSE){
  
  if (missing(data)) {
    mf <- model.frame(formula)
  } else {
    mf <- model.frame(formula, data)
  }
  cl <- match.call()
  idRandom <- strsplit(strsplit(as.character(random)[2], "[>|]")[[1]][2],"[/]")[[1]]
  
  print(head(mf))
  print(idRandom)
  
  idvar <- gsub("\\s", "", idRandom[1])
  ranvar <- gsub("\\s", "", idRandom[2])
  depvar <- colnames(mf)[1]
  if (colnames(mf)[2] == ranvar) fixvar <- colnames(mf)[3] else fixvar <- colnames(mf)[2]
  
  MC <- FALSE
  K <- length(table(mf[,ranvar])) 
  J <- length(table(mf[,fixvar]))
  JK <- J*K
  grp <- 1:JK
  est <- get(est)    ## converting string into function
   
  rowend <- table(mf[,ranvar])
  #data$rowind <- c(1:rowend[1], 1:rowend[2])
  mf$rowind <- sequence(rowend)
  rowend <- table(mf[, ranvar])
  rowend2 <- as.vector(sequence(rowend))
  maxrow <- max(rowend2)
  data.temp <- split(mf[,depvar], list(mf[,fixvar], mf[,ranvar]))
  data.temp1 <- data.temp[1:J]
  data.temp2 <- data.temp[(J+1):(J*K)]
  data.temp <- unlist(mapply(data.temp1, data.temp2, FUN = list, SIMPLIFY = FALSE), recursive = FALSE)
  
  data <- lapply(data.temp, function(xx) {
    yy <- rep(NA, maxrow)
    yy[1:length(xx)] <- xx
    return(yy)
  })
  x<-data
    
  jp<-1-K
  kv<-0
  kv2<-0
  for(j in 1:J){
    jp<-jp+K
    xmat<-matrix(NA,ncol=K,nrow=length(x[[jp]]))
    for(k in 1:K){
      kv<-kv+1
      xmat[,k]<-x[[kv]]
    }
    xmat<-elimna(xmat)
    for(k in 1:K){
      kv2<-kv2+1
      x[[kv2]]<-xmat[,k]
    }
  }
  xx<-x
  #if(SEED)set.seed(2) # set seed of random number generator so that
  #             results can be duplicated.
  # Next determine the n_j values
  nvec<-NA
  jp<-1-K
  for(j in 1:J){
    jp<-jp+K
    nvec[j]<-length(x[[jp]])
  }
  #
  # Now take bootstrap samples from jth level
  # of Factor A.
  #
  bloc<-matrix(NA,nrow=J,ncol=nboot)
  #print("Taking bootstrap samples. Please wait.")
  mvec<-NA
  ik<-0
  
  for(j in 1:J){
    #paste("Working on level ",j," of Factor A")
    x<-matrix(NA,nrow=nvec[j],ncol=K)
    #
    for(k in 1:K){
      ik<-ik+1
      x[,k]<-xx[[ik]]
      if(!avg)mvec[ik]<-est(xx[[ik]])
    }
    tempv<-apply(x,2,est)
    data<-matrix(sample(nvec[j],size=nvec[j]*nboot,replace=TRUE),nrow=nboot)
    bvec<-matrix(NA,ncol=K,nrow=nboot)
    for(k in 1:K){
      temp<-x[,k]
      bvec[,k]<-apply(data,1,rmanogsub,temp,est) # An nboot by K matrix
    }
    if(avg){
      mvec[j]<-mean(tempv)
      bloc[j,]<-apply(bvec,1,mean)
    }
    if(!avg){
      if(j==1)bloc<-bvec
      if(j>1)bloc<-cbind(bloc,bvec)
    }
  }
  
  
  if(avg){
    d<-(J^2-J)/2
    con<-matrix(0,J,d)
    id<-0
    Jm<-J-1
    for (j in 1:Jm){
      jp<-j+1
      for(k in jp:J){
        id<-id+1
        con[j,id]<-1
        con[k,id]<-0-1
      }}}
  if(!avg){
    MJK<-K*(J^2-J)/2 # NUMBER OF COMPARISONS
    JK<-J*K
    MJ<-(J^2-J)/2
    cont<-matrix(0,nrow=J,ncol=MJ)
    ic<-0
    for(j in 1:J){
      for(jj in 1:J){
        if(j<jj){
          ic<-ic+1
          cont[j,ic]<-1
          cont[jj,ic]<-0-1
        }}}
    tempv<-matrix(0,nrow=K-1,ncol=MJ)
    con1<-rbind(cont[1,],tempv)
    for(j in 2:J){
      con2<-rbind(cont[j,],tempv)
      con1<-rbind(con1,con2)
    }
    con<-con1
    if(K>1){
      for(k in 2:K){
        con1<-push(con1)
        con<-cbind(con,con1)
      }}}
  
  if(!avg)bcon<-t(con)%*%t(bloc) #C by nboot matrix
  if(avg)bcon<-t(con)%*%(bloc)
  tvec<-t(con)%*%mvec
  tvec<-tvec[,1]
  tempcen<-apply(bcon,1,mean)
  vecz<-rep(0,ncol(con))
  bcon<-t(bcon)
  temp=bcon
  for(ib in 1:nrow(temp))temp[ib,]=temp[ib,]-tempcen+tvec
  bcon<-rbind(bcon,vecz)
  if(!MDIS){
    if(!MC)dv=pdis(bcon,center=tvec)
    #if(MC)dv=pdisMC(bcon,center=tvec)
  }
  if(MDIS){
    smat<-var(temp)
    bcon<-rbind(bcon,vecz)
    chkrank<-qr(smat)$rank
    if(chkrank==ncol(smat))dv<-mahalanobis(bcon,tvec,smat)
    if(chkrank<ncol(smat)){
      smat<-ginv(smat)
      dv<-mahalanobis(bcon,tvec,smat,inverted=T)
    }}
  bplus<-nboot+1
  sig.level<-1-sum(dv[bplus]>=dv[1:nboot])/nboot
  
  result <- list(test = tvec, p.value = sig.level, call = cl)
  class(result) <- c("spp")
  result
}

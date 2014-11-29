sppbb <- function(formula, random, data, est = "mom", nboot = 500){
#
  if (missing(data)) {
    mf <- model.frame(formula)
  } else {
    mf <- model.frame(formula, data)
  }
  cl <- match.call()
  idRandom <- strsplit(strsplit(as.character(random)[2], "[>|]")[[1]][2],"[/]")[[1]]
  idvar <- gsub("\\s", "", idRandom[1])
  ranvar <- gsub("\\s", "", idRandom[2])
  depvar <- colnames(mf)[1]
  if (colnames(mf)[2] == ranvar) fixvar <- colnames(mf)[3] else fixvar <- colnames(mf)[2]
  K <- length(table(mf[, ranvar]))  ## number of repeated measurements
  J <- length(table(mf[, fixvar]))  ## number of levels
  rowend <- table(mf[, fixvar])/K
  rowend2 <- as.vector(sequence(rowend))
  maxrow <- max(rowend2)
  mf$rowind <- rowend2
  p <- J * K
  grp <- 1:p
  est <- get(est)    ## converting string into function
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
    }}
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
  # Now stack the data in an N by K matrix
  #
  x<-matrix(NA,nrow=nvec[1],ncol=K)
  #
  for(k in 1:K)x[,k]<-xx[[k]]
  kc<-K
  for(j in 2:J){
    temp<-matrix(NA,nrow=nvec[j],ncol=K)
    for(k in 1:K){
      kc<-kc+1
      temp[,k]<-xx[[kc]]
    }
    x<-rbind(x,temp)
  }
  # Now call function rmdzero to do the analysis
  temp<-rmdzero(x,est=est,nboot=nboot)
  result <- list(test = temp$center, p.value = temp$p.value, call = cl)
  class(result) <- c("spp")
  result
}

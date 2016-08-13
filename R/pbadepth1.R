pbadepth1<-function(x,est=onestep,con=0,alpha=.05,nboot=2000,grp=NA,op=1,allp=TRUE,
                   MM=FALSE,MC=FALSE,cop=3,SEED=TRUE,na.rm=FALSE,...){
  #
  #   Test the hypothesis that C linear contrasts all have a value of zero.
  #   By default, an M-estimator is used
  #
  #   Independent groups are assumed.
  #
  #   The data are assumed to be stored in x in list mode or in a matrix.
  #   If stored in list mode,
  #   x[[1]] contains the data for the first group, x[[2]] the data
  #   for the second group, etc. Length(x)=the number of groups = J, say.
  #   If stored in a matrix, columns correspond to groups.
  #
  #   By default, all pairwise differences are used, but contrasts
  #   can be specified with the argument con.
  #   The columns of con indicate the contrast coefficients.
  #   Con should have J rows, J=number of groups.
  #   For example, con[,1]=c(1,1,-1,-1,0,0) and con[,2]=c(,1,-1,0,0,1,-1)
  #   will test two contrasts: (1) the sum of the first
  #   two measures of location is
  #   equal to the sum of the second two, and (2) the difference between
  #   the first two is equal to the difference between the
  #   measures of location for groups 5 and 6.
  #
  #   The default number of bootstrap samples is nboot=2000
  #
  #   op controls how depth is measured
  #   op=1, Mahalanobis
  #   op=2, Mahalanobis based on MCD covariance matrix
  #   op=3, Projection distance
  #
  #   MC=TRUE, use a multicore processor when op=3
  #
  #   for arguments MM and cop, see pdis.
  #
  con<-as.matrix(con)
  if(is.matrix(x) || is.data.frame(x))x=listm(x)
  if(!is.list(x))stop("Data must be stored in list mode or in matrix mode.")
  if(!is.na(grp)){  # Only analyze specified groups.
    xx<-list()
    for(i in 1:length(grp))xx[[i]]<-x[[grp[i]]]
    x<-xx
  }
  J<-length(x)
  mvec<-NA
  nvec=NA
  for(j in 1:J){
    temp<-x[[j]]
    if(na.rm)temp<-temp[!is.na(temp)] # Remove missing values.
    x[[j]]<-temp
    mvec[j]<-est(temp,...)
    nvec[j]=length(temp)
  }
  Jm<-J-1
  d<-ifelse(con==0,(J^2-J)/2,ncol(con))
  if(sum(con^2)==0){
    if(allp){
      con<-matrix(0,J,d)
      id<-0
      for (j in 1:Jm){
        jp<-j+1
        for (k in jp:J){
          id<-id+1
          con[j,id]<-1
          con[k,id]<-0-1
        }}}
    if(!allp){
      con<-matrix(0,J,Jm)
      for (j in 1:Jm){
        jp<-j+1
        con[j,j]<-1
        con[jp,j]<-0-1
      }}}
  bvec<-matrix(NA,nrow=J,ncol=nboot)
  if(SEED)set.seed(2) # set seed of random number generator so that
  #             results can be duplicated.
  #print("Taking bootstrap samples. Please wait.")
  for(j in 1:J){
    #print(paste("Working on group ",j))
    data<-matrix(sample(x[[j]],size=length(x[[j]])*nboot,replace=TRUE),nrow=nboot)
    bvec[j,]<-apply(data,1,est,na.rm=na.rm,...) # J by nboot matrix, jth row contains
    #                          bootstrapped  estimates for jth group
  }
  chkna=sum(is.na(bvec))
  if(chkna>0){
    print("Bootstrap estimates of location could not be computed")
    print("This can occur when using an M-estimator")
    print("Might try est=tmean")
  }
  bcon<-t(con)%*%bvec #C by nboot matrix
  tvec<-t(con)%*%mvec
  tvec<-tvec[,1]
  tempcen<-apply(bcon,1,mean)
  vecz<-rep(0,ncol(con))
  bcon<-t(bcon)
  smat<-var(bcon-tempcen+tvec)
  temp<-bcon-tempcen+tvec
  bcon<-rbind(bcon,vecz)
  if(op==1)dv<-mahalanobis(bcon,tvec,smat)
  if(op==2){
    smat<-cov.mcd(temp)$cov
    dv<-mahalanobis(bcon,tvec,smat)
  }
  if(op==3){
    #print("Computing p-value. Might take a while with op=3")
    if(!MC)dv<-pdis(bcon,MM=MM,cop=cop)
    #if(MC)dv<-pdisMC(bcon,MM=MM,cop=cop)
  }
  bplus<-nboot+1
  sig.level<-1-sum(dv[bplus]>=dv[1:nboot])/nboot
  list(p.value=sig.level,psihat=tvec,con=con,n=nvec)
}

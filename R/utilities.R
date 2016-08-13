matl<-function(x){
  #
  # take data in list mode and store it in a matrix
  #
  J=length(x)
  nval=NA
  for(j in 1:J)nval[j]=length(x[[j]])
  temp<-matrix(NA,ncol=J,nrow=max(nval))
  for(j in 1:J)temp[1:nval[j],j]<-x[[j]]
  temp
}

list2mat=matl

list2vec<-function(x){
  if(!is.list(x))stop("x should have list mode")
  res=as.vector(matl(x))
  res
}

discANOVA.sub<-function(x){
  #
  #
  x=lapply(x,elimna)
  vals=lapply(x,unique)
  vals=sort(elimna(unique(list2vec(vals))))
  n=lapply(x,length)
  n=list2vec(n)
  K=length(vals)
  J=length(x)
  C1=matrix(0,nrow=K,ncol=J)
  for(j in 1:J){
    for(i in 1:K){
      C1[i,j]=C1[i,j]+sum(x[[j]]==vals[i])
    }
    C1[,j]=C1[,j]/n[j]
  }
  test=0
  for(i in 1:K)test=test+var(C1[i,])
  list(test=test,C1=C1)
}

modgen<-function(p,adz=FALSE){
  #
  # Used by regpre to generate all models
  # p=number of predictors
  # adz=T, will add the model where only a measure
  # of location is used.
  #
  #
  model<-list()
  if(p>5)stop("Current version is limited to 5 predictors")
  if(p==1)model[[1]]<-1
  if(p==2){
    model[[1]]<-1
    model[[2]]<-2
    model[[3]]<-c(1,2)
  }
  if(p==3){
    for(i in 1:3)model[[i]]<-i
    model[[4]]<-c(1,2)
    model[[5]]<-c(1,3)
    model[[6]]<-c(2,3)
    model[[7]]<-c(1,2,3)
  }
  if(p==4){
    for(i in 1:4)model[[i]]<-i
    model[[5]]<-c(1,2)
    model[[6]]<-c(1,3)
    model[[7]]<-c(1,4)
    model[[8]]<-c(2,3)
    model[[9]]<-c(2,4)
    model[[10]]<-c(3,4)
    model[[11]]<-c(1,2,3)
    model[[12]]<-c(1,2,4)
    model[[13]]<-c(1,3,4)
    model[[14]]<-c(2,3,4)
    model[[15]]<-c(1,2,3,4)
  }
  if(p==5){
    for(i in 1:5)model[[i]]<-i
    model[[6]]<-c(1,2)
    model[[7]]<-c(1,3)
    model[[8]]<-c(1,4)
    model[[9]]<-c(1,5)
    model[[10]]<-c(2,3)
    model[[11]]<-c(2,4)
    model[[12]]<-c(2,5)
    model[[13]]<-c(3,4)
    model[[14]]<-c(3,5)
    model[[15]]<-c(4,5)
    model[[16]]<-c(1,2,3)
    model[[17]]<-c(1,2,4)
    model[[18]]<-c(1,2,5)
    model[[19]]<-c(1,3,4)
    model[[20]]<-c(1,3,5)
    model[[21]]<-c(1,4,5)
    model[[22]]<-c(2,3,4)
    model[[23]]<-c(2,3,5)
    model[[24]]<-c(2,4,5)
    model[[25]]<-c(3,4,5)
    model[[26]]<-c(1,2,3,4)
    model[[27]]<-c(1,2,3,5)
    model[[28]]<-c(1,2,4,5)
    model[[29]]<-c(1,3,4,5)
    model[[30]]<-c(2,3,4,5)
    model[[31]]<-c(1,2,3,4,5)
  }
  if(adz){
    ic<-length(model)+1
    model[[ic]]<-0
  }
  model
}

pb2gen1<-function(x,y,alpha=.05,nboot=2000,est=onestep,SEED=TRUE,pr=FALSE,...){
  #
  #   Compute a bootstrap confidence interval for the
  #   the difference between any two parameters corresponding to
  #   independent groups.
  #   By default, M-estimators are compared.
  #   Setting est=mean, for example, will result in a percentile
  #   bootstrap confidence interval for the difference between means.
  #   Setting est=onestep will compare M-estimators of location.
  #   The default number of bootstrap samples is nboot=2000
  #
  x<-x[!is.na(x)] # Remove any missing values in x
  y<-y[!is.na(y)] # Remove any missing values in y
  if(SEED)set.seed(2) # set seed of random number generator so that
  #             results can be duplicated.
  datax<-matrix(sample(x,size=length(x)*nboot,replace=TRUE),nrow=nboot)
  datay<-matrix(sample(y,size=length(y)*nboot,replace=TRUE),nrow=nboot)
  bvecx<-apply(datax,1,est,...)
  bvecy<-apply(datay,1,est,...)
  bvec<-sort(bvecx-bvecy)
  low<-round((alpha/2)*nboot)+1
  up<-nboot-low
  temp<-sum(bvec<0)/nboot+sum(bvec==0)/(2*nboot)
  sig.level<-2*(min(temp,1-temp))
  se<-var(bvec)
  list(est.1=est(x,...),est.2=est(y,...),est.dif=est(x,...)-est(y,...),ci=c(bvec[low],bvec[up]),p.value=sig.level,sq.se=se,n1=length(x),n2=length(y))
}

bootdpci<-function(x,y,est=onestep,nboot=NA,alpha=.05,plotit=TRUE,dif=TRUE,BA=FALSE,SR=TRUE,...){
  #
  #   Use percentile bootstrap method,
  #   compute a .95 confidence interval for the difference between
  #   a measure of location or scale
  #   when comparing two dependent groups.
  #   By default, a one-step M-estimator (with Huber's psi) is used.
  #   If, for example, it is desired to use a fully iterated
  #   M-estimator, use fun=mest when calling this function.
  #
  okay=FALSE
  if(identical(est,onestep))okay=TRUE
  if(identical(est,mom))okay=TRUE
  if(!okay)SR=FALSE
  output<-rmmcppb(x,y,est=est,nboot=nboot,alpha=alpha,SR=SR,
                  plotit=plotit,dif=dif,BA=BA,...)$output
  list(output=output)
}



rmmcppb<-function(x,y=NULL,alpha=.05,
                  con=0,est=onestep,plotit=TRUE,dif=TRUE,grp=NA,nboot=NA,BA=FALSE,hoch=FALSE,xlab="Group 1",ylab="Group 2",pr=TRUE,SEED=TRUE,SR=FALSE,...){
  #
  #   Use a percentile bootstrap method to  compare dependent groups.
  #   By default,
  #   compute a .95 confidence interval for all linear contrasts
  #   specified by con, a J-by-C matrix, where  C is the number of
  #   contrasts to be tested, and the columns of con are the
  #   contrast coefficients.
  #   If con is not specified, all pairwise comparisons are done.
  #
  #   If est=onestep or mom, method SR (see my book on robust methods)
  #   is used to control the probability of at least one Type I error.
  #
  #   Otherwise, Hochberg is used.
  #
  #   dif=T indicates that difference scores are to be used
  #   dif=F indicates that measure of location associated with
  #   marginal distributions are used instead.
  #
  #   nboot is the bootstrap sample size. If not specified, a value will
  #   be chosen depending on the number of contrasts there are.
  #
  #   x can be an n by J matrix or it can have list mode
  #   for two groups, data for second group can be put in y
  #   otherwise, assume x is a matrix (n by J) or has list mode.
  #
  #   A sequentially rejective method is used to control alpha using method SR.
  #
  #   Argument BA: When using dif=F, BA=T uses a correction term
  #   when computing a p-value.
  #
  if(hoch)SR=FALSE #Assume Hochberg if hoch=TRUE even if SR=TRUE
  if(SR){
    okay=FALSE
    if(identical(est,onestep))okay=TRUE
    if(identical(est,mom))okay=TRUE
    SR=okay # 'Only use method SR (argument SR=T) when est=onestep or mom
  }
  if(dif){
    if(pr)print("dif=T, so analysis is done on difference scores")
    temp<-rmmcppbd(x,y=y,alpha=.05,con=con,est,plotit=plotit,grp=grp,nboot=nboot,
                   hoch=TRUE,...)
    output<-temp$output
    con<-temp$con
  }
  if(!dif){
    if(pr){
      print("dif=F, so analysis is done on marginal distributions")
      if(!BA){
        if(identical(est,onestep))print("With M-estimator or MOM, suggest using BA=T and hoch=T")
        if(identical(est,mom))print("With M-estimator or MOM, suggest using BA=T and hoch=T")
      }}
    if(!is.null(y[1]))x<-cbind(x,y)
    if(!is.list(x) && !is.matrix(x))stop("Data must be stored in a matrix or in list mode.")
    if(is.list(x)){
      if(is.matrix(con)){
        if(length(x)!=nrow(con))stop("The number of rows in con is not equal to the number of groups.")
      }}
    if(is.list(x)){
      # put the data in an n by J matrix
      mat<-matl(x)
    }
    if(is.matrix(x) && is.matrix(con)){
      if(ncol(x)!=nrow(con))stop("The number of rows in con is not equal to the number of groups.")
      mat<-x
    }
    if(is.matrix(x))mat<-x
    if(!is.na(sum(grp)))mat<-mat[,grp]
    mat<-elimna(mat) # Remove rows with missing values.
    x<-mat
    J<-ncol(mat)
    xcen<-x
    for(j in 1:J)xcen[,j]<-x[,j]-est(x[,j],...)
    Jm<-J-1
    if(sum(con^2)==0){
      d<-(J^2-J)/2
      con<-matrix(0,J,d)
      id<-0
      for (j in 1:Jm){
        jp<-j+1
        for (k in jp:J){
          id<-id+1
          con[j,id]<-1
          con[k,id]<-0-1
        }}}
    d<-ncol(con)
    if(is.na(nboot)){
      if(d<=4)nboot<-1000
      if(d>4)nboot<-5000
    }
    n<-nrow(mat)
    crit.vec<-alpha/c(1:d)
    connum<-ncol(con)
    if(SEED)set.seed(2) # set seed of random number generator so that
    #             results can be duplicated.
    xbars<-apply(mat,2,est,...)
    psidat<-NA
    for (ic in 1:connum)psidat[ic]<-sum(con[,ic]*xbars)
    psihat<-matrix(0,connum,nboot)
    psihatcen<-matrix(0,connum,nboot)
    bvec<-matrix(NA,ncol=J,nrow=nboot)
    bveccen<-matrix(NA,ncol=J,nrow=nboot)
    if(pr)print("Taking bootstrap samples. Please wait.")
    data<-matrix(sample(n,size=n*nboot,replace=TRUE),nrow=nboot)
    for(ib in 1:nboot){
      bvec[ib,]<-apply(x[data[ib,],],2,est,...)
      bveccen[ib,]<-apply(xcen[data[ib,],],2,est,...)
    }
    #
    # Now have an nboot by J matrix of bootstrap values.
    #
    test<-1
    bias<-NA
    for (ic in 1:connum){
      psihat[ic,]<-apply(bvec,1,bptdpsi,con[,ic])
      psihatcen[ic,]<-apply(bveccen,1,bptdpsi,con[,ic])
      bias[ic]<-sum((psihatcen[ic,]>0))/nboot-.5
      ptemp<-(sum(psihat[ic,]>0)+.5*sum(psihat[ic,]==0))/nboot
      if(BA)test[ic]<-ptemp-.1*bias[ic]
      if(!BA)test[ic]<-ptemp
      test[ic]<-min(test[ic],1-test[ic])
      test[ic]<-max(test[ic],0)  # bias corrected might be less than zero
    }
    test<-2*test
    ncon<-ncol(con)
    dvec<-alpha/c(1:ncon) # Assume Hochberg unless specified otherwise
    if(SR){
      if(alpha==.05){
        dvec<-c(.025,.025,.0169,.0127,.0102,.00851,.0073,.00639,.00568,.00511)
        dvecba<-c(.05,.025,.0169,.0127,.0102,.00851,.0073,.00639,.00568,.00511)
        if(ncon > 10){
          avec<-.05/c(11:ncon)
          dvec<-c(dvec,avec)
        }}
      if(alpha==.01){
        dvec<-c(.005,.005,.00334,.00251,.00201,.00167,.00143,.00126,.00112,.00101)
        dvecba<-c(.01,.005,.00334,.00251,.00201,.00167,.00143,.00126,.00112,.00101)
        if(ncon > 10){
          avec<-.01/c(11:ncon)
          dvec<-c(dvec,avec)
        }}
      if(alpha != .05 && alpha != .01){
        dvec<-alpha/c(1:ncon)
        dvecba<-dvec
        dvec[2]<-alpha
      }}
    if(hoch)dvec<-alpha/c(1:ncon)
    dvecba<-dvec
    if(plotit && ncol(bvec)==2){
      z<-c(0,0)
      one<-c(1,1)
      plot(rbind(bvec,z,one),xlab=xlab,ylab=ylab,type="n")
      points(bvec)
      totv<-apply(x,2,est,...)
      cmat<-var(bvec)
      dis<-mahalanobis(bvec,totv,cmat)
      temp.dis<-order(dis)
      ic<-round((1-alpha)*nboot)
      xx<-bvec[temp.dis[1:ic],]
      xord<-order(xx[,1])
      xx<-xx[xord,]
      temp<-chull(xx)
      lines(xx[temp,])
      lines(xx[c(temp[1],temp[length(temp)]),])
      abline(0,1)
    }
    temp2<-order(0-test)
    ncon<-ncol(con)
    zvec<-dvec[1:ncon]
    if(BA)zvec<-dvecba[1:ncon]
    sigvec<-(test[temp2]>=zvec)
    output<-matrix(0,connum,6)
    dimnames(output)<-list(NULL,c("con.num","psihat","p.value","p.sig","ci.lower","ci.upper"))
    tmeans<-apply(mat,2,est,...)
    psi<-1
    output[temp2,4]<-zvec
    for (ic in 1:ncol(con)){
      output[ic,2]<-sum(con[,ic]*tmeans)
      output[ic,1]<-ic
      output[ic,3]<-test[ic]
      temp<-sort(psihat[ic,])
      #icl<-round(output[ic,4]*nboot/2)+1
      icl<-round(alpha*nboot/2)+1
      icu<-nboot-(icl-1)
      output[ic,5]<-temp[icl]
      output[ic,6]<-temp[icu]
    }
  }
  num.sig<-sum(output[,3]<=output[,4])
  list(output=output,con=con,num.sig=num.sig)
}

bptdpsi<-function(x,con){
  # Used by bptd to compute bootstrap psihat values
  #
  bptdpsi<-sum(con*x)
  bptdpsi
}

rmmismcp<-function(x,y=NA,alpha=.05,con=0,est=tmean,plotit=TRUE,grp=NA,nboot=500,
                   SEED=TRUE,xlab="Group 1",ylab="Group 2",pr=FALSE,...){
  #
  #   Use a percentile bootstrap method to  compare  marginal measures of location for dependent groups.
  #   Missing values are allowed; vectors of observations that contain
  #   missing values are not simply removed as done by rmmcppb.
  #   Only marginal measures of location are compared,
  #   The function computes a .95 confidence interval for all linear contrasts
  #   specified by con, a J by C matrix, where  C is the number of
  #   contrasts to be tested, and the columns of con are the
  #   contrast coefficients.
  #   If con is not specified, all pairwise comparisons are done.
  #
  #   By default, a 20% trimmed is used and a sequentially rejective method
  #   is used to control the probability of at least one Type I error.
  #
  #   nboot is the bootstrap sample size.
  #
  #   x can be an n by J matrix or it can have list mode
  #   for two groups, data for second group can be put in y
  #   otherwise, assume x is a matrix (n by J) or has list mode.
  #
  #
  if(!is.na(y[1]))x<-cbind(x,y)
  if(is.list(x))x=matl(x)
  if(!is.list(x) && !is.matrix(x))stop("Data must be stored in a matrix or in list mode.")
  if(is.list(x)){
    if(is.matrix(con)){
      if(length(x)!=nrow(con))stop("The number of rows in con is not equal to the number of groups.")
    }}
  if(is.list(x)){
    # put the data in an n by J matrix
    mat<-matl(x)
  }
  if(is.matrix(x) && is.matrix(con)){
    if(ncol(x)!=nrow(con))stop("The number of rows in con is not equal to the number of groups.")
    mat<-x
  }
  J<-ncol(x)
  Jm<-J-1
  flag.con=F
  if(sum(con^2)==0){
    flag.con=T
    d<-(J^2-J)/2
    con<-matrix(0,J,d)
    id<-0
    for (j in 1:Jm){
      jp<-j+1
      for (k in jp:J){
        id<-id+1
        con[j,id]<-1
        con[k,id]<-0-1
      }}}
  d<-ncol(con)
  n<-nrow(x)
  crit.vec<-alpha/c(1:d)
  connum<-ncol(con)
  if(SEED)set.seed(2) # set seed of random number generator so that
  #             results can be duplicated.
  xbars<-apply(x,2,est,na.rm=TRUE)
  psidat<-NA
  bveccen<-matrix(NA,ncol=J,nrow=nboot)
  for (ic in 1:connum)psidat[ic]<-sum(con[,ic]*xbars)
  psihat<-matrix(0,connum,nboot)
  psihatcen<-matrix(0,connum,nboot)
  bvec<-matrix(NA,ncol=J,nrow=nboot)
  if(pr)print("Taking bootstrap samples. Please wait.")
  data<-matrix(sample(n,size=n*nboot,replace=TRUE),nrow=nboot)
  for(ib in 1:nboot){
    bvec[ib,]<-apply(x[data[ib,],],2,est,na.rm=TRUE,...)
  }
  #
  # Now have an nboot by J matrix of bootstrap measures of location.
  #
  test<-1
  for (ic in 1:connum){
    for(ib in 1:nboot){
      psihat[ic,ib]=sum(con[,ic]*bvec[ib,])
    }
    matcon=c(0,psihat[ic,])
    dis=mean((psihat[ic,]<0))+.5*mean((psihat[ic,]==0))
    test[ic]<-2*min(c(dis,1-dis)) # the p-value
  }
  ncon<-ncol(con)
  dvec<-alpha/c(1:ncon)
  if(plotit && ncol(bvec)==2){
    z<-c(0,0)
    one<-c(1,1)
    plot(rbind(bvec,z,one),xlab=xlab,ylab=ylab,type="n")
    points(bvec)
    totv<-apply(x,2,est,na.rm=TRUE,...)
    cmat<-var(bvec)
    dis<-mahalanobis(bvec,totv,cmat)
    temp.dis<-order(dis)
    ic<-round((1-alpha)*nboot)
    xx<-bvec[temp.dis[1:ic],]
    xord<-order(xx[,1])
    xx<-xx[xord,]
    temp<-chull(xx)
    lines(xx[temp,])
    lines(xx[c(temp[1],temp[length(temp)]),])
    abline(0,1)
  }
  temp2<-order(0-test)
  ncon<-ncol(con)
  zvec<-dvec[1:ncon]
  sigvec<-(test[temp2]>=zvec)
  output<-matrix(0,connum,6)
  dimnames(output)<-list(NULL,c("con.num","psihat","p.value",
                                "crit.sig","ci.lower","ci.upper"))
  tmeans<-apply(x,2,est,na.rm=TRUE,...)
  psi<-1
  output[temp2,4]<-zvec
  for (ic in 1:ncol(con)){
    output[ic,2]<-sum(con[,ic]*tmeans)
    output[ic,1]<-ic
    output[ic,3]<-test[ic]
    temp<-sort(psihat[ic,])
    icl<-round(output[ic,4]*nboot/2)+1
    icu<-nboot-(icl-1)
    output[ic,5]<-temp[icl]
    output[ic,6]<-temp[icu]
  }
  if(!flag.con){
  }
  if(flag.con){
    CC=(J^2-J)/2
    test<-matrix(NA,CC,7)
    dimnames(test)<-list(NULL,c("Group","Group","psi.hat","p.value","p.crit",
                                "ci.low","ci.upper"))
    jcom<-0
    for (j in 1:J){
      for (k in 1:J){
        if (j < k){
          jcom<-jcom+1
          test[jcom,1]=j
          test[jcom,2]=k
          test[jcom,3:5]=output[jcom,2:4]
          test[jcom,6:7]=output[jcom,5:6]
          con=NULL
        }}}}
  if(!flag.con)test=output
  #num.sig<-sum(output[,4]<=output[,5])
  if(flag.con)num.sig<-sum(test[,4]<=test[,5])
  if(!flag.con)num.sig<-sum(test[,3]<=test[,4])
  list(output=test,con=con,num.sig=num.sig)
}


rmmcppbd<-function(x,y=NULL,alpha=.05,con=0,est=onestep,plotit=TRUE,grp=NA,nboot=NA,
                   hoch=TRUE,SEED=TRUE,...){
  #
  #   Use a percentile bootstrap method to  compare dependent groups
  #   based on difference scores.
  #   By default,
  #   compute a .95 confidence interval for all linear contrasts
  #   specified by con, a J by C matrix, where  C is the number of
  #   contrasts to be tested, and the columns of con are the
  #   contrast coefficients.
  #   If con is not specified, all pairwise comparisons are done.
  #
  #   By default, one-step M-estimator is used
  #    and a sequentially rejective method
  #   is used to control the probability of at least one Type I error.
  #
  #   nboot is the bootstrap sample size. If not specified, a value will
  #   be chosen depending on the number of contrasts there are.
  #
  #   x can be an n by J matrix or it can have list mode
  #   for two groups, data for second group can be put in y
  #   otherwise, assume x is a matrix (n by J) or has list mode.
  #
  #   A sequentially rejective method is used to control alpha.
  #   If n>=80, hochberg's method is used.
  #
  if(!is.null(y[1]))x<-cbind(x,y)
  if(!is.list(x) && !is.matrix(x))stop("Data must be stored in a matrix or in list mode.")
  if(is.list(x)){
    if(is.matrix(con)){
      if(length(x)!=nrow(con))stop("The number of rows in con is not equal to the number of groups.")
    }}
  if(is.list(x)){
    # put the data in an n by J matrix
    mat<-matl(x)
  }
  if(is.matrix(x) && is.matrix(con)){
    if(ncol(x)!=nrow(con))stop("The number of rows in con is not equal to the number of groups.")
    mat<-x
  }
  if(is.matrix(x))mat<-x
  if(!is.na(sum(grp)))mat<-mat[,grp]
  x<-mat
  mat<-elimna(mat) # Remove rows with missing values.
  x<-mat
  J<-ncol(mat)
  n=nrow(mat)
  if(n>=80)hoch=T
  Jm<-J-1
  if(sum(con^2)==0){
    d<-(J^2-J)/2
    con<-matrix(0,J,d)
    id<-0
    for (j in 1:Jm){
      jp<-j+1
      for (k in jp:J){
        id<-id+1
        con[j,id]<-1
        con[k,id]<-0-1
      }}}
  d<-ncol(con)
  if(is.na(nboot)){
    nboot<-5000
    if(d<=10)nboot<-3000
    if(d<=6)nboot<-2000
    if(d<=4)nboot<-1000
  }
  n<-nrow(mat)
  crit.vec<-alpha/c(1:d)
  connum<-ncol(con)
  # Create set of differences based on contrast coefficients
  xx<-x%*%con
  xx<-as.matrix(xx)
  if(SEED)set.seed(2) # set seed of random number generator so that
  #             results can be duplicated.
  psihat<-matrix(0,connum,nboot)
  bvec<-matrix(NA,ncol=connum,nrow=nboot)
  data<-matrix(sample(n,size=n*nboot,replace=TRUE),nrow=nboot)
  # data is an nboot by n matrix
  if(ncol(xx)==1){
    for(ib in 1:nboot)psihat[1,ib]<-est(xx[data[ib,]],...)
  }
  if(ncol(xx)>1){
    for(ib in 1:nboot)psihat[,ib]<-apply(xx[data[ib,],],2,est,...)
  }
  #
  # Now have an nboot by connum matrix of bootstrap values.
  #
  test<-1
  for (ic in 1:connum){
    test[ic]<-(sum(psihat[ic,]>0)+.5*sum(psihat[ic,]==0))/nboot
    test[ic]<-min(test[ic],1-test[ic])
  }
  test<-2*test
  ncon<-ncol(con)
  if(alpha==.05){
    dvec<-c(.025,.025,.0169,.0127,.0102,.00851,.0073,.00639,.00568,.00511)
    if(ncon > 10){
      avec<-.05/c(11:ncon)
      dvec<-c(dvec,avec)
    }}
  if(alpha==.01){
    dvec<-c(.005,.005,.00334,.00251,.00201,.00167,.00143,.00126,.00112,.00101)
    if(ncon > 10){
      avec<-.01/c(11:ncon)
      dvec<-c(dvec,avec)
    }}
  if(alpha != .05 && alpha != .01){
    dvec<-alpha/c(1:ncon)
    dvec[2]<-alpha/2
  }
  if(hoch)dvec<-alpha/(2*c(1:ncon))
  dvec<-2*dvec
  if(plotit && connum==1){
    plot(c(psihat[1,],0),xlab="",ylab="Est. Difference")
    points(psihat[1,])
    abline(0,0)
  }
  temp2<-order(0-test)
  ncon<-ncol(con)
  zvec<-dvec[1:ncon]
  sigvec<-(test[temp2]>=zvec)
  output<-matrix(0,connum,6)
  dimnames(output)<-list(NULL,c("con.num","psihat","p.value","p.crit","ci.lower","ci.upper"))
  tmeans<-apply(xx,2,est,...)
  psi<-1
  icl<-round(dvec[ncon]*nboot/2)+1
  icu<-nboot-icl-1
  for (ic in 1:ncol(con)){
    output[ic,2]<-tmeans[ic]
    output[ic,1]<-ic
    output[ic,3]<-test[ic]
    output[temp2,4]<-zvec
    temp<-sort(psihat[ic,])
    output[ic,5]<-temp[icl]
    output[ic,6]<-temp[icu]
  }
  num.sig<-sum(output[,3]<=output[,4])
  list(output=output,con=con,num.sig=num.sig)
}

HuberTun=function(kappa,p){
  # Tunning parameter when use Huber type weight
  #------------------------------------------------------------
  # Input:
  #kappa: the proportion of cases to be controlled
  #p: the number of variables
  # Output
  # r: the critical value of Mahalalanobis distance, as defined in (20)
  # tau: the constant to make the robust estimator of Sigma to be unbiased, as defined in (20)
  
  prob=1-kappa
  chip=qchisq(prob,p)
  r=sqrt(chip)
  tau=(p*pchisq(chip,p+2)+ chip*(1-prob))/p	
  Results=list(r=r,tau=tau)
  return(Results)	
}


robEst=function(Z,r,tau,ep){
  
  p=ncol(Z)
  n=nrow(Z)
  # Starting values 	
  mu0=MeanCov(Z)$zbar
  Sigma0=MeanCov(Z)$S
  Sigin=solve(Sigma0)
  
  diverg=0 # convergence flag
  
  for (k in 1:200) {   	
    sumu1=0
    mu=matrix(0,p,1)
    Sigma=matrix(0,p,p)
    d=rep(NA,n)
    u1=rep(NA,n)
    u2=rep(NA,n)
    
    for (i in 1:n) {			zi=Z[i,]
    zi0=zi-mu0
    di2=t(zi0)%*%Sigin%*%zi0
    di=as.numeric(sqrt(di2))
    d[i]=di
    
    #get u1i,u2i
    if (di<=r) {
      u1i=1.0
      u2i=1.0/tau
    }else {
      u1i=r/di
      u2i=u1i^2/tau
    }
    u1[i]=u1i
    u2[i]=u2i
    
    sumu1=sumu1+u1i
    mu=mu+u1i*zi
    Sigma=Sigma+u2i*zi0%*%t(zi0)
    
    } # end of loop i
    
    mu1=mu/sumu1
    Sigma1=Sigma/n
    Sigdif=Sigma1-Sigma0
    dt=sum(Sigdif^2)
    
    mu0=mu1
    Sigma0=Sigma1
    Sigin=solve(Sigma0)
    if (dt<ep) {break}
    
  } # end of loop k
  
  
  if (k==200) {
    diverg=1
    mu0=rep(0,p)
    sigma0=matrix(NA,p,p)
    
  }
  
  theta=MLEst(Sigma0)
  
  Results=list(mu=mu0,Sigma=Sigma0,theta=theta,d=d,u1=u1,u2=u2,diverg=diverg)
  return(Results)
}

SErob=function(Z,mu,Sigma,theta,d,r,tau){
  n=nrow(Z)
  p=ncol(Z)
  ps=p*(p+1)/2
  q=6
  Dup=Dp(p)
  DupPlus=solve(t(Dup)%*%Dup)%*%t(Dup)
  
  InvSigma=solve(Sigma)
  sigma=vech(Sigma)
  W=0.5*t(Dup)%*%(InvSigma%x%InvSigma)%*%Dup
  
  Zr=matrix(NA,n,p) # robustly transformed data
  A=matrix(0,p+q,p+q)
  B=matrix(0,p+q,p+q)
  sumg=rep(0,p+q)
  
  for (i in 1:n) {
    zi=Z[i,]
    zi0=zi-mu
    di=d[i]
    
    if (di<=r) {
      u1i=1.0
      u2i=1.0/tau
      du1i=0
      du2i=0
    }else {
      u1i=r/di
      u2i=u1i^2/tau
      du1i=-r/di^2
      du2i=-2*r^2/tau/di^3
    }
    
    #robust transformed data
    Zr[i,]=sqrt(u2i)*t(zi0)
    
    ####	gi
    
    g1i=u1i*zi0	# defined in (24)
    vTi=vech(zi0%*%t(zi0))
    g2i=u2i*vTi-sigma	# defined in (25)
    gi=rbind(g1i,g2i)
    sumg=gi+sumg
    
    B=B+gi%*%t(gi)
    
    ####	gdoti
    
    #	derivatives of di with respect to mu and sigma
    ddmu=-1/di*t(zi0)%*%InvSigma
    ddsigma=-t(vTi)%*%W/di
    
    #	
    dg1imu=-u1i*diag(p)+du1i*zi0%*%ddmu
    dg1isigma=du1i*zi0%*%ddsigma
    dg2imu=-u2i*DupPlus%*%(zi0%x%diag(p)+diag(p)%x%zi0)+du2i*vTi%*%ddmu
    dg2isigma=du2i*vTi%*%ddsigma-diag(q)
    
    dgi=rbind(cbind(dg1imu,dg1isigma),cbind(dg2imu,dg2isigma))
    A=A+dgi
  } # end of loop n
  
  A=-1*A/n
  B=B/n
  invA=solve(A)
  OmegaSW=invA%*%B%*%t(invA)
  OmegaSW=OmegaSW[(p+1):(p+q),(p+1):(p+q)]
  
  
  SEsw=getSE(theta,OmegaSW,n)
  SEinf=SEML(Zr,theta)$inf
  
  Results=list(inf=SEinf,sand=SEsw,Zr=Zr)  	
  return(Results)
  
}

MeanCov=function(Z){
  n=nrow(Z)
  p=ncol(Z)
  
  zbar=t(Z)%*%matrix(1/n,n,1)
  S=t(Z)%*%(diag(n)-matrix(1/n,n,n))%*%Z/n
  
  Results=list(zbar=zbar,S=S)
  return(Results)
}

#-----------------------------------------------------------------------
# Obtaining normal-theory MLE of parameters in the mediation model
#-----------------------------------------------------------------------
# Input:
# S: sample covariance
# Output:
# thetaMLE: normal-theory MLE of theta. theta is defined in the subsection: MLEs of a,b, and c

MLEst=function(S){
  ahat=S[1,2]/S[1,1]
  vx=S[1,1]
  # M on X
  Sxx=S[1:2,1:2]
  sxy=S[1:2,3]
  vem=S[2,2]-S[2,1]*S[1,2]/S[1,1]
  
  # Y on X and M
  invSxx=solve(Sxx)
  beta.v=invSxx%*%sxy # chat, bhat
  vey=S[3,3]-t(sxy)%*%invSxx%*%sxy
  thetaMLE=c(ahat,beta.v[2],beta.v[1],vx,vem,vey)
  return(thetaMLE)
}

Dp=function(p){
  p2=p*p
  ps=p*(p+1)/2
  Dup=matrix(0,p2,ps)
  count=0
  for (j in 1:p){
    for (i in j:p){
      count=count+1
      if (i==j){
        Dup[(j-1)*p+j, count]=1
      }else{
        Dup[(j-1)*p+i, count]=1
        Dup[(i-1)*p+j, count]=1
      }
    }
  }	
  
  return(Dup)	
}

vech=function(A){
  l=0
  p=nrow(A)
  ps=p*(p+1)/2
  vA=matrix(0,ps,1)
  for (i in 1:p) {
    for (j in i:p) {
      l=l+1
      vA[l,1]=A[j,i]			
    }	
  }
  
  return(vA)	
}


getSE=function(theta,Omega,n){
  
  hdot=gethdot(theta)
  COV=hdot%*%(Omega/n)%*%t(hdot)  # delta method
  se.v=sqrt(diag(COV)) # se.v of theta
  
  a=theta[1]
  b=theta[2]
  SobelSE=sqrt(a^2*COV[2,2]+b^2*COV[1,1])
  
  se.v=c(se.v,SobelSE) # including Sobel SE
  
  return(se.v)
  
}
gethdot=function(theta){
  
  p=3
  ps=p*(p+1)/2
  q=ps
  
  a=theta[1]
  b=theta[2]
  c=theta[3]
  #ab=a*b
  vx=theta[4]
  vem=theta[5]
  vey=theta[6]
  
  sigmadot=matrix(NA,ps,q)
  sigmadot[1,]=c(0,0,0,1,0,0)
  sigmadot[2,]=c(vx,0,0,a,0,0)
  sigmadot[3,]=c(b*vx,a*vx,vx,a*b+c,0,0)
  sigmadot[4,]=c(2*a*vx,0,0,a^2,1,0)
  sigmadot[5,]=c((2*a*b+c)*vx,a^2*vx+vem,a*vx,a^2*b+a*c,b,0)
  sigmadot[6,]=c((2*b*c+2*a*b^2)*vx,(2*c*a+2*a^2*b)*vx+2*b*vem,(2*a*b+2*c)*vx,c^2+2*c*a*b+a^2*b^2,b^2,1)
  
  hdot=solve(sigmadot)
  
  return(hdot)
}
SEML=function(Z,thetaMLE){
  n=nrow(Z)
  p=ncol(Z)
  ps=p*(p+1)/2
  q=ps
  zbar=MeanCov(Z)$zbar
  S=MeanCov(Z)$S
  Dup=Dp(p)
  InvS=solve(S)
  W=0.5*t(Dup)%*%(InvS%x%InvS)%*%Dup
  OmegaInf=solve(W)  # only about sigma, not mu
  
  
  # Sandwich-type Omega
  S12=matrix(0,p,ps)
  S22=matrix(0,ps,ps)
  
  for (i in 1:n){
    zi0=Z[i,]-zbar
    difi=zi0%*%t(zi0)-S
    vdifi=vech(difi)   		
    S12=S12+zi0%*%t(vdifi)	
    S22=S22+vdifi%*%t(vdifi)		
  }
  
  OmegaSW=S22/n # only about sigma, not mu
  
  SEinf=getSE(thetaMLE,OmegaInf,n)
  SEsw=getSE(thetaMLE,OmegaSW,n)
  
  Results=list(inf=SEinf,sand=SEsw)
  return(Results)
  
}
BCI=function(Z,Zr,ab=NULL,abH,B,level){
  p=ncol(Z)
  n=nrow(Z)
  #	abhat.v=rep(NA,B) # save MLEs of a*b in the B bootstrap samples
  abhatH.v=matrix(NA,B)
  Index.m=matrix(NA,n,B)
  
  t1=0
  t2=0	
  for(i in 1:B){
    U=runif(n,min=1,max=n+1)
    index=floor(U)	
    Index.m[,i]=index
    #H(.05)
    Zrb=Zr[index,]
    SH=MeanCov(Zrb)$S
    thetaH=MLEst(SH)
    abhatH=thetaH[1]*thetaH[2]	
    abhatH.v[i]=abhatH
    if (abhatH<abH){
      t2=t2+1	
    }
    
  } # end of B loop
  
  abhatH.v=abhatH.v[!is.na(abhatH.v)]
  SEBH=sd(abhatH.v)
  
  # bootstrap confidence intervals using robust method
  CI2 =BpBCa(Zr,abhatH.v,t2,level)
  #    Results=list(CI=CI2)
  Results=list(CI=CI2[[1]],pv=CI2[[2]])
  return(Results)
  
}# end of function

BpBCa=function(Z,abhat.v,t,level){
  # Bootstrap percentile
  oab.v=sort(abhat.v)
  B=length(abhat.v)
  
  ranklowBp=round(B*level/2)
  
  if(ranklowBp==0){
    ranklowBp=1
  }
  
  Bpl=oab.v[ranklowBp]
  Bph=oab.v[round(B*(1-level/2))]	
  BP=c(Bpl,Bph)
  pstar=mean(oab.v>0)
  pv=2*min(c(pstar,1-pstar))
  #	Results=list(BP=BP)
  #    return(Results)
  list(BP,pv)
}

idealf<-function(x,na.rm=FALSE){
  #
  # Compute the ideal fourths for data in x
  #
  if(na.rm)x<-x[!is.na(x)]
  j<-floor(length(x)/4 + 5/12)
  y<-sort(x)
  g<-(length(x)/4)-j+(5/12)
  ql<-(1-g)*y[j]+g*y[j+1]
  k<-length(x)-j+1
  qu<-(1-g)*y[k]+g*y[k-1]
  list(ql=ql,qu=qu)
}

ifmest<-function(x,bend=1.28,op=2){
  #
  #   Estimate the influence function of an M-estimator, using
  #   Huber's Psi, evaluated at x.
  #
  #   Data are in the vector x, bend is the percentage bend
  #
  #  op=2, use adaptive kernel estimator
  #  otherwise use Rosenblatt's shifted histogram
  #
  tt<-mest(x,bend)  # Store M-estimate in tt
  s<-mad(x)*qnorm(.75)
  
#   if(op==2){
#     val<-akerd(x,pts=tt,plotit=FALSE,pyhat=T)
#     val1<-akerd(x,pts=tt-s,plotit=FALSE,pyhat=T)
#     val2<-akerd(x,pts=tt+s,plotit=FALSE,pyhat=T)
#   }
#   if(op!=2){
   val<-kerden(x,0,tt)
    val1<-kerden(x,0,tt-s)
    val2<-kerden(x,0,tt+s)
#  }
  ifmad<-sign(abs(x-tt)-s)-(val2-val1)*sign(x-tt)/val
  ifmad<-ifmad/(2*.6745*(val2+val1))
  y<-(x-tt)/mad(x)
  n<-length(x)
  b<-sum(y[abs(y)<=bend])/n
  a<-hpsi(y,bend)*mad(x)-ifmad*b
  ifmest<-a/(length(y[abs(y)<=bend])/n)
  ifmest
}

kerden<-function(x,q=.5,xval=0){
  #   Compute the kernel density estimator of the
  #   probability density function evaluated at the qth quantile.
  #
  #   x contains vector of observations
  #   q is the quantile of interest, the default is the median.
  #   If you want to evaluate f hat at xval rather than at the
  #   q th quantile, set q=0 and xval to desired value.
  #
  y<-sort(x)
  n<-length(x)
  temp<-idealf(x)
  h<-1.2*(temp$qu-temp$ql)/n^(.2)
  iq<-floor(q*n+.5)
  qhat<-y[iq]
  if (q==0) qhat<-xval
  xph<-qhat+h
  A<-length(y[y<=xph])
  xmh<-qhat-h
  B<-length(y[y<xmh])
  fhat<-(A-B)/(2*n*h)
  fhat
}

hpsi<-function(x,bend=1.28){
  #
  #   Evaluate Huber`s Psi function for each value in the vector x
  #   The bending constant defaults to 1.28.
  #
  hpsi<-ifelse(abs(x)<=bend,x,bend*sign(x))
  hpsi
}


yuenv2<-function(x,y=NULL,tr=.2,alpha=.05,plotit=FALSE,op=TRUE,VL=TRUE,cor.op=FALSE,loc.fun=median,
                 xlab="Groups",ylab="",PB=FALSE,nboot=100,SEED=FALSE){
  #plotfun=splot,
  #  Perform Yuen's test for trimmed means on the data in x and y.
  #  The default amount of trimming is 20%
  #  Missing values (values stored as NA) are automatically removed.
  #
  #  A confidence interval for the trimmed mean of x minus the
  #  the trimmed mean of y is computed and returned in yuen$ci.
  #  The significance level is returned in yuen$siglevel
  #
  #  For an omnibus test with more than two independent groups,
  #  use t1way.
  #
  #   Unlike the function yuen, a robust heteroscedastic measure
  #   of effect size is returned.
  #  PB=FALSE means that a Winsorized variation of prediction error is used to measure effect size.
  #  PB=TRUE:  A percentage bend measure of variation is used instead.
  #
  #if(tr==.5)stop("Use medpb to compare medians.")
  #if(tr>.5)stop("Can't have tr>.5")
  if(is.null(y)){
    if(is.matrix(x) || is.data.frame(x)){
      y=x[,2]
      x=x[,1]
    }
    if(is.list(x)){
      y=x[[2]]
      x=x[[1]]
    }
  }
  #library(MASS)
  #if(SEED)set.seed(2)
  x<-x[!is.na(x)]  # Remove any missing values in x
  y<-y[!is.na(y)]  # Remove any missing values in y
  n1=length(x)
  n2=length(y)
  h1<-length(x)-2*floor(tr*length(x))
  h2<-length(y)-2*floor(tr*length(y))
  q1<-(length(x)-1)*winvar(x,tr)/(h1*(h1-1))
  q2<-(length(y)-1)*winvar(y,tr)/(h2*(h2-1))
  df<-(q1+q2)^2/((q1^2/(h1-1))+(q2^2/(h2-1)))
  crit<-qt(1-alpha/2,df)
  m1=mean(x,tr)
  m2=mean(y,tr)
  mbar=(m1+m2)/2
  dif=m1-m2
  low<-dif-crit*sqrt(q1+q2)
  up<-dif+crit*sqrt(q1+q2)
  test<-abs(dif/sqrt(q1+q2))
  yuen<-2*(1-pt(test,df))
  xx=c(rep(1,length(x)),rep(2,length(y)))
  if(h1==h2){
    pts=c(x,y)
    top=var(c(m1,m2))
    #
    if(!PB){
      if(tr==0)cterm=1
      if(tr>0)cterm=area(dnormvar,qnorm(tr),qnorm(1-tr))+2*(qnorm(tr)^2)*tr
      bot=winvar(pts,tr=tr)/cterm
      e.pow=top/bot
      if(e.pow>1){
        x0=c(rep(1,length(x)),rep(2,length(y)))
        y0=c(x,y)
        e.pow=wincor(x0,y0,tr=tr)$cor^2
      }}
    #
#     if(PB){
#       bot=pbvar(pts)
#       e.pow=top/bot
#     }
    #
  }
  if(n1!=n2){
    N=min(c(n1,n2))
    vals=0
    for(i in 1:nboot)vals[i]=yuen.effect(sample(x,N),sample(y,N),tr=tr)$Var.Explained
    e.pow=loc.fun(vals)
  }
#   if(plotit){
#     plot(xx,pts,xlab=xlab,ylab=ylab)
#     if(op)
#       points(c(1,2),c(m1,m2))
#     if(VL)lines(c(1,2),c(m1,m2))
#   }
  list(ci=c(low,up),n1=n1,n2=n2,
       p.value=yuen,dif=dif,se=sqrt(q1+q2),teststat=test,
       crit=crit,df=df,Var.Explained=e.pow,Effect.Size=sqrt(e.pow))
}

yuen.effect<-function(x,y,tr=.2,alpha=.05,plotit=FALSE,
                      op=TRUE,VL=TRUE,cor.op=FALSE,
                      xlab="Groups",ylab="",PB=FALSE){
  #plotfun=splot,
  #  Same as yuen, only it computes explanatory power and the related
  # measure of effect size. Only use this with n1=n2. Called by yuenv2
  # which allows n1!=n2.
  #
  #
  #  Perform Yuen's test for trimmed means on the data in x and y.
  #  The default amount of trimming is 20%
  #  Missing values (values stored as NA) are automatically removed.
  #
  #  A confidence interval for the trimmed mean of x minus the
  #  the trimmed mean of y is computed and returned in yuen$ci.
  #  The significance level is returned in yuen$siglevel
  #
  #  For an omnibus test with more than two independent groups,
  #  use t1way.
  #  This function uses winvar from chapter 2.
  #
 # if(tr==.5)stop("Use medpb to compare medians.")
#  if(tr>.5)stop("Can't have tr>.5")
  x<-x[!is.na(x)]  # Remove any missing values in x
  y<-y[!is.na(y)]  # Remove any missing values in y
  h1<-length(x)-2*floor(tr*length(x))
  h2<-length(y)-2*floor(tr*length(y))
  q1<-(length(x)-1)*winvar(x,tr)/(h1*(h1-1))
  q2<-(length(y)-1)*winvar(y,tr)/(h2*(h2-1))
  df<-(q1+q2)^2/((q1^2/(h1-1))+(q2^2/(h2-1)))
  crit<-qt(1-alpha/2,df)
  m1=mean(x,tr)
  m2=mean(y,tr)
  mbar=(m1+m2)/2
  dif=m1-m2
  low<-dif-crit*sqrt(q1+q2)
  up<-dif+crit*sqrt(q1+q2)
  test<-abs(dif/sqrt(q1+q2))
  yuen<-2*(1-pt(test,df))
  xx=c(rep(1,length(x)),rep(2,length(y)))
  pts=c(x,y)
  top=var(c(m1,m2))
  #
  if(!PB){
    if(tr==0)cterm=1
    if(tr>0)cterm=area(dnormvar,qnorm(tr),qnorm(1-tr))+2*(qnorm(tr)^2)*tr
    bot=winvar(pts,tr=tr)/cterm
  }
  #if(PB)bot=pbvar(pts)/1.06
  #
  e.pow=top/bot
  if(e.pow>1){
    x0=c(rep(1,length(x)),rep(2,length(y)))
    y0=c(x,y)
    e.pow=wincor(x0,y0,tr=tr)$cor^2
  }
 
  list(ci=c(low,up),p.value=yuen,dif=dif,se=sqrt(q1+q2),teststat=test,
       crit=crit,df=df,Var.Explained=e.pow,Effect.Size=sqrt(e.pow))
}


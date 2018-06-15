t1way <- function(formula, data, tr = 0.2) {

  if (missing(data)) {
    mf <- model.frame(formula)
  } else {
    mf <- model.frame(formula, data)
  }
  cl <- match.call()
  
  x <- split(model.extract(mf, "response"), mf[,2])   
  
  if (tr==0.5) warning("Comparing medians should not be done with this function!")

  grp <- 1:length(x) 
  
  J <-length(grp)
  h <- vector("numeric",J)
  w <- vector("numeric",J)
  xbar <- vector("numeric",J)
  nv <- NA
  
  pts=NULL
  for(j in 1:J){
    xx <- !is.na(x[[j]])
    val <- x[[j]]
    x[[j]] <- val[xx]  # Remove missing values
    nv[j] <- length(x[[j]])
    h[j] <- length(x[[grp[j]]])-2*floor(tr*length(x[[grp[j]]]))
    w[j] <- h[j]*(h[j]-1)/((length(x[[grp[j]]])-1)*winvar(x[[grp[j]]],tr))
    if (winvar(x[[grp[j]]],tr) == 0) stop("Standard error cannot be computed because of Winsorized variance of 0 (e.g., due to ties). Try do decrease the trimming level.") 
    xbar[j] <- mean(x[[grp[j]]],tr)
    val<-elimna(val)
    pts=c(pts,val)
 }
 u <- sum(w)
 xtil <- sum(w*xbar)/u
 A <- sum(w*(xbar-xtil)^2)/(J-1)
 B <- 2*(J-2)*sum((1-w/u)^2/(h-1))/(J^2-1)
 TEST <- A/(B+1)
 nu1 <- J-1
 nu2 <- 1./(3*sum((1-w/u)^2/(h-1))/(J^2-1))
 sig <- 1-pf(TEST,nu1,nu2)
 
 
 
 chkn=var(nv)
 if(chkn==0){
   top=var(xbar)
   bot=winvarN(pts,tr=tr)
   e.pow=top/bot
 }
 if(chkn!=0){
   vals=0
   N=min(nv)
   xdat=list()
   nboot = 100
   for(i in 1:nboot){
     for(j in 1:J){
       xdat[[j]]=sample(x[[j]],N)
     }
     vals[i]=t1way.effect(xdat,tr=tr)$Var.Explained
   }
   loc.fun=median
   e.pow=loc.fun(vals,na.rm=TRUE)
 }
 
 result <- list(test = TEST, df1 =nu1, df2 = nu2, p.value = sig, effsize = sqrt(e.pow), call = cl)
 class(result) <- c("t1way")
 result
}

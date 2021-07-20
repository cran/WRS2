t1way <- function(formula, data, tr = 0.2, alpha = 0.05, nboot = 100, ...) {

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

 ## effect size
 chkn=var(nv)
 if(chkn==0){                       ## equal sample size
   top=var(xbar)
   bot=winvarN(pts,tr=tr)
   e.pow = sqrt(top/bot)            ## exact computation
 }

## unequal sample size and bootstrap CI
 vals=0
 N=min(nv)
 xdat=list()
 for(i in 1:nboot){
   for(j in 1:J){
     xdat[[j]] <- sample(x[[j]], N, replace = TRUE)
   }
   #vals[i]=t1way.effect(xdat,tr=tr)$Var.Explained
   vals[i] <- t1wayv2(xdat, tr = tr, nboot = 5, SEED = FALSE)$Effect.Size
 }
 loc.fun <- median
 if(chkn!=0) e.pow <- loc.fun(vals,na.rm=TRUE)    ## unequal sample size effect size (bootstrap)

## CI computation
 ilow <- round((alpha/2) * nboot)
 ihi <- nboot - ilow
 ilow <- ilow+1
 val <- sort(vals)
 ci <- val[ilow]
 ci[2] <- val[ihi]
   #if(chk$p.value>alpha)ci[1]=0


 result <- list(test = TEST, df1 =nu1, df2 = nu2, p.value = sig, effsize = e.pow, effsize_ci = ci, alpha = alpha,
                call = cl)
 class(result) <- c("t1way")
 result
}

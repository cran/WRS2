print.bwtrim <-
function(x,...)
{
  cat("Call:\n")
  print(x$call)
  cat("\n")
  
 dfx <- data.frame(value = c(x$Qa, x$Qb, x$Qab), df1 = c(x$A.df[1], x$B.df[1], x$AB.df[1]), df2 = c(x$A.df[2], x$B.df[2], x$AB.df[2]), 
                  p.value = c(x$A.p.value, x$B.p.value, x$AB.p.value))
 rownames(dfx) <- c(x$varnames[2], x$varnames[3], paste0(x$varnames[2], ":", x$varnames[3]))
 dfx <- round(dfx, 4)
 print(dfx)
 cat("\n")
}



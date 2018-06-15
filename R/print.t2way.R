print.t2way <-
function(x,...)
{
  cat("Call:\n")
  print(x$call)
  cat("\n")
  
  if (length(grep("med2way", x$call)) == 0) {
    if (!is.na(x$Qa)) {
      df <- data.frame(value = c(x$Qa, x$Qb, x$Qab), p.value = c(x$A.p.value, x$B.p.value, x$AB.p.value))
    } else {
      df <- data.frame(p.value = c(x$A.p.value, x$B.p.value, x$AB.p.value))
    }
    rownames(df) <- c(x$varnames[2], x$varnames[3], paste0(x$varnames[2], ":", x$varnames[3]))
    print(round(df,4))
    cat("\n")
  } else {
    df <- data.frame(value = c(x$Qa, x$Qb, x$Qab), p.value = c(x$A.p.value, x$B.p.value, x$AB.p.value))
    rownames(df) <- c(x$varnames[2], x$varnames[3], paste0(x$varnames[2], ":", x$varnames[3]))
    df <- round(df, 4)
    free <- c(paste0("F(", x$dim[1]-1, ",", Inf, ")"), paste0("F(", x$dim[2]-1, ",", Inf, ")"), 
              paste0("Chisq(", prod(x$dim-1), ")"))
    df1 <- data.frame(value = df[,1], df = free, p.value = df[,2])
    print(df1)
    cat("\n")
  }
}

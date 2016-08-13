print.pbcor <- function(x,...)
  {
    cat("Call:\n")
    print(x$call)
    
    cat("\nRobust correlation coefficient:", round(x$cor, 4))
    cat("\nTest statistic:", round(x$test, 4)) 
    cat("\np-value:", round(x$p.value, 5), "\n\n")
  }

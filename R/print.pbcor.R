print.pbcor <- function(x,...)
  {
    cat("Call:\n")
    print(x$call)
    
    cat("\nRobust correlation coefficient:", round(x$cor, 4))
    cat("\nTest statistic:", round(x$test, 4)) 
    cat("\np-value:", round(x$p.value, 5), "\n")
    
    if (!(is.na(x$cor_ci[1]))) cat("\nBootstrap CI: [", round(x$cor_ci[1],4), "; ", round(x$cor_ci[2],4), "]\n\n", sep = "")
  }

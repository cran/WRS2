print.twocor <- function(x,...)
  {
    cat("Call:\n")
    print(x$call)
    
    cat("\nFirst correlation coefficient:", round(x$r1, 4))
    cat("\nSecond correlation coefficient:", round(x$r2, 4))
    cat("\nConfidence interval (difference):", round(x$ci, 4)) 
    if (!is.null(x$p.value)) cat("\np-value:", round(x$p.value, 5), "\n\n") else cat("\n\n")
  }

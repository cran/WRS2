print.pball <- function(x,...)
  {
    cat("Call:\n")
    print(x$call)
    
    cat("\nRobust correlation matrix:\n")
    if (!is.null(x$pbcorm)) print(round(x$pbcorm, 4)) else print(round(x$cor, 4))
    cat("\np-values:\n")
    print(round(x$p.values, 5))
    
    if (!is.null(x$H)) cat("\n\nTest statistic H: ", round(x$H, 4), ", p-value = ", round(x$H.p.value, 5), "\n\n", sep = "") 
  }

print.robmed <- function(x,...)
{
  cat("Call:\n")
  print(x$call)
  
  cat("\nMediated effect:", round(x$ab.est, 4), "\n")
  cat("Confidence interval:", round(x$CI.ab, 4), "\n")
  cat("p-value:", round(x$p.value, 4), "\n")
  cat("\n")
}

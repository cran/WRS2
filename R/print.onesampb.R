print.onesampb <- function(x,...) {
  cat("Call:\n")
  print(x$call)
  cat("\n")
  cat("Robust location estimate:", round(x$estimate, 4), "\n")
  cat("95% confidence interval:", round(x$ci, 4), "\n")
  cat("p-value:", round(x$p.value, 4), "\n\n")
}
  
    
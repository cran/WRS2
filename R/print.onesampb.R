print.onesampb <- function(x,...) {
  cat("Call:\n")
  print(x$call)
  cat("\n")
  cat("Robust location estimate:", round(x$estimate, 4), "\n")
  cistr <- paste0(1-x$alpha, "% confidence interval:")
  cat(cistr, round(x$ci, 4), "\n")
  cat("p-value:", round(x$p.value, 4), "\n\n")
}
  
    
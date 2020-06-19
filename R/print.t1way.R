print.t1way <-
function(x,...)
{
  cat("Call:\n")
  print(x$call)
  
  cat("\nTest statistic: F =", round(x$test, 4),"\n")
  cat("Degrees of freedom 1:", round(x$df1, 2),"\n")
  cat("Degrees of freedom 2:", round(x$df2, 2) ,"\n")
  cat("p-value:", round(x$p.value, 5), "\n")
  cat("\n")
  if (!is.null(x$effsize)) {
    cat("Explanatory measure of effect size:", round(x$effsize, 2), "\n")
    cat("Bootstrap CI: [", round(x$effsize_ci[1], 2), "; ", round(x$effsize_ci[2], 2), "]\n\n", sep = "")
  }
}

print.med1way <-
function(x,...)
{
  cat("Call:\n")
  print(x$call)
  cat("\nTest statistic F:", round(x$test, 4),"\n")
# cat("Degrees of freedom 1:", x$df1, "\n")
# cat("Degrees of freedom 2:", x$df2, "\n")
  if (!is.na(x$crit.val)) cat("Critical value:", round(x$crit.val, 4),"\n")
  cat("p-value:", round(x$p.value, 5), "\n")
  cat("\n")
}

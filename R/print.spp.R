print.spp <-
function(x,...)
{
  cat("Call:\n")
  print(x$call)
  
  cat("\nTest statistics:\n")
  #cat(round(x$test, 4), "\n")
  print(x$test, digits = 4)
  cat("\nTest whether the corrresponding population parameters are the same:\n")
  cat("p-value:", round(x$p.value, 5), "\n")
  cat("\n")
}

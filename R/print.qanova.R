print.qanova <-
function(x,...)
{
  cat("Call:\n")
  print(x$call)
  
  #cat("\nTest statistics: ", round(x$psihat, 4))
  #cat("\np-value: ", round(x$p.value, 5), "\n\n")
  
  #partable <- data.frame(x$psihat, x$p.value)
  partable <- data.frame(x$p.value)
  cat("\n")
  print(round(partable, 4))
  cat("\n")
}

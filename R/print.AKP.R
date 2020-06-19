print.AKP <-
function(x,...)
{
  cat("Call:\n")
  print(x$call)
  
  cat("\nAKP effect size:", round(x$AKPeffect, 2), "\n")
  cat("Bootstrap CI: [", round(x$AKPci[1], 2), "; ", round(x$AKPci[2], 2), "]\n\n", sep = "")
}

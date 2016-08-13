print.robtab <-
function(x,...)
{
  cat("Call:\n")
  print(x$call)
  
  cat("\nParameter table: \n")
  if (ncol(x$partable) >= 8) {
    print(round(x$partable, 4))
  } else {  
    print(x$partable)
  }
  cat("\n")
}

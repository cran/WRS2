trimcibt <- function(x,
                     nv = 0,
                     tr = 0.2,
                     alpha = 0.05,
                     nboot = 200,
                     ...) {
  test <- (mean(x, tr) - nv) / WRS2::trimse(x, tr)
  data <- matrix(sample(x, size = length(x) * nboot, replace = TRUE), nrow = nboot) - mean(x, tr)
  top <- apply(data, 1, mean, tr)
  bot <- apply(data, 1, WRS2::trimse, tr)
  tval <- sort(abs(top / bot))
  icrit <- round((1 - alpha) * nboot)
  ibot <- round(alpha * nboot / 2)
  itop <- nboot - ibot
  trimcibt <- mean(x, tr) - tval[icrit] * WRS2::trimse(x, tr)
  trimcibt[2] <- mean(x, tr) + tval[icrit] * WRS2::trimse(x, tr)
  p.value <- (sum(abs(test) <= abs(tval))) / nboot

  result <- list(
    estimate = mean(x, tr),
    ci = trimcibt,
    test.stat = test,
    tr = tr,
    p.value = p.value,
    n = length(x)
  )

  class(result) <- "trimcibt"
  result
}

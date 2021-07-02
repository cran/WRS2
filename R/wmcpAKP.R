wmcpAKP <- function(x, tr = 0.2, nboot = 200, ...) {
  #  Compute Algina et al measure of effect size for pairs of J dependent groups
  if (is.matrix(x) || is.data.frame(x)) x <- listm(x)
  J <- length(x)
  C <- (J^2 - J) / 2
  A <- matrix(NA, nrow = C, ncol = 5)
  dimnames(A) <- list(NULL, c("Group", "Group", "Effect.Size", "ci.low", "ci.up"))
  ic <- 0
  for (j in 1:J) {
    for (k in 1:J) {
      if (j < k) {
        ic <- ic + 1
        A[ic, 1] <- j
        A[ic, 2] <- k
        A[ic, 3] <- dep.effect(x[[j]], x[[k]], tr = tr, nboot = nboot)[5]
        A[ic, 4] <- dep.effect(x[[j]], x[[k]], tr = tr, nboot = nboot)[21]
        A[ic, 5] <- dep.effect(x[[j]], x[[k]], tr = tr, nboot = nboot)[25]
      }
    }
  }

  res <- c("Effect.Size" = mean(A[, 3]), "ci.low" = mean(A[, 4]), "ci.up" = mean(A[, 5]))

  class(res) = c("numeric", "wmcpAKP")
  res
}





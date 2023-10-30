rmultinomial <- function (n, size, prob) {
    if (length(n) > 1)
        n <- length(n)
    if (length(n) == 0 || as.integer(n) == 0)
        return(numeric(0))
    n <- as.integer(n)
    if (n < 0)
        stop("integer(n) can not be negative in rmultinomial")
    if (is.vector(prob) || (dim(prob)[1]) == 1) {
        if (length(size) == 1)
            return(t(rmultinom(n, size, prob)))
        prob <- matrix(prob, nrow = 1)
    }
    nrp <- nrow(prob)
    mnr <- min(max(nrp, length(size)), n)
    ss <- rep(size, length.out = mnr)
    if (nrp != mnr)
        prob <- matrix(t(prob), ncol = ncol(prob), nrow = mnr, byrow = TRUE)
    n1 <- n%/%mnr
    n2 <- n%%mnr
    res <- sapply(1:mnr, function(x) rmultinom(n1, ss[x], prob[x,]))
    res <- matrix(res, ncol = ncol(prob), byrow = TRUE)
    index <- as.vector(matrix(1:(mnr * n1), ncol = mnr, byrow = TRUE))
    res <- res[index, ]
    if (n2 != 0) {
        res <- rbind(res, t(sapply(1:n2, function(x) rmultinom(1, ss[x], prob[x, ]))))
    }
    return(res)
}

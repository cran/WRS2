Qanova <- function(formula, data, q = 0.5, nboot = 600){
  #
  # Test global hypothesis that J independent groups
  # have equal medians.
  # Performs well when there are tied values.
  #
  # Basically, use pbadepth in conjunction with the Harrell--Davis
  # estimator.
  #
  #
  if (missing(data)) {
    mf <- model.frame(formula)
  } else {
    mf <- model.frame(formula, data)
  }
  cl <- match.call()
  
  x <- split(model.extract(mf, "response"), mf[,2])   
  
  op=3
  MC <- FALSE
  
  chkcar=NA
  for(j in 1:length(x))chkcar[j]=length(unique(x[[j]]))
  #if(min(chkcar<20)) warning("Cardinality of sample space is less than 20 for one more groups. Type I error might not be controlled!")

  nc <- length(pbadepth1(x, est = hd, q = q[1], allp = TRUE, SEED = FALSE, op = op, nboot = nboot, MC = MC, na.rm = TRUE)$psihat)
  psimat <- matrix(NA, length(q), nc)
  pvals <- NA
  for (i in 1:length(q)) {
     output <- pbadepth1(x, est = hd, q = q[i], allp = TRUE, SEED = FALSE, op = op, nboot = nboot, MC = MC, na.rm = TRUE)
     psimat[i,] <- output$psihat
     pvals[i] <- output$p.value
  }
  
  psidf <- as.data.frame(psimat)
  rownames(psidf) <- paste0("q = ",q)
  colnames(psidf) <- paste0("con", 1:nc)
  
  pvalues <- as.data.frame(cbind(pvals, p.adjust(pvals, method = 'hochberg')))
  rownames(pvalues) <- paste0("q = ",q)
  colnames(pvalues) <- c("p-value", "p-adj")
  
  result <- list(psihat = psidf, p.value = pvalues, contrasts = output$con, call = cl)
  class(result) <- "qanova"
  result
}
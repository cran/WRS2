disc2.chi.sq<-function(x,y,simulate.p.value=FALSE,B=2000){
  #
  # Test the hypothesis of identical discrete distributions
  # using a chi-squared test and a simulated p.value
  #
  n1 = length(x)
  n2 = length(y)
  g = c(rep(1,n1), rep(2,n2))
  d = c(x,y)
  df = data.frame(d, g)
  res=chisq.test(df$d, df$g, simulate.p.value=simulate.p.value, B=B)
  list(X.squared=res[1]$statistic,p.value=res[3]$p.value)
}

disc2com <- disc2.chi.sq

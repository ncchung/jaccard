getp <- function(lr,lr0,ties.method = "random") {
  # Get resampled p-values, pulling across variables (e.g., genes)
  # lr: observed statistics
  # lr0: null statistics (i.e. from resampled residuals)

  m = length(lr)
  v = c(rep(TRUE,m),rep(FALSE,length(lr0)))
  v = v[rev(order(c(lr,lr0)))]
  u = 1:length(v)
  w = 1:m
  p = ((u[v==TRUE]-w)+1)/(length(lr0)+2)
  #to account for the extreme observed values
  #p = ((u[v==TRUE]-w)+1)/(length(lr0)+2)
  p = p[rank(-lr, ties.method = ties.method)]
  return(p)
}

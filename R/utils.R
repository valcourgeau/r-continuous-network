RowNormalised <- function(adj){
  # zero on diag, row sum to 1 with positive weights
  
  diag(adj) <- 0.0
  adj[adj!=0.0] <- 1.0
  divisors <- pmax(rowSums(adj), 1.0)
  
  return(diag(1.0/divisors) %*% adj)
}

ConcatenateCol <- function(col_vec, n){
  return(matrix(rep(col_vec, n), ncol = n))
}

CleanData <- function(data, frequency=24, s.window=24, t.window=24, ...){
  stl_cleaned <- lapply(1:ncol(data),
                        FUN = function(i){
                          stl(ts(data[,i], frequency = frequency), s.window = s.window, t.window = t.window, ...)
                        }
  )
  names(stl_cleaned) <- colnames(data)
  remainders <- vapply(stl_cleaned, function(x){as.numeric(x$time.series[,3])}, rep(0, length(data[,1])))
  colnames(remainders) <- colnames(data)
  std.dev <- apply(remainders, MARGIN = 2, sd)
  # remainders <- t(apply(remainders, 1, function(x){x/std.dev}))
  return(list(stl_obj=stl_cleaned, remainders=remainders, std.dev=std.dev))
}
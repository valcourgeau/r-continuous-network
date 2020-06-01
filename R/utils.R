RowNormalised <- function(adj){
  # zero on diag, row sum to 1 with positive weights
  
  diag(adj) <- 0.0
  adj[adj!=0.0] <- 1.0
  divisors <- pmax(rowSums(adj), 1.0)
  
  return(diag(1.0/divisors) %*% adj)
}


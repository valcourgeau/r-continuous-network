setwd('~/GitHub/r-continuous-network/R/')
source('package_to_load.R')
source('adjacency_generation.R')

#' @example ConstructPath(diag(10), matrix(rnorm(100*10,0,1), ncol=10), 10, 0.01)
ConstructPath <- function(nw_topo, noise, y_init, delta_time){
  # generates path given noise increments and fixed delta_time.
  assertthat::assert_that(delta_time > 0)
  assertthat::are_equal(ncol(nw_topo), nrow(nw_topo))
  
  n <- nrow(noise)
  d <- ncol(nw_topo)
  d_y <- matrix(0, nrow=n, ncol=d)
  d_y[1,] <- y_init
  for(i in 2:n){
    d_y[i,] <-  d_y[i-1,] - as.vector(nw_topo %*% d_y[i-1,]) * delta_time + noise[i,]
  }
  
  return(d_y)  
}


CorrelatedBrownianNoise <- function(sigma_matrix, n_rows, delta_time){
  n_nodes <- ncol(sigma_matrix)
  assertthat::assert_that(n_nodes > 1)
  sigma_matrix <- stats::cov2cor(sigma_matrix) # standard
  uncorrelated_gaussian <- matrix(rnorm(n_rows*n_nodes, mean = 0, sd = sqrt(delta_time)), nrow=n_nodes)
  print('ok')
  sigma_matrix_copy <- sigma_matrix
  diag(sigma_matrix_copy) <- 0.0
  diag(sigma_matrix) <- sqrt(1-rowSums(sigma_matrix_copy^2))
  sigma_matrix[upper.tri(sigma_matrix, diag = FALSE)] <- 0.0
  corr_bm <- t(sigma_matrix %*% uncorrelated_gaussian)
  print(corr_bm %>% dim)
  integrated_bm <- apply(corr_bm, 2, cumsum)
  return(integrated_bm)
}

set.seed(42)
corr_bw <- CorrelatedBrownianNoise(sigma_matrix = matrix(c(1, 0, 0, 1), 2, 2), n_rows = 5e4, delta_time = 0.01)
plot(corr_bw[,1], type='l')
lines(corr_bw[,2])


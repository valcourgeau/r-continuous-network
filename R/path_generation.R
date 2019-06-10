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

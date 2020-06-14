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

CorrelatedBrownianNoise <- function(sigma_matrix, n, delta_time){
  n_nodes <- ncol(sigma_matrix)
  assertthat::assert_that(n_nodes > 1)
  sigs <- diag(sigma_matrix)
  sigma_matrix <- stats::cov2cor(sigma_matrix) # standard
  uncorrelated_gaussian <- matrix(rnorm(n*n_nodes, mean = 0, sd = sqrt(delta_time)), nrow=n_nodes)
  sqrt_sigma <- Matrix::chol(sigma_matrix) #pracma::sqrtm(sigma_matrix)$B# chol(sigma_matrix)#
  corr_bm <- t(sqrt_sigma %*% uncorrelated_gaussian * as.vector(sigs))
  return(corr_bm)
}

# Covariance/Correlation test
set.seed(42)
corr_bw <- lapply(
  1:10000,
  function(i){CorrelatedBrownianNoise(sigma_matrix = matrix(c(1.0, 0.2, 0.2, 0.2, 1.0, 0.2, 0.2, 0.2, 1.0), 3, 3), n = 1e3, delta_time = 0.01)}
)
corr_bw <- do.call(rbind, lapply(corr_bw, function(x){apply(x, 2, sum)})) # comparing only the last element across samples
cor(corr_bw)

CompoundPoissonJumps <- function(d, n, delta_time, jump_values){
  if(d==1){
    n_jumps <- length(jump_values)
  }else{
    n_jumps <- nrow(jump_values)
  }
  if(n < n_jumps) warning('n < n_jumps in CompoundPoissonJumps')
  results <- matrix(0, nrow = n, ncol = d)
  assertthat::are_equal(d, ncol(jump_values))
  
  horizon <- n * delta_time
  jump_times <- lapply(1:d, function(x){sort(runif(min = 0, max = horizon, n = n_jumps))})
  time_grid <- seq(0, to = delta_time*n, by=delta_time)
  idx_jump_injection <- lapply(
    jump_times,
    function(times){
      vapply(times, 
             FUN = function(jump_time){
                   pmin(which.max(jump_time < time_grid), n)
               }, 3)
    }
  )
  
  for(i in 1:length(idx_jump_injection)){
    jump_value_idx <- 1
    for(idx in idx_jump_injection[[i]]){
      results[idx, i] <- results[idx, i] + jump_values[jump_value_idx, i]
      jump_value_idx <- jump_value_idx + 1
    }
  }
  return(list(noise=results, jump_times=do.call(cbind, jump_times)))
}

# n_jumps <- 1000
# cmpnd_poisson <- CompoundPoissonJumps(
#   d = 3, n = 1000, delta_time = 0.001,
#   jump_values=matrix(rep(rnorm(n_jumps, 0, 0.1), n_jumps*3), nrow = n_jumps, ncol = 3))
# cmpnd_poisson_noise <-  apply(cmpnd_poisson$noise, 2, cumsum)
# sig_matrix <- matrix(-.2, 3, 3)
# diag(sig_matrix) <- 1
# corr_bw <- CorrelatedBrownianNoise(sigma_matrix = sig_matrix, n = 1000, delta_time = 0.001)
# corr_bw <-  apply(corr_bw, 2, cumsum)
# noise_add <- cmpnd_poisson_noise + corr_bw
# 
# plot(noise_add[,1], col='darkred', type='l', ylim=c(-3,3))
# lines(corr_bw[,1], col='red')
# abline(v=cmpnd_poisson$jump_times[,1]/0.001, col='darkred')
# lines(noise_add[,2], col='darkblue', type='l')
# lines(corr_bw[,2], col='blue')
# abline(v=cmpnd_poisson$jump_times[,2]/0.001, col='darkblue')


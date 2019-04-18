source("package_to_load.R")
LOAD_EXAMPLES <- FALSE


### Functions to compute cts-time NAR-equivalent MLE and drift MLE

EnsureTopo <- function(nw_topo){
  # creates matrix with entries either 0 or 1.
  nw_topo[nw_topo != 0] <- 1
  return(nw_topo)
}

# Example
nw_topo <- c(2,1,1,
             0,2,3,
             0,0,-1)
nw_topo <- matrix(nw_topo, ncol = 3, nrow = 3)
nw_topo
EnsureTopo(nw_topo = nw_topo)

test_topo <- matrix(1, ncol=3, nrow=3)
EnsureTopo(nw_topo = test_topo)

StdTopo <- function(nw_topo){
  nw_topo <- EnsureTopo(nw_topo = nw_topo)
  diag(nw_topo) <- 0 
  std_coeff <- rowSums(nw_topo)
  std_coeff[std_coeff == 0] <- 1

  return(nw_topo/std_coeff)
}

# Example
nw_topo <- c(2,1,1,
             0,2,3,
             0,0,1)
nw_topo <- matrix(nw_topo, ncol = 3, nrow = 3)
StdTopo(nw_topo = nw_topo)

test_topo <- matrix(1, ncol=3, nrow=3)
StdTopo(nw_topo = test_topo)

DataFiltering <- function(data, thresholds, diff_values=F, one_d=F){
  # Shall we have the first or second point when diff_values=F?
  filtered_data <- diff(data)
 
  if(one_d){
    N <- length(data)
    indicator_data <- (abs(filtered_data) <= thresholds)
    if(diff_values){
      return(filtered_data * indicator_data)
    }else{
      return(data[-N,] * indicator_data)
    }
  }else{
    N <- length(data[,1])
    under_thres <- t(apply(abs(filtered_data), 
                           MARGIN = 1, '<=', thresholds))
    if(diff_values){
      return(filtered_data*under_thres)
    }else{
      return(data[-N,] * under_thres)
    }
  }
}

# Example
data <- c(1,2,5,
          1,3,6,
          1,4,7)
data <- t(matrix(data, ncol=3))
data
DataFiltering(data, c(2,3,6))
DataFiltering(data, c(2,3,6), diff = T)

TimeMatrix <- function(times, ncol, one_d=F){
  # remove one_d just usel ncol
  times <- diff(times)
  
  if(one_d){
    return(times)
  }else{
    return(t(matrix(rep(times, each=ncol), nrow=ncol)))
  }
}

# Example
times <- c(1,2,4,7,11)
TimeMatrix(times, ncol = 2)
TimeMatrix(times, ncol = 2)

NOUfit1D <- function(times, data, threshold){
    N <- length(data)
    if(threshold <= 0){
      stop('threshold should be positive.')
    }
    
    diff_filtered <- DataFiltering(data,
                                   thresholds = threshold,
                                   diff_values=T, one_d = T)
    diff_times <- TimeMatrix(times = times, ncol=1, one_d=T)

    mle_estimate_up <- data[-N] * diff_filtered 
    mle_estimate_down <- data[-N]^2 * diff_times
    return(-sum(mle_estimate_up) / sum(mle_estimate_down))
}

# Example
if(LOAD_EXAMPLES)
{
  source("network_generation.R")
  d <- 1 #dims
  N <- 24*1500 #numbers of points
  M <- 10 # number of simulations
  Y0 <- 1 #start point
  delta_t <- 1/24
  draw_n <- seq(from=8, to=N, by=1000)
  mle_fit_example <- matrix(0, ncol=length(draw_n),
                            nrow=M)
  
  
  sigma <- 0.1
  set.seed(42)
  nw_topo <- genRdmAssymetricGraphs(d = d, p.link = 0.25,
                                    theta_1 = 1, theta_2 = 0)
  nw_topo <- StdTopo(nw_topo)
  
  times <- seq(from = 0, by = delta_t, length.out = N)
  for(i in 1:M){
    nw_data_bm <- matrix(rnorm(n = d*N, mean = 0, sd = sigma*sqrt(delta_t)), ncol = d)
    nw_data <- matrix(0, ncol=d, nrow=N)
    nw_data[1,] <- nw_data_bm[1,]*sqrt(delta_t)
    
    nw_q <- 0.1*nw_topo
    diag(nw_q) <- 1
    
    for(index in 2:N){
      nw_data[index,] <-  nw_data[index-1,] - (nw_q %*% nw_data[index-1,]) * (times[index]-times[index-1]) + nw_data_bm[index,]
    }
    
    #plot(nw_data[1:300,1], type="l")
    
    j <- 0
    for(n_temp in draw_n){
      mle_fit_example[i,j] <- NOUfit1D(times = times[1:n_temp], data = nw_data[1:n_temp], threshold = 5)
      j <- j + 1
    }
  }
  
  plot(draw_n[-length(draw_n)], colMeans(mle_fit_example[,-length(draw_n)]), 
       xlab="Sample size", ylab="1D MLE", main="1D MLE as sample size",
       ylim=c(0,2), type="l", col="red")  
  lines(draw_n[-length(draw_n)], apply(mle_fit_example[,-length(draw_n)],
                                      MARGIN = 2, 
                                      FUN = function(x){quantile(x,0.95)}), 
       xlab="Sample size", ylab="1D MLE", main="1D MLE as sample size",
       ylim=c(0,2), type="l", lty=2, col="red") 
  lines(draw_n[-length(draw_n)], apply(mle_fit_example[,-length(draw_n)],
                                       MARGIN = 2, 
                                       FUN = function(x){quantile(x,0.05)}), 
        xlab="Sample size", ylab="1D MLE", main="1D MLE as sample size",
        ylim=c(0,2), type="l", lty=2, col="red") 
  
  for(i in 1:M){
    lines(draw_n[-length(draw_n)], mle_fit_example[i,-length(draw_n)], 
          col=rgb(red = 0, green = 0, blue = 0, alpha = 0.02))  
  }
}


NOUfit <- function(nw_topo, times, data, thresholds){
  # fitting MLE to the network
  # TODO add support for jumps by using s- instead of s
  
  # Check 
  nw_topo <- StdTopo(nw_topo = nw_topo)
  nw_test <- nw_topo
  diag(nw_test) <- 1 
  if(max(Re(eigen(nw_test)$values)) <= 0){
    stop("Q matrix should be positive definite.")
  }
  
  remove(nw_test) # delete test matrix
  
  d <- length(data[1,])
  N <- length(data[,1])

  if(length(thresholds) != d){
    stop('Wrong dimension between data and thresholds.')
  }
  if(length(times) != N){
    stop('Wrong dimensions between data and times.')
  }
  
  diff_filtered <- DataFiltering(data, thresholds, diff_values=T)
  diff_times <- TimeMatrix(times = times, ncol=d)
  a_bar_y_t <- t(nw_topo %*% t(data))
  a_bar_y_t <- a_bar_y_t[-N,]
  y_t <- data[-N,]

  wide_mle_xi <- matrix(c(
      sum(a_bar_y_t^2*diff_times), sum(a_bar_y_t*y_t*diff_times),
      sum(a_bar_y_t*y_t*diff_times), sum(y_t^2*diff_times)
    ), ncol=2)
  wide_mle_c_filtered <- c(
      sum(a_bar_y_t * diff_filtered),
      sum(y_t * diff_filtered)
    ) 
  
  #return(solve(a = wide_mle_xi, b = wide_mle_c_filtered))
  results <- list("MLE_wide"=
                    -as.vector(solve(a = wide_mle_xi, b = wide_mle_c_filtered)))
  
  nw_topo <- results[["MLE_wide"]][1] * nw_topo
  diag(nw_topo) <- results[["MLE_wide"]][2]
  results[["Q"]] <- nw_topo
  
  return(results)
} 

# Example
#source("network_generation.R")
if(LOAD_EXAMPLES)
{
  source("network_generation.R")
  d <- 5 #dims
  N <- 24*1500 #numbers of points
  M <- 10 # number of simulations
  Y0 <- 1 #start point
  delta_t <- 1/24
  draw_n <- seq(from=8, to=N, by=1000)
  mle_fit_example <- rep(0, length(draw_n)*M*2)
  mle_fit_example <- array(data = mle_fit_example,
                           dim=c(M,length(draw_n),2))
  
  sigma <- 0.1
  set.seed(42)
  nw_topo <- genRdmAssymetricGraphs(d = d, p.link = 0.25,
                                    theta_1 = 1, theta_2 = 0)
  nw_topo <- StdTopo(nw_topo)
  
  times <- seq(from = 0, by = delta_t, length.out = N)
  thresholds <- rep(0, d)
  for(i in 1:M){
    nw_data_bm <- matrix(rnorm(n = d*N, mean = 0, sd = sigma*sqrt(delta_t)), ncol = d)
    nw_data <- matrix(0, ncol=d, nrow=N)
    nw_data[1,] <- nw_data_bm[1,]*sqrt(delta_t)
    
    nw_q <- 0.1*nw_topo
    diag(nw_q) <- 1
    
    for(index in 2:N){
      nw_data[index,] <-  nw_data[index-1,] - (nw_q %*% nw_data[index-1,]) * (times[index]-times[index-1]) + nw_data_bm[index,]
    }
    
    #plot(nw_data[1:300,1], type="l")
    
    j <- 0
    for(n_temp in draw_n){
      res <- NOUfit(nw_topo = nw_topo, times = times[1:n_temp],
                    data = nw_data[1:n_temp,], thresholds = thresholds)
      mle_fit_example[i,j,] <- res$MLE_wide

      j <- j + 1
    }
  }
  mle_true <- c(0.1, 1.0)
  par(mfrow=(c(1,2)))
  index_mle_param <- 1
  
  plot(draw_n[-length(draw_n)], colMeans(mle_fit_example[,-length(draw_n),index_mle_param])
       - mle_true[index_mle_param], main="First para MLE / Error and Empirical var")
  abline(h=0)
  lines(draw_n[-length(draw_n)], colMeans((mle_fit_example[,-length(draw_n),index_mle_param]
         -colMeans(mle_fit_example[,-length(draw_n),index_mle_param]))^2),
        col="red", lty=2, lwd=2)
  
  index_mle_param <- 2
  
  plot(draw_n[-length(draw_n)], colMeans(mle_fit_example[,-length(draw_n),index_mle_param])
       - mle_true[index_mle_param], main="First para MLE / Error and Var", ylim=c(-0.01,0.05))
  abline(h=0)
  lines(draw_n[-length(draw_n)], colMeans((mle_fit_example[,-length(draw_n),index_mle_param]
                                           -colMeans(mle_fit_example[,-length(draw_n),index_mle_param]))^2),
        col="red", lty=2, lwd=2)
}

rep.col <- function(matrix, n_copy){
  n_row <- dim(matrix)[1] 
  return(matrix(apply(matrix, MARGIN = 2,
                 function(x){replicate(n_copy,x)}), nrow=n_row))
}

# Example
rep.col(matrix(1:9, nrow=3), 2)
rep.col(matrix(1:6, nrow=2), 3)

rep.mat <- function(matrix, n_copy){
  return(matrix(rep(matrix, n_copy), nrow=nrow(matrix)))
}

# Example
rep.mat(matrix(1:9, nrow=3), 2)
rep.mat(matrix(1:6, nrow=2), 3)

NOUmatrix <- function(nw_topo, times, data, thresholds){
  # TODO check sum with difference time series
  
  # Check 
  nw_topo <- StdTopo(nw_topo = nw_topo)
  nw_test <- nw_topo
  diag(nw_test) <- 1 
  if(max(Re(eigen(nw_test)$values)) <= 0){
    stop("Q matrix should be positive definite.")
  }
  
  remove(nw_test) # delete test matrix
  
  d <- ncol(data)
  N <- nrow(data)
  
  if(length(thresholds) != d){
    stop('Wrong dimension between data and thresholds.')
  }
  if(length(times) != N){
    stop('Wrong dimensions between data and times.')
  }
  
  diff_times <- TimeMatrix(times = times, ncol = d)
  diff_filtered <- DataFiltering(data = data, diff_values = T,
                                 thresholds = thresholds)
  a_n <- colSums(rep.col(data[-N,], d) * 
                   rep.mat(diff_filtered, d))
  k_n <- matrix(colSums(rep.col(data[-N,], n_copy = d) * 
                         rep.mat(data[-N,], n_copy = d) * 
                         rep.col(diff_times[-N,], d)), d)

  a_n <- matrix(a_n, nrow=d)
  #s_n <- rep.mat(matrix = solve(k_n), n_copy = d)
  
  # TODO change with LU decomposition
  return(-apply(X = a_n, MARGIN = 2, function(x){solve(a = k_n, b = x)}))
}

# Example
if(LOAD_EXAMPLES)
{
  source("network_generation.R")
  d <- 5 #dims
  N <- 24*1500 #numbers of points
  M <- 10 # number of simulations
  Y0 <- 1 #start point
  delta_t <- 1/24
  draw_n <- seq(from=8, to=N, by=1000)

  
  sigma <- 0.1
  set.seed(42)
  nw_topo <- genRdmAssymetricGraphs(d = d, p.link = 0.25,
                                    theta_1 = 1, theta_2 = 0)
  nw_topo <- StdTopo(nw_topo)
  
  times <- seq(from = 0, by = delta_t, length.out = N)
  thresholds <- rep(0, d)

  nw_data_bm <- matrix(rnorm(n = d*N, mean = 0, sd = sigma*sqrt(delta_t)), ncol = d)
  nw_data <- matrix(0, ncol=d, nrow=N)
  nw_data[1,] <- nw_data_bm[1,]*sqrt(delta_t)
  
  nw_q <- 3*nw_topo
  diag(nw_q) <- 10
    
  for(index in 2:N){
    nw_data[index,] <-  nw_data[index-1,] - (nw_q %*% nw_data[index-1,]) * (times[index]-times[index-1]) + nw_data_bm[index,]
  }
    
  plot(nw_data[1:300,1], type="l")
  thresholds <- rep(10, d)
  mle_matrix_test <- NOUmatrix(nw_topo = nw_q, times = times, data = nw_data, thresholds = thresholds)
  
  
  # repeated testing
  M <- 1000 # number of simulations
  sigma <- 0.1
  set.seed(42)
  nw_topo <- genRdmAssymetricGraphs(d = d, p.link = 0.25,
                                    theta_1 = 1, theta_2 = 0)
  nw_topo <- StdTopo(nw_topo)
  
  times <- seq(from = 0, by = delta_t, length.out = N)
  thresholds <- rep(10, d)
  mle_mat_example <- matrix(0, ncol=d^2, nrow=M)
  for(i in 1:M){
    nw_data_bm <- matrix(rnorm(n = d*N, mean = 0, sd = sigma*sqrt(delta_t)), ncol = d)
    nw_data <- matrix(0, ncol=d, nrow=N)
    nw_data[1,] <- nw_data_bm[1,]*sqrt(delta_t)
    
    nw_q <- 3*nw_topo
    diag(nw_q) <- 10
    
    for(index in 2:N){
      nw_data[index,] <-  nw_data[index-1,] - (nw_q %*% nw_data[index-1,]) * (times[index]-times[index-1]) + nw_data_bm[index,]
    }
    
    mle_mat_example[i,] <- as.vector(
      NOUmatrix(nw_topo = nw_q, times = times, data = nw_data, thresholds = thresholds)
    )
  }
}


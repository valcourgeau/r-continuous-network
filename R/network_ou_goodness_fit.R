setwd("~/GitHub/r-continuous-network/R/")
source("R/package_to_load.R")

### Function to assess goodness-of-fit of MLE estimations
# based on Renyi entropy and statistics (Morales et al. (2000))

NOUGof <- function(theta, theta_true, data, times, quiet=F){
  # TODO test computation of \Sigma with non Id case
  N <- nrow(data)
  d <- ncol(data)
  diff_times <- TimeMatrix(times, d)
  
  theta <- as.vector(theta)
  theta_true <- as.vector(theta_true)
  
  # trying to do kronecker product with d*d matrices
  k_n <- matrix(colSums(rep.col(data[-N,], n_copy = d) * 
                          rep.mat(data[-N,], n_copy = d) * 
                          rep.col(diff_times[-N,], d)), d)
  theta_mat <- matrix(theta-theta_true, nrow = d, byrow = T)
  rv_value <- as.vector(t(theta_mat %*% k_n)) %*%(theta-theta_true)
  
  eigen_vals <- eigen(k_n, only.values = T)$value
  first <- suppressWarnings(imhof(q=rv_value, lambda = eigen_vals))
  abserr <- first$abserr
  first <- 1-first$Qq
  last <- suppressWarnings(imhof(q=rv_value, lambda = eigen_vals)$Qq)
  
  temp <- min(first, last)
  last <- max(first, last)
  first <- max(temp, 0.0) 

  probability <- abs(last-first)
  
  if(probability > 0.975){
    probability <- '***'
  }else{
    if(probability > 0.95){
      probability <- '**'
    }else{
      if(probability > 0.90){
        probability <- '*'
      }else{
        probability <- 'x'
      }
    }
  }
  test_values <- list()
  test_values[["lower"]] <- first
  test_values[["upper"]] <-last
  test_values[["rating"]] <-as.character(probability)
  test_values[["abs.error"]] <-abserr
                     
  if(!quiet){
    cat("Goodness-of-test based on Renyi statistics\n")
    print("97.5%='***' / 95%='**' / 90%='*'")
  }
  
  return(test_values)
}

# https://cran.r-project.org/web/packages/CompQuadForm/CompQuadForm.pdf
{
  source("network_generation.R")
  d <- 5 #dims
  N <- 24*1500 #numbers of points
  M <- 1 # number of simulations
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
  
  ans <- NOUGof(theta=mle_matrix_test, theta_true=nw_q, data=nw_data, times=times, quiet=T)
  ans
}
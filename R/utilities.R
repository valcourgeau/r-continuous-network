LOAD_EXAMPLES <- FALSE

# Functions to help performs technical tasks such as:
#   - Integration with discrete data;
#   - Jump tests
#   - Levy increments recovery;

JumpTest1D <- function(data, method){
  # Wrapper of JumpTest::jumptestperiod
  # TODO diff or no diff for "Amed" or "Amin" tests
  
  data <- apply(data, MARGIN = 2, diff)
  return(JumpTest::jumptestperiod(retmat = data,
                                  method = method))
}

# Example

if(LOAD_EXAMPLES){
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
    nw_data_poisson <- matrix(rpois(n = d*N, lambda = 0.1/sigma)*sqrt(delta_t), ncol = d)
    
    nw_data <- matrix(0, ncol=d, nrow=N)
    nw_data[1,] <- nw_data_bm[1,]*sqrt(delta_t)
    
    nw_q <- 0.1*nw_topo
    diag(nw_q) <- 1
    
    for(index in 2:N){
      nw_data[index,] <-  nw_data[index-1,] - 
        (nw_q %*% nw_data[index-1,]) * (times[index]-times[index-1]) +
        nw_data_bm[index,] +
        nw_data_poisson[index,]
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
  par(mfrow=c(1,1))
  plot(nw_data[,1], type='l')
  JumpTest1D(data = nw_data, method = "BNS")
}


#' Returns an offset diagonal matrix

AugmentedDiag <- function(d, offset){
  assertthat::assert_that(d > offset)
  
  mat_temp <- cbind(
    matrix(0, nrow=d, ncol=abs(offset)),
    rbind(diag(d-abs(offset)), matrix(0, ncol=d-abs(offset), nrow=abs(offset)))
  )
  if(offset < 0){
    return(t(mat_temp))
  }else{
    return(mat_temp)
  }
}

AugmentedDiag(d=10, -2)

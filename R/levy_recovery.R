# Script to implement levy increment recovery given network topology or
CHECK_EXAMPLE <- FALSE

LevyRecovery <- function(nw_topo, data, times, m=1, fitted=FALSE, on_matrix=FALSE){
   # only for OU processes
  assertthat::are_equal(length(times), nrow(data))

  if(!fitted){
    if(on_matrix){
      Q_hat <- NOUmatrix(nw_topo = nw_topo, data = data, times = times, thresholds = rep(1000, ncol(data)))
    }else{
      Q_hat <- NOUfit(nw_topo = nw_topo, data = data, times = times, thresholds = rep(1000, ncol(data)))
      Q_hat <- Q_hat[['Q']]
    }
  }else{
    Q_hat <- nw_topo
    warning('Given fitted Q, ignoring flag on_matrix;')
  }
  
  if(m!=1){
    stop('m > 1: not implemented')
  }
  
  diff_x <- apply(data, MARGIN = 2, FUN = function(x){diff(x, lag = m)})
  diff_times <- diff(times, lag =  m)
  integrated_x <- zoo::rollapply(t(t(data[-((nrow(data)-m+1):nrow(data)),])*diff_times), width = m, FUN = sum) # create integral component-wise

  # few dimension checks
  assertthat::are_equal(nrow(integrated_x), nrow(diff_x))
  assertthat::are_equal(nrow(data)-m, nrow(diff_x))
  assertthat::are_equal(dim(Q_hat)[1], ncol(integrated_x))
  
  q_integrated_x <- apply(integrated_x, MARGIN = 1, function(x){return(Q_hat %*% x)})
  if(class(Q_hat) == 'dgCMatrix'){
    q_integrated_x <- lapply(q_integrated_x, as.matrix) # getting rid of potential sparse matrix 
    q_integrated_x <- t(as.matrix(as.data.frame(q_integrated_x)))
  }else{
    q_integrated_x <- t(q_integrated_x)
  }
 
  recover <- diff_x + q_integrated_x
  return(recover)
}

# after doing re_europe_inference.R
data_tmp <- as.matrix(df_load[,-1])[1:1000,]
adj_tmp <- igraph::as_adjacency_matrix(topo_graph)
lr <- LevyRecovery(nw_topo = adj_tmp, data = data_tmp, times = (1:nrow(data_tmp))/24, m = 1, fitted = F, on_matrix = T)
hist(lr[1:800,2], breaks = 50)
hist(lr[1:800,2], breaks = 50, add=T)



if(CHECK_EXAMPLE){
  source("R/network_generation.R")
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

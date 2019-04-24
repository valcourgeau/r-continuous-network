setwd('~/GitHub/r-continuous-network/R/')
source('package_to_load.R')
source('utilities.R')
source('adjacency_generation.R')
source('path_generation.R')
source('network_ou_mle.R')
source('utilities.R')
source('bivariate_params_generation.R')

library(parallel)

set.seed(42)

### LAYOUT
n_paths <- 1000
d <- 50
theta_2 <- 0.8
theta_1 <- 0.2

# Fixing starting points etc
N <- 1000
delta_t <- 1/100 # ratio N * delta_t = 10 here
y_0 <- 0.5

# generating noise
l_data <- list()
for(i in 1:n_paths){
  l_data[[i]] <- matrix(rnorm(N*d,c(0,0.1)*sqrt(delta_t)), nrow=N)
  #l_data[[i]] <- matrix(rghyp(N*d, NIG(chi = 5, psi = 3)), nrow=N)
}

# generating theta grid
n_x <- 16
n_y <- 16
theta_grid <- CornerGrid(top_x = 4, top_y = 4, n_x = n_x, n_y = n_y, 'linear')



# choosing network type
network_types <- c(IsolatedNetwork, PolymerNetwork, LatticeNetwork, FullyConnectedNetwork)
f_network <- PolymerNetwork

PLOT_IT <- FALSE

sim_res <- list()

total <- nrow(theta_grid)
# create progress bar
pb <- txtProgressBar(min = n_x+1, max = total, style = 3, width = 50)


for(index_path in (n_x+1):nrow(theta_grid)){
  # for each couple of thetas, perform study on all paths
  theta_1 <- theta_grid[index_path,1]
  theta_2 <- theta_grid[index_path,2]
  setTxtProgressBar(pb, index_path)
  
  #cat('theta_1', theta_1, 'theta_2', theta_2, '\n')
  network_topo <- f_network(d = d, theta_1 = theta_1, theta_2 = theta_2)
  
  res <- mclapply(X = l_data, FUN = 
                    function(x){
                      ConstructPath(nw_topo = network_topo, noise = x, delta_time = delta_t, y_0)
                    }
  )
  
  res_nou <- mclapply(X=res, FUN = (function(x){NOUfit(nw_topo = network_topo, 
                                                       times = seq(0, length.out = N, by = delta_t), 
                                                       data = x, thresholds = rep(10, d))}))
  res_mles <- vapply(res_nou, function(x){x$MLE_wide}, c(1,1)) %>% t
  
  sim_res[[index_path]] <- res_mles
  
  if(PLOT_IT){
    vioplot(res_mles[,1], res_mles[,2], names=c(expression(theta[1]), expression(theta[2])))
    abline(h=theta_grid[index_path,1], col='red')
    abline(h=theta_grid[index_path,2], col='blue')
  
    cat('MSE Theta_1', sum((res_mles[,1] - theta_1)^2),' \n')
    cat('MSE Theta_2', sum((res_mles[,2] - theta_2)^2),' \n')
  }
}
close(pb)

sim_res$N <- N
sim_res$delta_t <- delta_t
sim_res$n_x <- n_x
sim_res$n_y <- n_y
sim_res$y_0 <- y_0
sim_res$n_paths <- n_paths

rlist::list.save(sim_res, 'Polymer_4_4.RData')
sim_mean <- matrix(0, ncol=2, nrow=(n_x-1)*n_y)
sim_bias <- matrix(0, ncol=2, nrow=(n_x-1)*n_y)
sim_bias_percent <- matrix(0, ncol=2, nrow=(n_x-1)*n_y)
sim_sd <- matrix(0, ncol=2, nrow=(n_x-1)*n_y)

theta_mat <- as.matrix(theta_grid)

for(index in (n_x+1):nrow(theta_mat)){
  sim_mean[index-n_x,] <- apply(as.matrix(sim_res[[index]]), FUN = mean, MARGIN = 2)
  sim_bias[index-n_x,1:2] <- sim_mean[index-n_x,1:2] - theta_mat[index,1:2]
  sim_bias_percent[index-n_x,] <- sim_bias[index-n_x,1:2] / theta_mat[index,] * 100
  sim_sd[index-n_x,] <- apply(as.matrix(sim_res[[index]]), FUN = function(x){sqrt(sum(x^2)/(sim_res$n_paths-1))}, MARGIN = 2)
}

plot(theta_grid$x[order(theta_grid$x[-(1:n_x)])], sim_bias_percent[,1][order(theta_grid$x[-(1:n_x)])], 
     ylab = 'Relative bias (in %)', xlab=expression(theta[1]))
plot(theta_grid$y[order(theta_grid$y[-(1:n_x)])], sim_bias_percent[,2][order(theta_grid$y[-(1:n_y)])], 
     ylab = 'Relative bias (in %)', xlab=expression(theta[2]))

plot(theta_grid$x[order(theta_grid$x[-(1:n_x)])], sim_sd[,1][order(theta_grid$x[-(1:n_x)])], 
     ylab = 'Relative bias (in %)', xlab=expression(theta[1]))
plot(theta_grid$y[order(theta_grid$y[-(1:n_x)])], sim_sd[,2][order(theta_grid$y[-(1:n_y)])], 
     ylab = 'MC Standard deviation', xlab=expression(theta[2]))


plot(theta_grid$x[order(theta_grid$x[-(1:n_x)])], sim_bias[,1][order(theta_grid$x[-(1:n_x)])], 
     ylab = 'Relative bias (in %)', xlab=expression(theta[1]))
plot(theta_grid$y[order(theta_grid$y[-(1:n_x)])], sim_bias[,2][order(theta_grid$y[-(1:n_y)])], 
     ylab = 'MC Standard deviation', xlab=expression(theta[2]))

contour_mat <- matrix(0, ncol=n_x, nrow=n_y)
for(i in 1:(n_y-1)){
  contour_mat[i,] <- sim_bias[((i-1)*n_x+1):(i*n_x),2]
}
contour(contour_mat[-n_y,], xlab = 't2', ylab='t1')











     
par(mfrow=c(1,2))
plot((cumsum(res_mles[,1])/1:n_paths), log='xy', type='l', ylim=c(0.15,0.25),
     main=paste(n_paths, 'paths;', N, 'points;', d, 'nodes;', 'FC type;'), xlab=expression(theta[1]))
abline(h=theta_1, col = 'red', lty = 2)

plot((cumsum(res_mles[,2])/1:n_paths), log='xy', type='l', ylim=c(0.9*theta_2,1.1*theta_2),
     main=paste(n_paths, 'paths;', N, 'points;', d, 'nodes;', 'FC type;'), xlab=expression(theta[2]))
abline(h=theta_2, col = 'red', lty = 2)



line(log(1:n_paths), cumsum(res_mles[,1])/1:n_paths)
line(log(1:n_paths), cumsum(res_mles[,2])/1:n_paths)

library(ghyp)


library(vioplot)

par(mfrow=c(1,1))
vioplot(res_mles[,1], res_mles[,2], names=c(expression(theta[1]), expression(theta[2])))
title('0.2                                                                 0.8')



library(gplots)
filled.contour(res_nou[[1]]$Q, main=paste('Q MLE with theta_2 =', theta_2,' theta_1', theta_1))

network_topo


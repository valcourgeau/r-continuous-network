setwd('~/GitHub/r-continuous-network/R/')
source('package_to_load.R')
source('utilities.R')
source('adjacency_generation.R')
source('path_generation.R')
source('network_ou_mle.R')
source('utilities.R')
source('bivariate_params_generation.R')

load('~/GitHub/r-continuous-network/data/standard-networks/adj_20_nodes')
load('~/GitHub/r-continuous-network/data/standard-networks/adj_39_nodes')

set.seed(42)

### LAYOUT
n_paths <- 2
network_topo <- adj_20
d <- nrow(network_topo)
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

PLOT_IT <- FALSE

# generating theta grid
n_x <- 16
n_y <- 16
theta_grid <- CornerGrid(top_x = 4, top_y = 4, n_x = n_x, n_y = n_y, 'linear')

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
  #network_topo <- f_network(d = d, theta_1 = theta_1, theta_2 = theta_2)
  diag(network_topo) <- 0
  network_topo <- network_topo / rowSums(network_topo)
  network_topo[is.nan(network_topo)] <- 0.0
  network_topo <- theta_1 * network_topo
  diag(network_topo) <- theta_2
  
  res <- mclapply(X = l_data, FUN = 
                    function(x){
                      ConstructPath(nw_topo = network_topo, noise = x, delta_time = delta_t, y_0)
                    })
  
  res_nou <- mclapply(X=res, FUN = (function(x){NOUfit(nw_topo = network_topo, 
                                                       times = seq(0, length.out = N, by = delta_t), 
                                                       data = x, thresholds = rep(10, d))}))
  res_mles <- vapply(res_nou, function(x){x$MLE_wide}, c(1,1)) %>% t
  
  sim_res[[index_path]] <- res_mles
  
  if(PLOT_IT){
    vioplot::vioplot(res_mles[,1], res_mles[,2], names=c(expression(theta[1]), expression(theta[2])))
    abline(h=theta_grid[index_path,1], col='red')
    abline(h=theta_grid[index_path,2], col='blue')
    
    cat('MSE Theta_1', sum((res_mles[,1] - theta_1)^2),' \n')
    cat('MSE Theta_2', sum((res_mles[,2] - theta_2)^2),' \n')
  }
}
close(pb)

sim_res$thetas <- theta_grid
sim_res$N <- N
sim_res$delta_t <- delta_t
sim_res$n_x <- n_x
sim_res$n_y <- n_y
sim_res$y_0 <- y_0
sim_res$n_paths <- n_paths

rlist::list.save(sim_res, 'IEEE_20_nodes_BrownianMotion_4_4.RData')


set.seed(42)

### LAYOUT
n_paths <- 500
network_topo <- adj_20
d <- nrow(network_topo)
theta_2 <- 0.8
theta_1 <- 0.2

n_substeps <- 12 # number of sub-training performed

# Fixing starting points etc
N <- 2^n_substeps
delta_t <- 1/100 # ratio N * delta_t = 10 here
y_0 <- 0.5



# generating noise
l_data <- list()
for(i in 1:n_paths){
  l_data[[i]] <- matrix(rnorm(N*d,c(0,0.1)*sqrt(delta_t)), nrow=N)
  #l_data[[i]] <- matrix(rghyp(N*d, NIG(chi = 5, psi = 3)), nrow=N)
}

PLOT_IT <- FALSE

# generating theta grid
n_x <- 16
n_y <- 16
theta_grid <- CornerGrid(top_x = 4, top_y = 4, n_x = n_x, n_y = n_y, 'linear')

sim_res <- list()

total <- nrow(theta_grid)
# create progress bar
pb <- txtProgressBar(min = n_x+1, max = total, style = 3, width = 50)

for(index_path in (n_x+1):nrow(theta_grid)){
  tt <- proc.time()
  # for each couple of thetas, perform study on all paths
  theta_1 <- theta_grid[index_path,1]
  theta_2 <- theta_grid[index_path,2]
  setTxtProgressBar(pb, index_path)
  
  #cat('theta_1', theta_1, 'theta_2', theta_2, '\n')
  #network_topo <- f_network(d = d, theta_1 = theta_1, theta_2 = theta_2)
  network_topo <- adj_20
  diag(network_topo) <- 0
  network_topo <- network_topo / rowSums(network_topo)
  network_topo[is.nan(network_topo)] <- 0.0
  network_topo <- theta_1 * network_topo
  diag(network_topo) <- theta_2
  
  res <- mclapply(X = l_data, FUN = 
                    function(x){
                      ConstructPath(nw_topo = network_topo, noise = x, delta_time = delta_t, y_0)
                    })
  
  #res <- res[[1]]
  subs <- 2^(1:n_substeps)
  sim_res[[index_path]] <- vapply(
           res, 
           function(x){
                res_nou <- mclapply(subs, 
                                    FUN = (function(k_final){
                                                nou_fit <- NOUfit(nw_topo = network_topo, 
                                                                  times = seq(0, length.out = k_final, by = delta_t), 
                                                                  data = x[1:k_final,],
                                                                  thresholds = rep(10, d))
                                                return(nou_fit)}
                                           ))
                res_mles <- vapply(res_nou, function(x){x$MLE_wide}, c(1,1)) %>% t
                return(res_mles)
              },
              matrix(0, ncol=2, nrow=n_substeps))
  #sim_res[[index_path]] <- res_mles
  cat('Iteration in', (proc.time() - tt)[3], 's.\n')
}
close(pb)

sim_res$substeps <- 2^(1:n_substeps)
sim_res$thetas <- theta_grid
sim_res$N <- N
sim_res$delta_t <- delta_t
sim_res$n_x <- n_x
sim_res$n_y <- n_y
sim_res$y_0 <- y_0
sim_res$n_paths <- n_paths
sim_res$n_substeps <- n_substeps

rlist::list.save(sim_res, 'IEEE_20_nodes_BrownianMotion_4_4_subtraining.RData')

pb <- txtProgressBar(min = 1, max = sim_res$n_paths, style = 3, width = 50)
bias_array_full <- array(dim = c(nrow(sim_res$thetas), sim_res$n_paths, sim_res$n_substeps, 2))
for(path_number in 1:sim_res$n_paths){
  bias_array <- array(dim = c(nrow(sim_res$thetas), sim_res$n_substeps, 2))
  for(index_path in (n_x+1):nrow(theta_grid)){
    bias_list <- apply(sim_res[[index_path]][,,path_number], function(x){abs(x - sim_res$thetas[index_path,1:2])}, MARGIN = 1)
    bias_mat <- matrix(0, nrow = sim_res$n_substeps, ncol = 2)
    i <- 1
    for(bias_indiv in bias_list){
      bias_mat[i,] <- c(bias_indiv$x, bias_indiv$y)
      i <- i + 1
    }
    bias_array[index_path,,] <- bias_mat
  }
  bias_array_full[,path_number,,] <- bias_array[,,]
  setTxtProgressBar(pb, path_number)
}

# bias_array_full is thetas, path, subsampling, choice of theta
# write.table(bias_array_full, 'bias_array_full')

regression_res <- rep(0, length(n_x*n_y))
for(i in (n_x+1):(n_x*n_y)){
  # bias_array_full[20,1,2:12,1] %>% plot
  sd_ba <- vapply(1:12, function(slice){bias_array_full[i,,slice,1] %>% sd}, 1)
  mean_ba <- vapply(1:12, function(slice){bias_array_full[i,,slice,1] %>% mean}, 1)
  
  plot(mean_ba[2:12]  %>% log, type= 'b', xaxt="n", xlab='Sample size', ylab='log-Bias', 
       cex.lab=1.5, cex.axis=1.3, main='IEEE 20-Nodes bias regression')
  axis(1, at=1:11, labels=c(expression('2'^1), expression('2'^2), expression('2'^3), expression('2'^4), expression('2'^5),
                         expression('2'^6), expression('2'^7), expression('2'^8), expression('2'^9), expression('2'^10),
                         expression('2'^11)), cex.axis=1.3)
  
  lines((mean_ba[2:12] + 1.96*sd_ba[2:12] / sqrt(sim_res$n_paths)) %>% log, type= 'b', lty=2)
  lines((mean_ba[2:12] - 1.96*sd_ba[2:12] / sqrt(sim_res$n_paths)) %>% log, type= 'b', lty=2)
  regression_res[i] <- (line(mean_ba[2:12])$coefficients[2])
}

hist(regression_res[17:256], breaks=100)

for(i in 2:500){
  lines(bias_array_full[20,i,2:12,1] )
}


plot(2^(2:sim_res$n_substeps), bias_array[sim_res$n_x+1,2:sim_res$n_substeps,2], log='xy', type='b', ylab='Absolute bias',
     xlab='N', cex.lab = 1.4, las=1,xaxt="n",yaxt="n")

df <- cbind(2^(2:sim_res$n_substeps), bias_array[index_path, 2:sim_res$n_substeps,2]) %>% as.data.frame()
library('scales')
library('ggplot2')
ggplot(df, aes(x = V1, y = V2)) + geom_point() + scale_x_continuous(trans = log2_trans()) + scale_y_continuous(trans = log2_trans())

pp <- ggplot()
df <- cbind(2^(2:17), bias_array[, 2:sim_res$n_substeps,2]) %>% as.data.frame()
ggplot(df, aes(x = V1, y = V2)) + geom_point() + scale_x_continuous(trans = log2_trans()) + scale_y_continuous(trans = log2_trans())

coeffs <- rep(0, nrow(sim_res$thetas))
for(index_path in (n_x+1):nrow(sim_res$thetas)){
  coeffs[index_path] <- line(x=1:sim_res$n_substeps, y = bias_array[index_path,1:sim_res$n_substeps,1] %>% log)$coefficients[2]
  #lines(2:sim_res$n_substeps, bias_array[index_path, 2:sim_res$n_substeps,2] %>% (function(x){log(x,2)}), type='b')
}

plot(theta_grid$x[order(theta_grid$x[-(1:n_x)])], coeffs[order(theta_grid$x[-(1:n_x)])], xlab=expression(theta[1]))
plot(theta_grid$y[order(theta_grid$y[-(1:n_x)])], coeffs[order(theta_grid$y[-(1:n_x)])], xlab=expression(theta[2]))

plot(theta_grid[(4*n_x+1):(5*n_x),1],coeffs[(4*n_x+1):(5*n_x)])


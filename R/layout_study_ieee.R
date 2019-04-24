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
n_paths <- 1000
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

rlist::list.save(sim_res, 'IEEE_20_nodes_BrownianMotion_4_4.RData')

setwd('~/GitHub/r-continuous-network/R/')
source('package_to_load.R')
source('utilities.R')
source('adjacency_generation.R')
source('path_generation.R')
source('network_ou_mle.R')

set.seed(42)

### LAYOUT
n_paths <- 500
d <- 50
theta_2 <- 0.8
theta_1 <- 0.2
network_types <- c(IsolatedNetwork, PolymerNetwork, FullyConnectedNetwork)
f_network <- PolymerNetwork

network_topo <- f_network(d = d, theta_1 = theta_1, theta_2 = theta_2)
N <- 1000
delta_t <- 1/100
y_0 <- 0.5

l_data <- list()
for(i in 1:n_paths){
   l_data[[i]] <- matrix(rnorm(N*d,c(0,0.1)*sqrt(delta_t)), nrow=N)
  #l_data[[i]] <- matrix(rghyp(N*d, NIG(chi = 5, psi = 3)), nrow=N)
}

library(parallel)
res <- mclapply(X = l_data, FUN = 
                  function(x){
                    ConstructPath(nw_topo = network_topo, noise = x, delta_time = delta_t, y_0)
                  }
                )
# res <- lapply(X = l_data,  function(x){
#   ConstructPath(nw_topo = network_topo, noise = x, delta_time = 1/24, 0.5)
# })

res_nou <- mclapply(X=res, FUN = (function(x){NOUfit(nw_topo = network_topo, 
                                             times = seq(0, length.out = N, by = delta_t), 
                                             data = x, thresholds = rep(10, d))}))
res_mles <- vapply(res_nou, function(x){x$MLE_wide}, c(1,1)) %>% t

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


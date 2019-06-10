setwd('~/GitHub/r-continuous-network/R/')
source('package_to_load.R')
source('utilities.R')
source('adjacency_generation.R')
source('path_generation.R')
source('network_ou_mle.R')
source('utilities.R')
source('bivariate_params_generation.R')

### LAYOUT
setwd('~/GitHub/r-continuous-network/data/')
data_50_re <- rlist::list.load('data_re_europe_50.RData')$data
adj_50_re <- rlist::list.load('adjacency_re_europe_50.RData')$adj
Q_used <- rlist::list.load('mle_re_reurope_50.RData')
ghyp_fit_re_europe_50 <- rlist::list.load('50_nodes_ghyp_fit_v3.RData')

setwd('~/GitHub/r-continuous-network/data/sim-study/')
n_paths <- 1000
d <- 50
theta_2 <- Q_used$MLE_wide[2]
theta_1 <- Q_used$MLE_wide[1]

# Fixing starting points etc
N <- 3000
delta_t <- 2/24 # ratio N * delta_t = 10 here
y_0 <- (data_50_re)[1,]

ghyp_full <- rlist::list.load('~/GitHub/r-continuous-network/data/ghyp_in_full_0_re.RData')
ghyp_full <- ghyp_full[[1]]
# generating noise
l_data <- list()#array(0, c(n_paths, N, d))
l_data_finite <- list()
l_data_infinite <- list()


par(mfrow=c(2,5), mar=c(4,4.5,2,4))
cex_value <- 1.8
for(i in c(10,40)){
  i <- round(i)
  print(i)
  
  quantil_vals <- 1:4999/5000
  
  val_x_qq <- quantile(x = lr$increments[,i], quantil_vals)
  set.seed(42)
  val_y_qq <- ghyp::rghyp(object = ghyp_fit_re_europe_50$NIG, n = length(lr$increments[,i]))[, i]
  val_y_qq <- quantile(x = val_y_qq, quantil_vals)
  qqplot(x = val_x_qq, y = val_y_qq, 
         ylab = 'NIG', xlab= paste('Levy increments', i), 
         cex.lab = cex_value, cex.axis=cex_value)

  if(i == 10){title('NIG', cex.main=cex_value, line = 0.8)}
  
  set.seed(42)
  val_y_qq <- ghyp::rghyp(object = ghyp_fit_re_europe_50$GAUSS, n = length(lr$increments[,i]))[, i]
  val_y_qq <- quantile(x = val_y_qq, quantil_vals)
  qqplot(x = val_x_qq, y = ghyp::rghyp(object = ghyp_fit_re_europe_50$GAUSS, n = length(lr$increments[,i]))[, i], 
         ylab = 'Gaussian', xlab= paste('Levy increments', i), 
         cex.lab = cex_value, cex.axis=cex_value)
  if(i == 10){title('Gaussian', cex.main=cex_value, line = 0.8)}
  
  set.seed(42)
  val_y_qq <- ghyp::rghyp(object = ghyp_fit_re_europe_50$VG, n = length(lr$increments[,i]))[, i]
  val_y_qq <- quantile(x = val_y_qq, quantil_vals)
  qqplot(x = val_x_qq, y = ghyp::rghyp(object = ghyp_fit_re_europe_50$VG, n = length(lr$increments[,i]))[, i], 
         ylab = 'VG', xlab= paste('Levy increments', i), 
         cex.lab = cex_value, cex.axis=cex_value)
  if(i == 10){title('VG', cex.main=cex_value, line = 0.8)}
  
  set.seed(42)
  val_y_qq <- ghyp::rghyp(object = ghyp_fit_re_europe_50$T, n = length(lr$increments[,i]))[, i]
  val_y_qq <- quantile(x = val_y_qq, quantil_vals)
  qqplot(x = val_x_qq, y = ghyp::rghyp(object = ghyp_fit_re_europe_50$T, n = length(lr$increments[,i]))[, i], 
         ylab = 'Student\'s t', xlab= paste('Levy increments', i),
         cex.lab = cex_value, cex.axis=cex_value)
  if(i == 10){title('Student\'s t', cex.main=cex_value, line = 0.8)}
  
  set.seed(42)
  val_y_qq <- ghyp::rghyp(object = ghyp_full, n = length(lr$increments[,i]))[, i]
  val_y_qq <- quantile(x = val_y_qq, quantil_vals)
  qqplot(x = val_x_qq, y = ghyp::rghyp(object = ghyp_full, n = length(lr$increments[,i]))[, i], 
         ylab = 'GHYP', xlab= paste('Levy increments', i), 
         cex.lab = cex_value, cex.axis=cex_value)
  if(i == 10){title('GHYP', cex.main=cex_value, line = 0.8)}
}


set.seed(40)
for(i in 1:n_paths){
  # l_data[[i]] <- matrix(rnorm(N*d,c(0,0.1)*sqrt(delta_t)), nrow=N)
  l_data[[i]]<- matrix(ghyp::rghyp(n = N, object = ghyp_full), nrow=N)
  # l_data_finite[[i]] <- matrix(rpois(N*d, lambda = 2/24) * rnorm(N*d, mean = 0, sd = 0.2), ncol=d)
  # l_data_infinite <-
}

load('~/GitHub/r-continuous-network/data/standard-networks/adj_20_nodes')
load('~/GitHub/r-continuous-network/data/standard-networks/adj_39_nodes')

# choosing network type
network_types <- c(PolymerNetwork, LatticeNetwork, FullyConnectedNetwork, adj_50_re)
network_types_name <- c('polymer_50_mles', 'lattice_50_mles', 'fc_50_mles', '50re_50_mles')
index_network <- 1
network_study <- list()

for(f_network in network_types){
  print(index_network)
  if(index_network < 4){
    network_topo <- f_network(d = d, theta_1 = theta_1, theta_2 = theta_2)
  }else{
    network_topo <- network_types[index_network][[1]] %>% as.matrix
    diag(network_topo) <- 0
    network_topo <- network_topo / rowSums(network_topo)
    network_topo[is.nan(network_topo)] <- 0
    network_topo <- theta_1 * network_topo
    diag(network_topo) <- theta_2
  }
  if(index_network != 3){
    network_topo <- as(network_topo, 'CsparseMatrix')
  }
  generated_paths <- lapply(X = l_data, FUN =
                              function(x){
                                ConstructPath(nw_topo = network_topo, noise = x, delta_time = delta_t, y_0)
                              })
  generated_paths<- rlist::list.save(generated_paths, paste(network_types_name[index_network], '.RData', sep = ''))
  # generated_paths <- rlist::list.load(paste(network_types_name[index_network], '.RData', sep = ''))
  print('paths generated')
  
  # starting from topology
  #network_topo[which(abs(network_topo) > 1e-16)] <- 1
  generated_fit <- lapply(X = generated_paths, FUN = 
                            function(x){
                              NOUfit(nw_topo = network_topo, times = seq(0, length.out = nrow(generated_paths[[1]]), by = delta_t), data = x, thresholds = rep(1000, d))$MLE_wide
                            })
  print('fit generated')
  gen_fit_matrix <- matrix(unlist(generated_fit), ncol = 2, byrow = T)
  network_study[[paste('network_',index_network,'_t1', sep = '')]] <- gen_fit_matrix[,1]
  network_study[[paste('network_',index_network,'_t2', sep = '')]] <- gen_fit_matrix[,2]
  index_network <- index_network + 1
}

rlist::list.save(network_study, 'network_study_results.RData')
colors <- c('#08605F', '#177E89', '#598381', '#8E936D', '#A2AD59')
colors2 <- c('#E88D67', '#CA7AB8', '#9999C3', '#7B8CDE')
par(mfrow=c(2,1), mar=c(2,3,1.8,0.5))
vioplot::vioplot(network_study[c(7,1,3,5)], names = c('RE-Europe 50', 'Polymer', 'Lattice', 'Complete'),
                 col = colors2,
                 # col = c('#BE2B2B', '#2BBEBE', '#2BBEBE', '#2BBEBE'),
                 main=expression(theta[1]),
                 cex.main=1.8, cex.names=1.5, cex.axis=1.5)
par(xpd = F) 
abline(h=theta_1, lwd=2, lty=2)
vioplot::vioplot(network_study[c(8,2,4,6)], names = c('RE-Europe 50', 'Polymer', 'Lattice', 'Complete'),
                 col =colors2,
                 # col = c('#BE2B2B', '#2BBEBE', '#2BBEBE', '#2BBEBE'),
                 main=expression(theta[2]),
                 cex.main=1.8, cex.names=1.5, cex.axis=1.5)
abline(h=theta_2, lwd=2, lty=2)


# estimate MC SD


list(, gen_fit_matrix[,2])
index_network <- index_network + 1

ghyp::lik.ratio.test(ghyp_full, ghyp_fit_re_europe_50$NIG, conf.level = 0.99)

tseries::kpss.test(x = lr$increments[, 49])
tseries::kpss.test(x = lr$increments[, 9])
tseries::kpss.test(x = lr$increments[, 20])

tseries::kpss.test(x = lr$increments[, 49], null = 'Trend')
tseries::kpss.test(x = lr$increments[, 9], null = 'Trend')
tseries::kpss.test(x = lr$increments[, 20], null = 'Trend')

Box.test(x = lr$increments[, 49], lag = 10, type='Lju')
Box.test(x = lr$increments[, 9], lag = 1)
Box.test(x = lr$increments[, 20], lag = 1)

tseries::adf.test(lr$increments[,49])
tseries::adf.test(data_tmp[, 30])

par(mfrow=c(1,1))
vioplot::vioplot()

par(mfrow=c(1,2), mar=c(4,4,2,2))
forecast::Acf(lr$increments[1:500,10], 10, main='', xlab='', ylab='',
              lwd = 3, cex = 1.3, cex.axis = 1.0)
title("Node 10", line = 1)
title(xlab="Lag", line=2.3, cex.lab=1.4)
title(ylab="ACF", line=2.3, cex.lab=1.4)
forecast::Acf(lr$increments[1:500,40], 10, main='', xlab='', ylab='',
              lwd = 3, cex = 1.3, cex.axis = 1.0)
title("Node 40", line = 1)
title(xlab="Lag", line=2.3, cex.lab=1.4)
title(ylab="ACF", line=2.3, cex.lab=1.4)

# # generating theta grid
# n_x <- 16
# n_y <- 16
# theta_grid <- CornerGrid(top_x = 4, top_y = 4, n_x = n_x, n_y = n_y, 'linear')

PLOT_IT <- FALSE

sim_res <- list()
sim_tmp <- list()


set.seed(40)
for(i in 1:n_paths){
  # l_data[[i]] <- matrix(rnorm(N*d,c(0,0.1)*sqrt(delta_t)), nrow=N)
  l_data[[i]]<- matrix(rghyp(n = N, object = ghyp_full), nrow=N)
  # l_data_finite[[i]] <- matrix(rpois(N*d, lambda = 2/24) * rnorm(N*d, mean = 0, sd = 0.2), ncol=d)
  # l_data_infinite <-
}

n_paths <- 1000
d <- 50
theta_2 <- Q_used$MLE_wide[2]
theta_1 <- Q_used$MLE_wide[1]

# Fixing starting points etc
N <- 2^16
delta_t <- 2/24 # ratio N * delta_t = 10 here
y_0 <- (data_50_re)[1,]


network_topo <- adj_50_re %>% as.matrix
diag(network_topo) <- 0
network_topo <- network_topo / rowSums(network_topo)
network_topo[is.nan(network_topo)] <- 0.0
network_topo <- theta_1 * network_topo
diag(network_topo) <- theta_2

set.seed(41)
sample_study <- list()
generated_paths_N <- lapply(1:500,
                          function(i){
                              x <- rghyp(n = N, object = ghyp_full)
                              constr_path <- ConstructPath(nw_topo = network_topo, noise = x, delta_time = delta_t, y_0)
                              if(i %% 100 == 0){
                                print(i)
                              }
                              return(lapply(7:log(N, base = 2),
                                     FUN =  function(k){
                                       return(NOUfit(nw_topo = network_topo, 
                                                  times = seq(0, length.out = 2^k, by = delta_t),
                                                  data = constr_path[1:(2^k),],
                                                  thresholds = rep(1000, d))$MLE_wide)
                                     }))
                            })

rlist::list.save(generated_paths_N, paste('sample_study_adj_50_re_europe.RData'))
generated_paths_N <- rlist::list.load('~/GitHub/r-continuous-network/data/sample_study_adj_50_re_europe.RData')

generated_results_matrices_N <- lapply(generated_paths_N, 
       FUN = function(x){x %>% unlist %>% (function(x){matrix(x, ncol=2, byrow = T)})})

stack_results_matrices <- vapply(FUN = function(x){array(x, dim = c(1, dim(x)))}, 
                                 X = generated_results_matrices_N, 
                                 FUN.VALUE = matrix(data = generated_results_matrices_N[[1]], 
                                                    ncol = 2))
stack_results_matrices <- aperm(stack_results_matrices, c(3,1,2))
stack_results_matrices[,,]
mean_results <- apply(stack_results_matrices, MARGIN = 2,
                      FUN = function(x){apply(x, MARGIN = 2,
                                              FUN = mean)}) %>% t
mean_results
sd_results <- apply(stack_results_matrices, MARGIN = 2,
                    FUN = function(x){apply(x, MARGIN = 2,
                                            FUN = sd)}) %>% t
print(sd_results)

sample_study$mean_results <- mean_results
sample_study$sd_results <- sd_results
sample_study$stack_results_matrices <- stack_results_matrices
sample_study$call <- 'apply(stack_results_matrices, MARGIN = 2,
                      FUN = function(x){apply(x, MARGIN = 2,
                                              FUN = mean)}) %>% t'
sample_study$logN <- 7:log(N, base = 2)
rlist::list.save(sample_study, 'sample_study_matrices.RData')
sample_study <- rlist::list.load('~/GitHub/r-continuous-network/data/sample_study_matrices.RData')

# Fitting lines
par(mfrow=c(1,1), mar=c(4,5,2,1))
log_sd <- sample_study$sd_results %>% (function(x){log(x,base=2)})
ll <- line(x=sample_study$logN , y = log_sd[,1], iter = 1)
plot(sample_study$logN, log_sd[,1], main='',
     ylab='Parameter value (log)', xlab='', ylim=c(-10,-1),
     cex.axis=1.4, cex.lab = 1.4, cex = 2, cex.main=1.8)
title(xlab=paste('Sample Size (log)'), line=2.3, cex.lab=1.4)
lines(sample_study$logN, ll$fitted.values, lwd=3, col='#177E89', lty=4)

ll <- line(x=sample_study$logN , y = log_sd[,2], iter = 1)
points(sample_study$logN, log_sd[,2], main=expression(theta[2]),
     ylab='Parameter value (log)', xlab='',
     cex.axis=1.4, cex.lab = 1.4, cex = 2, cex.main=1.8)
title(xlab=paste('Sample Size (log)'), line=2.3, cex.lab=1.4)
lines(sample_study$logN, ll$fitted.values, lwd=3, col='#BD632F',lty=2)

legend(14.9, -1, legend=c(expression(paste(theta[1], ' ')), expression(paste(theta[2], ' '))),
       col=c('#177E89', '#BD632F'), lty=c(4,2), lwd=4, cex=2)




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
  
  sim_tmp[[index_path]] <- res_mles
  
  if(PLOT_IT){
    vioplot(res_mles[,1], res_mles[,2], names=c(expression(theta[1]), expression(theta[2])))
    abline(h=theta_grid[index_path,1], col='red')
    abline(h=theta_grid[index_path,2], col='blue')
  
    cat('MSE Theta_1', sum((res_mles[,1] - theta_1)^2),' \n')
    cat('MSE Theta_2', sum((res_mles[,2] - theta_2)^2),' \n')
  }
}
close(pb)

sim_res$sims <- sim_tmp
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


setwd('~/GitHub/r-continuous-network/R/')
polym <- rlist::list.load('Polymer_4_4.RData')
lat <- rlist::list.load('Lattice_4_4.RData')
fc <- rlist::list.load('FullyConnected_4_4.RData')

sim_res_temp <- polym
par(mfrow=c(6,2))

sim_mean <- matrix(0, ncol=2, nrow=(n_x-1)*n_y)
sim_bias <- matrix(0, ncol=2, nrow=(n_x-1)*n_y)
sim_bias_sd <- matrix(0, ncol=2, nrow=(n_x-1)*n_y)
sim_bias_percent <- matrix(0, ncol=2, nrow=(n_x-1)*n_y)
sim_sd <- matrix(0, ncol=2, nrow=(n_x-1)*n_y)

theta_mat <- as.matrix(theta_grid)

for(index in (n_x+1):nrow(theta_mat)){
  sim_mean[index-n_x,] <- apply(as.matrix(sim_res_temp[[index]]), FUN = mean, MARGIN = 2)
  sim_bias[index-n_x,1:2] <- apply(as.matrix(sim_res_temp[[index]]) - 
                                matrix(rep(theta_mat[index,1:2], times=nrow(as.matrix(sim_res_temp[[index]]))), byrow = T, ncol=2), mean, MARGIN = 2)
  sim_bias_sd[index-n_x,1:2] <- apply(as.matrix(sim_res_temp[[index]]) - 
                                     matrix(rep(theta_mat[index,1:2], times=nrow(as.matrix(sim_res_temp[[index]]))), byrow = T, ncol=2), sd, MARGIN = 2)
    # abs(sim_mean[index-n_x,1:2] - theta_mat[index,1:2])
  sim_bias_percent[index-n_x,] <- sim_bias[index-n_x,1:2] / theta_mat[index,] * 100
  sim_sd[index-n_x,] <- apply(as.matrix(sim_res_temp[[index]]), FUN = function(x){sd(x)}, MARGIN = 2)
}

plot(theta_grid[17:256,1][order(theta_grid[17:256,1])], sim_bias[,1])

par(mar=c(4, 5, 1, 1), mfrow=c(3,2) )
# layout(matrix(c(1,3,2,4,5,6), 3, 2, byrow = TRUE))

# plot(theta_grid$x[order(theta_grid$x[-(1:n_x)])], sim_bias_percent[,1][order(theta_grid$x[-(1:n_x)])], 
#      ylab = 'Relative bias (in %)', xlab=expression(theta[1]), cex.lab=1.5)
# plot(theta_grid$y[order(theta_grid$y[-(1:n_x)])], sim_bias_percent[,2][order(theta_grid$y[-(1:n_y)])], 
#      ylab = 'Relative bias (in %)', xlab=expression(theta[2]), cex.lab=1.5)

plot(theta_grid$x[order(theta_grid$x[-(1:n_x)])], sim_sd[,1][order(theta_grid$x[-(1:n_x)])], 
     ylab = 'MC std dev', xlab=expression(theta[1]), cex.lab=1.5)
plot(theta_grid$y[order(theta_grid$y[-(1:n_x)])], sim_sd[,2][order(theta_grid$y[-(1:n_y)])], 
     ylab = 'MC std dev', xlab=expression(theta[2]), cex.lab=1.5)


plot(theta_grid$x[order(theta_grid$x[-(1:n_x)])], sim_bias[,1][order(theta_grid$x[-(1:n_x)])], 
     ylab = 'Bias', xlab=expression(theta[1]), cex.lab=1.5)
plot(theta_grid$y[order(theta_grid$y[-(1:n_x)])], sim_bias[,2][order(theta_grid$y[-(1:n_y)])], 
     ylab = 'Bias', xlab=expression(theta[2]), cex.lab=1.5)


order_theta_1 <- order(theta_grid[,1])
order_theta_2 <- order(theta_grid[,2])

plot(o_theta_val, mu_centered, xlim=c(-4,4), ylim=c(-1e-2,1e-2), col = 'white')
already_done <- rep(0 , length(order_theta_1))
for(o_theta in order_theta_1){
  o_theta_val <- theta_grid[o_theta,1]
  mu_centered <- mean(sim_res_temp[[o_theta]][,1]) - o_theta_val
  sigma <- sd(sim_res_temp[[o_theta]][,1])
  if(!(o_theta_val %in% already_done)){
    points(o_theta_val, mu_centered)
    already_done[o_theta] <- o_theta_val
  }
}

ggplot(sim_bias[,1] %>% as.data.frame, aes(x=theta_grid[17:256,1], y=theta_grid[17:256,2]) ) +
  geom_hex(bins = 70) +
  theme_bw()


data_to_plot <- data.frame(zz=sim_bias[,1], xx=theta_grid[17:256,1]+ rnorm(mean=0,sd = 0.001, n = 240), yy=theta_grid[17:256,2] + rnorm(mean=0,sd = 0.001, n = 240))

# Area + contour
ggplot(data_to_plot, aes(xx,yy, z=zz)) +
  # stat_contour(aes(fill = ..level..), geom = "polygon")+
  labs(fill = "Bias", x = expression(theta[1]), y= expression(theta[2]), cex=2) + 
  theme_classic() +
  scale_fill_viridis_c(alpha = 1, begin = 0.2) +
  ylim(0,4)+
  theme(
    axis.title.x = element_text(size = 18),
    axis.text.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text = element_text(size = 18),
    legend.text=element_text(size=18),
    legend.title = element_text(size=17))


data_to_plot <- data.frame(zz=sim_bias[,2])
ggplot(data_to_plot, aes(theta_grid[17:256,1], theta_grid[17:256,2])) +
  geom_contour(aes(z = zz))
  

geom_raster(aes(fill = data_to_plot$zz))
  
  
  stat_contour() +
  labs(fill = "Bias", x = expression(theta[1]), y= expression(theta[2]), cex=2) + 
  theme_classic() +
  scale_fill_viridis_c(alpha = 1, begin = 0.2) +
  ylim(0,4) +
  theme(
    axis.title.x = element_text(size = 18),
    axis.text.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text = element_text(size = 18),
    legend.text=element_text(size=18),
    legend.title = element_text(size=17))


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


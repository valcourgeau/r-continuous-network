# loading the packages
source("~/GitHub/r-continuous-network/R/package_to_load.R")
source("~/GitHub/r-continuous-network/R/utils.R")
source("~/GitHub/r-continuous-network/R/adjacency_generation.R")

source("~/GitHub/r-continuous-network/R/levy_recovery_v2.R")
source("~/GitHub/r-continuous-network/R/utils.R")
source("~/GitHub/r-continuous-network/R/path_generation.R")
source("~/GitHub/r-continuous-network/R/grou_mle.R")

# We do not use the load data but wind since load is too correlated and not very random
AS_SPARSE <- FALSE

#######################################################################
##################### DATA PREPARATION ################################
#######################################################################

# Functions and procedures to clean the data
data_path <- "~/GitHub/r-continuous-network/data/re-europe/"
n_df_load <- 25000
n_nodes <- 50
df_load <- data.table::fread(paste(data_path, "Nodal_TS/wind_signal_COSMO.csv", sep=""), nrows = n_df_load+10)[,2:(n_nodes+1)]
df_load <- df_load[-c(1:10),]
df_load <- as.matrix(df_load)

clean_wind_data <- CleanData(df_load, frequency = 24, s.window = 24, t.window = 24*7*4)
core_wind <- clean_wind_data$remainders
plot(clean_wind_data$stl_obj$V2)
plot(core_wind[,1])

# Network topology
load_nodes <- read.csv(file=paste(data_path, "Static_data/network_nodes.csv", sep=""))
load_nodes <- load_nodes[1:n_nodes,]
topo_nodes <- data.frame("name"=load_nodes$ID, 
                         "lon"= load_nodes$longitude,
                         "lat"=load_nodes$latitude)
load_edges <- read.csv(file=paste(data_path, "Static_data/network_edges.csv", sep=""))
load_edges <- load_edges[which(load_edges$fromNode %in% 1:n_nodes & 
                                 load_edges$toNode %in% 1:n_nodes),]
topo_edges <- data.frame("from" = load_edges$fromNode, 
                         "to" = load_edges$toNode)
topo_graph <- igraph::graph.data.frame(d = topo_edges, directed = FALSE, vertices = topo_nodes)

adj_grid <- igraph::as_adjacency_matrix(topo_graph, sparse = AS_SPARSE)
adj_grid <- as.matrix(adj_grid)
mesh_size <- 2/24
observed_times <- seq(0, by=mesh_size, length.out = n_df_load)
cat('Asymptotic horizon:', mesh_size * n_df_load)
mle_theta_matrix <- GrouMLE(times=observed_times,
                            data=core_wind, adj = adj_grid, div = 1e3,
                            mode="network", output = "matrix")
mle_theta_vector <- GrouMLE(times=observed_times,
                            data=core_wind, adj = adj_grid, div = 1e3,
                            mode="network", output = "vector")
mle_theta_vector

recovery_times <- observed_times
levy_increments_recovery <- LevyRecovery(fitted_adj = mle_theta_matrix, data = core_wind, times = recovery_times, look_ahead = 1)
ghyp_levy_recovery_fit <- FitLevyRecoveryDiffusion(levy_increments_recovery$increments)
warning('Set FitLevyRecoveryDiffusion to original version')

#######################################################################
###################### LEVY FIT PLOTS #################################
#######################################################################
par(mfrow=c(2,5), mar=c(4,4.5,2,4))
cex_value <- 1.8
for(i in c(3)){ # i is the index of the plotted node
  i <- round(i)
  print(i)
  
  quantil_vals <- 1:4999/5000
  
  val_x_qq <- quantile(x = levy_increments_recovery$increments[,i], quantil_vals)
  set.seed(42)
  val_y_qq <- ghyp::rghyp(object = ghyp_levy_recovery_fit$NIG, n = length(levy_increments_recovery$increments[,i]))[, i]
  val_y_qq <- quantile(x = val_y_qq, quantil_vals)
  qqplot(x = val_x_qq, y = val_y_qq, 
         ylab = 'NIG', xlab= paste('Lévy incr. (Node ', i, ')', sep=''), 
         cex.lab = cex_value, cex.axis=cex_value)
  
  if(i == 10){title('NIG', cex.main=cex_value, line = 0.8)}
  
  set.seed(42)
  val_y_qq <- ghyp::rghyp(object = ghyp_levy_recovery_fit$GAUSS, n = length(levy_increments_recovery$increments[,i]))[, i]
  val_y_qq <- quantile(x = val_y_qq, quantil_vals)
  qqplot(x = val_x_qq, y = ghyp::rghyp(object = ghyp_levy_recovery_fit$GAUSS, n = length(levy_increments_recovery$increments[,i]))[, i], 
         ylab = 'Gaussian', xlab= paste('Lévy incr. (Node ', i, ')', sep=''), 
         cex.lab = cex_value, cex.axis=cex_value)
  if(i == 10){title('Gaussian', cex.main=cex_value, line = 0.8)}
  
  set.seed(42)
  val_y_qq <- ghyp::rghyp(object = ghyp_levy_recovery_fit$VG, n = length(levy_increments_recovery$increments[,i]))[, i]
  val_y_qq <- quantile(x = val_y_qq, quantil_vals)
  qqplot(x = val_x_qq, y = ghyp::rghyp(object = ghyp_levy_recovery_fit$VG, n = length(levy_increments_recovery$increments[,i]))[, i], 
         ylab = 'VG', xlab= paste('Lévy incr. (Node ', i, ')', sep=''), 
         cex.lab = cex_value, cex.axis=cex_value)
  if(i == 10){title('VG', cex.main=cex_value, line = 0.8)}
  
  set.seed(42)
  val_y_qq <- ghyp::rghyp(object = ghyp_levy_recovery_fit$T, n = length(levy_increments_recovery$increments[,i]))[, i]
  val_y_qq <- quantile(x = val_y_qq, quantil_vals)
  qqplot(x = val_x_qq, y = ghyp::rghyp(object = ghyp_levy_recovery_fit$T, n = length(levy_increments_recovery$increments[,i]))[, i], 
         ylab = 'Student\'s t', xlab= paste('Lévy incr. (Node ', i, ')', sep=''),
         cex.lab = cex_value, cex.axis=cex_value)
  if(i == 10){title('Student\'s t', cex.main=cex_value, line = 0.8)}
  
  set.seed(42)
  val_y_qq <- ghyp::rghyp(object = ghyp_levy_recovery_fit$FULL, n = length(levy_increments_recovery$increments[,i]))[, i]
  val_y_qq <- quantile(x = val_y_qq, quantil_vals)
  qqplot(x = val_x_qq, y = ghyp::rghyp(object = ghyp_levy_recovery_fit$FULL, n = length(levy_increments_recovery$increments[,i]))[, i], 
         ylab = 'GHYP', xlab= paste('Lévy incr. (Node ', i, ')', sep=''), 
         cex.lab = cex_value, cex.axis=cex_value)
  if(i == 10){title('GHYP', cex.main=cex_value, line = 0.8)}
}

#######################################################################
################ SIMULATION STUDY - ARTIFICIAL DATA ###################
#######################################################################

set.seed(42)
n_paths <- 20
N <- 10000
levy_increment_sims <- list()
for(i in 1:n_paths){
  levy_increment_sims[[i]] <- matrix(ghyp::rghyp(n = N, object = ghyp_levy_recovery_fit$FULL), nrow=N)
}

DO_PARALLEL <- TRUE

# choosing network type
network_types <- list(PolymerNetwork, LatticeNetwork, FullyConnectedNetwork, adj_grid)
network_types_name <- c('polymer_mles', 'lattice_mles', 'fc_mles', 're_europe_mles')

index_network <- 1
network_study <- list()

theta_1 <- mle_theta_vector[1]
theta_2 <- mle_theta_vector[2]

for(f_network in network_types){
  warning('TODO: implement file saves')
  print(index_network)
  if(index_network < 4){
    network_topo <- f_network(d = n_nodes, theta_1 = theta_1, theta_2 = theta_2) %>% as.matrix
  }else{
    network_topo <- network_types[index_network][[1]] %>% as.matrix
  }
  network_topo <- RowNormalised(network_topo) * theta_1
  diag(network_topo) <- theta_2
  
  network_topo_raw <- RowNormalised(network_topo)
  diag(network_topo_raw) <- 1.0
  
  first_point <- head(core_wind, 1)
  # first_point <- rep(0, n_nodes)
  time_init <- Sys.time() 
  if(DO_PARALLEL){
    num_cores <- detectCores()-1
    cl <- makeCluster(num_cores)
    clusterExport(cl, varlist=c("ConstructPath", "network_topo", "mesh_size", "first_point", "levy_increment_sims", "network_topo_raw",
                                "beta", "mesh_size", "n_nodes", "GrouMLE", "NodeMLELong", "CoreNodeMLE", "RowNormalised"
    ))
    generated_paths <- parLapply(
      cl,
      levy_increment_sims,
      function(x){ConstructPath(nw_topo = network_topo, noise = x, delta_time = mesh_size, first_point)}
    )
  }else{
    generated_paths <- lapply(
      X = levy_increment_sims, 
      FUN =
        function(x){
          ConstructPath(nw_topo = network_topo, noise = x, delta_time = mesh_size, first_point)
        }
    )
  }
  # generated_paths<- rlist::list.save(generated_paths, paste(network_types_name[index_network], '.RData', sep = ''))
  # generated_paths <- rlist::list.load(paste(network_types_name[index_network], '.RData', sep = ''))
  print('paths generated in')
  print(Sys.time()-time_init)
  
  # starting from topology
  #network_topo[which(abs(network_topo) > 1e-16)] <- 1
  beta <- 0.001
  n_row_generated <- nrow(generated_paths[[1]])
  
  if(DO_PARALLEL){
    clusterExport(cl, varlist=c("n_row_generated"))
    generated_fit <- parLapply(
      cl,
      generated_paths,
      function(x){
        GrouMLE(times = seq(0, length.out = n_row_generated, by = mesh_size),
                adj = as.matrix(network_topo_raw),
                data = x,
                thresholds = rep(10000, n_nodes), # mesh_size^beta
                mode = 'network',
                output = 'vector')
      }
    )
    stopCluster(cl)
  }else{
    generated_fit <- lapply(X = generated_paths,
                            FUN =
                              function(x){
                                GrouMLE(times = seq(0, length.out = nrow(generated_paths[[1]]), by = mesh_size),
                                        adj = as.matrix(network_topo_raw),
                                        data = x,
                                        thresholds = rep(mesh_size^beta, n_nodes), # mesh_size^beta
                                        mode = 'network',
                                        output = 'vector')
                              }
    )
  }
  print('fit generated')
  gen_fit_matrix <- do.call(rbind, generated_fit)
  network_study[[paste( network_types_name[index_network], '_t1', sep = '')]] <- gen_fit_matrix[,1]
  network_study[[paste(network_types_name[index_network], '_t2', sep = '')]] <- gen_fit_matrix[,2]
  index_network <- index_network + 1
}

vioplot::vioplot(network_study$polymer_mles_t2/theta_2,
                 network_study$lattice_mles_t2/theta_2,
                 network_study$fc_mles_t2/theta_2,
                 network_study$re_europe_mles_t2/theta_2)
theta_2

vioplot::vioplot(network_study$polymer_mles_t1/theta_1,
                 network_study$lattice_mles_t1/theta_1,
                 network_study$fc_mles_t1/theta_1,
                 network_study$re_europe_mles_t1/theta_1)
theta_1

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
 
#######################################################################
###################### SIMULATION STUDY - IEEE ########################
#######################################################################

pipe_layout_39 <- read.csv(file="~/GitHub/r-continuous-network/data/standard-networks/39_pipe_layout.csv")
pipe_layout_39_dt <- data.frame("from" = pipe_layout_39$Bus_i, "to" = pipe_layout_39$Bus_j)
pipe_layout_39_topo <- igraph::graph.data.frame(d = pipe_layout_39_dt, directed = FALSE)
pipe_layout_39_adj <- igraph::as_adjacency_matrix(pipe_layout_39_topo, sparse = AS_SPARSE)

pipe_layout_23 <- read.csv(file="~/GitHub/r-continuous-network/data/standard-networks/pipe_layout.csv")
pipe_layout_23_dt <- data.frame("from" = pipe_layout_23$Gnode_m, "to" = pipe_layout_23$Gnode_n)
pipe_layout_23_topo <- igraph::graph.data.frame(d = pipe_layout_23_dt, directed = FALSE)
pipe_layout_23_adj <- igraph::as_adjacency_matrix(pipe_layout_23_topo, sparse = AS_SPARSE)
pipe_layout_23_adj[pipe_layout_23_adj!= 0] <- 1 # rm those counted twice

set.seed(42)
n_paths <- 500
N <- 10000
levy_increment_sims <- list()
for(i in 1:n_paths){
  levy_increment_sims[[i]] <- matrix(ghyp::rghyp(n = N, object = ghyp_levy_recovery_fit$FULL), nrow=N)
}

# choosing network type
network_types <- list(pipe_layout_39_adj, pipe_layout_23_adj)
network_types_name <- c('pipe_layout_39_mles', 'pipe_layout_23_mles')
network_nodes <- c(39, 23)

index_network <- 1
# network_study <- list()

theta_1 <- mle_theta_vector[1]
theta_2 <- mle_theta_vector[2]

for(f_network in network_types){
  warning('TODO: implement file saves')
  print(index_network)
  n_nodes_network <- network_nodes[index_network]
  network_topo <- network_types[index_network][[1]] %>% as.matrix
  network_topo <- RowNormalised(network_topo)
  diag(network_topo) <- theta_2
  
  first_point <- head(core_wind, 1)
  first_point <- rep(0, n_nodes_network)
  time_init <- Sys.time() 
  # generated_paths <- lapply(X = levy_increment_sims, FUN =
  #                             function(x){
  #                               ConstructPath(nw_topo = network_topo, noise = x, delta_time = mesh_size, first_point)
  #                             })
  num_cores <- detectCores()-1
  cl <- makeCluster(num_cores)
  clusterExport(cl, varlist=c("ConstructPath", "network_topo", "mesh_size", "first_point", "levy_increment_sims",
                              "beta", "mesh_size", "n_nodes_network", "GrouMLE", "NodeMLELong", "CoreNodeMLE", "RowNormalised"))
  cat('levy_increment_sims', levy_increment_sims[[1]] %>% dim, '\n')
  generated_paths <- parLapply(
    cl,
    levy_increment_sims,
    function(x){ConstructPath(nw_topo = network_topo, noise = x[,1:n_nodes_network], delta_time = mesh_size, first_point)}
  )
  
  # generated_paths<- rlist::list.save(generated_paths, paste(network_types_name[index_network], '.RData', sep = ''))
  # generated_paths <- rlist::list.load(paste(network_types_name[index_network], '.RData', sep = ''))
  print('paths generated in')
  print(Sys.time()-time_init)
  
  # starting from topology
  #network_topo[which(abs(network_topo) > 1e-16)] <- 1
  beta <- 0.001
  n_row_generated <- nrow(generated_paths[[1]])
  
  clusterExport(cl, varlist=c("n_row_generated"))
  # generated_fit <- lapply(X = generated_paths,
  #   FUN =
  #     function(x){
  #       GrouMLE(times = seq(0, length.out = nrow(generated_paths[[1]]), by = mesh_size),
  #               adj = as.matrix(network_topo),
  #               data = x,
  #               thresholds = rep(mesh_size^beta, n_nodes), # mesh_size^beta
  #               mode = 'network',
  #               output = 'vector')
  #     }
  # )
  generated_fit <- parLapply(
    cl,
    generated_paths,
    function(x){
      GrouMLE(times = seq(0, length.out = n_row_generated, by = mesh_size),
              adj = as.matrix(network_topo),
              data = x,
              thresholds = rep(mesh_size^beta, n_nodes_network), # mesh_size^beta
              mode = 'network',
              output = 'vector')
    }
  )
  stopCluster(cl)
  print('fit generated')
  gen_fit_matrix <- matrix(unlist(generated_fit), ncol = 2, byrow = T)
  network_study[[paste('network_', network_types_name[index_network], '_t1', sep = '')]] <- gen_fit_matrix[,1]
  network_study[[paste('network_', network_types_name[index_network], '_t2', sep = '')]] <- gen_fit_matrix[,2]
  index_network <- index_network + 1
}


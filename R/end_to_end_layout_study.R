# loading the packages
source("~/GitHub/r-continuous-network/R/package_to_load.R")
source("~/GitHub/r-continuous-network/R/utils.R")
source("~/GitHub/r-continuous-network/R/adjacency_generation.R")

source("~/GitHub/r-continuous-network/R/levy_recovery_v2.R")
source("~/GitHub/r-continuous-network/R/utils.R")
source("~/GitHub/r-continuous-network/R/path_generation.R")
source("~/GitHub/r-continuous-network/R/grou_mle.R")

# We do not use the load data but wind since load is too correlated and not very random
AS_SPARSE <- F

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
print(adj_grid)
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
finite_levy_recovery_fit <- FitBrownianMotionCompoundPoisson(data = levy_increments_recovery$increments, mesh_size = mesh_size, thresholds = rep(mesh_size^{.499}, n_nodes))

median_n_jumps <- round(finite_levy_recovery_fit$n_jumps %>% unlist %>% median)


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
n_paths <- 100
N <- 10000
levy_increment_sims <- list()
finite_increment_sims <- list()

for(i in 1:n_paths){
  levy_increment_sims[[i]] <- matrix(ghyp::rghyp(n = N, object = ghyp_levy_recovery_fit$FULL), nrow=N)
}

set.seed(42)
for(i in 1:n_paths){
  finite_increment_sims[[i]] <- BrownianMotionCompoundPoisson(
    n = N,
    n_jumps = median_n_jumps,
    sigma = finite_levy_recovery_fit$sigma,
    jump_sigma = finite_levy_recovery_fit$jump_sigma,
    delta_time = mesh_size)
}

DO_PARALLEL <- T

# choosing network type
network_types <- list(PolymerNetwork, LatticeNetwork, FullyConnectedNetwork, adj_grid)
network_types_name <- c('polymer_mles', 'lattice_mles', 'fc_mles', 're_europe_mles')
# network_types <- list(PolymerNetwork, adj_grid)
# network_types_name <- c('polymer_mles', 're_europe_mles')

index_network <- 1
network_study <- list()
network_study_finite <- list()

theta_1 <- mle_theta_vector[1]
theta_2 <- mle_theta_vector[2]

for(f_network in network_types){
  warning('TODO: implement file saves')
  print(index_network)
  if(index_network < length(network_types_name)){ # if not adj_grid, we generate an adj matrix
    network_topo <- f_network(d = n_nodes, theta_1 = theta_1, theta_2 = theta_2) %>% as.matrix
  }else{
    network_topo <- network_types[index_network][[1]] %>% as.matrix
  }
  network_topo <- RowNormalised(network_topo) * theta_1
  diag(network_topo) <- theta_2
  
  # print('network_topo')
  # print(network_topo)
  
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
  if(DO_PARALLEL){
    num_cores <- detectCores()-1
    cl <- makeCluster(num_cores)
    clusterExport(cl, varlist=c("ConstructPath", "network_topo", "mesh_size", "first_point", "levy_increment_sims", "network_topo_raw",
                                "beta", "mesh_size", "n_nodes", "GrouMLE", "NodeMLELong", "CoreNodeMLE", "RowNormalised"
    ))
    generated_paths_finite <- parLapply(
      cl,
      finite_increment_sims,
      function(x){ConstructPath(nw_topo = network_topo, noise = x, delta_time = mesh_size, first_point)}
    )
  }else{
    generated_paths_finite <- lapply(
      X = finite_increment_sims, 
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
  beta <- 0.4999
  n_row_generated <- nrow(generated_paths[[1]])
  time_init <- Sys.time() 
  if(DO_PARALLEL){
    clusterExport(cl, varlist=c("n_row_generated"))
    generated_fit <- parLapply(
      cl,
      generated_paths,
      function(x){
        GrouMLE(times = seq(0, length.out = n_row_generated, by = mesh_size),
                adj = as.matrix(network_topo_raw),
                data = x,
                thresholds = rep(mesh_size^beta, n_nodes), # mesh_size^beta
                mode = 'network',
                output = 'vector')
      }
    )
  }else{
    generated_fit <- lapply(
      X = generated_paths,
      FUN =
        function(x){
          vect <- GrouMLE(times = seq(0, length.out = nrow(generated_paths[[1]]), by = mesh_size),
                         adj = as.matrix(network_topo_raw),
                         data = x,
                         thresholds = rep(mesh_size^beta, n_nodes), # mesh_size^beta
                         mode = 'network',
                         output = 'vector')
         return(vect)
        }
    )
  }

    # finite
  if(DO_PARALLEL){
    clusterExport(cl, varlist=c("n_row_generated"))
    generated_fit_finite <- parLapply(
      cl,
      generated_paths_finite,
      function(x){
        vect <- GrouMLE(times = seq(0, length.out = n_row_generated, by = mesh_size),
                adj = as.matrix(network_topo_raw),
                data = x,
                thresholds = rep(mesh_size^beta, n_nodes), # mesh_size^beta
                mode = 'network',
                output = 'vector')
        return(vect)
      }
    )
    stopCluster(cl)
  }else{
    generated_fit_finite <- lapply(
      X = generated_paths_finite,
      FUN =
        function(x){
          vect <- GrouMLE(times = seq(0, length.out = nrow(generated_paths_finite[[1]]), by = mesh_size),
                          adj = as.matrix(network_topo_raw),
                          data = x,
                          thresholds = rep(mesh_size^beta, n_nodes), # mesh_size^beta
                          mode = 'network',
                          output = 'vector')
          return(vect)
        }
    )
  }
  print('fit generated')
  print(Sys.time()-time_init)
  gen_fit_matrix <- do.call(rbind, generated_fit)
  network_study[[paste(network_types_name[index_network], '_t1', sep = '')]] <- gen_fit_matrix[,1]
  network_study[[paste(network_types_name[index_network], '_t2', sep = '')]] <- gen_fit_matrix[,2]
  
  gen_fit_matrix_finite <- do.call(rbind, generated_fit_finite)
  network_study_finite[[paste(network_types_name[index_network], '_t1', sep = '')]] <- gen_fit_matrix_finite[,1]
  network_study_finite[[paste(network_types_name[index_network], '_t2', sep = '')]] <- gen_fit_matrix_finite[,2]
  
  index_network <- index_network + 1
}

par(mfrow=c(1,2))
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
layout_names <- c('RE-Europe 50', 'Polymer', 'Lattice', 'Complete')
par(mfrow=c(2,1), mar=c(2,3,1.8,0.5))
vioplot::vioplot(network_study[c(7,1,3,5)], names = c('RE-Europe 50', 'Polymer', 'Lattice', 'Complete'),
                 col = colors,
                 # col = c('#BE2B2B', '#2BBEBE', '#2BBEBE', '#2BBEBE'),
                 main=expression(theta[1]),
                 cex.main=1.8, cex.names=1.5, cex.axis=1.5)
par(xpd = F) 
abline(h=theta_1, lwd=2, lty=2)
vioplot::vioplot(network_study[c(8,2,4,6)], names = c('RE-Europe 50', 'Polymer', 'Lattice', 'Complete'),
                 col =colors,
                 # col = c('#BE2B2B', '#2BBEBE', '#2BBEBE', '#2BBEBE'),
                 main=expression(theta[2]),
                 cex.main=1.8, cex.names=1.5, cex.axis=1.5, ylim=theta_2*c(0.97,1.03))
abline(h=theta_2, lwd=2, lty=2)

dt_layout <- data.frame(
  layout=rep(layout_names, each=n_paths),
  theta_1 = c(network_study$re_europe_mles_t1, network_study$polymer_mles_t1, network_study$lattice_mles_t1, network_study$fc_mles_t1),
  theta_2 = c(network_study$re_europe_mles_t2, network_study$polymer_mles_t2, network_study$lattice_mles_t2, network_study$fc_mles_t2),
  theta_1_finite = c(network_study_finite$re_europe_mles_t1, network_study_finite$polymer_mles_t1, network_study_finite$lattice_mles_t1, network_study_finite$fc_mles_t1),
  theta_2_finite = c(network_study_finite$re_europe_mles_t2, network_study_finite$polymer_mles_t2, network_study_finite$lattice_mles_t2, network_study_finite$fc_mles_t2)
)

factor_topo <- factor(dt_layout$layout, level = c('RE-Europe 50', 'Polymer', 'Lattice', 'Complete'))
multiplicator <- .4
ggplot_theta_1_finite <- ggplot2::ggplot(data=dt_layout, ggplot2::aes(x=factor_topo, y=theta_1_finite, fill=layout)) +
  ggplot2::geom_violin() +
  ggplot2::ggtitle(' ') +
  ggplot2::geom_boxplot(width=0.1, fill='white') +
  ggplot2::geom_hline(yintercept=theta_1, linetype="dashed", color = "black", size=1.1) +
  ggplot2::scale_fill_manual(values=colors[4:1]) + 
  ggplot2::ylab(expression(paste(theta[1]))) +
  ggplot2::xlab('') +
  ggplot2::theme(legend.position = "none", 
                 panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                 panel.background = ggplot2::element_blank(), 
                 strip.background = ggplot2::element_blank(),
                 axis.title.y = ggplot2::element_text(size = ggplot2::rel(multiplicator*1.5)),
                 axis.title.x = ggplot2::element_text(size = ggplot2::rel(multiplicator*1.2)),
                 strip.text.x = ggplot2::element_text(size = ggplot2::rel(multiplicator*.8), color = "black"),
                 # strip.text.x = ggplot2::element_text(size = ggplot2::rel(1.3), color = "black"),
                 axis.text.x = ggplot2::element_text(size = ggplot2::rel(multiplicator*1.4), color = "black"),
                 # legend.title = ggplot2::element_text(size = ggplot2::rel(multiplicator*1.2), color = "black"),
                 # legend.text = ggplot2::element_text(size = ggplot2::rel(multiplicator*1.2), color = "black"),
                 strip.text =  ggplot2::element_text(size = ggplot2::rel(multiplicator*1.0), color = "black"),
                 axis.text.y = ggplot2::element_text(size = ggplot2::rel(multiplicator*1.5), color = "black"),
                 plot.title = ggplot2::element_text(size = ggplot2::rel(multiplicator*1.5), color = "black", hjust =-.8))
ggplot_theta_2_finite <- ggplot2::ggplot(data=dt_layout, ggplot2::aes(x=factor_topo, y=theta_2_finite, fill=layout)) +
  ggplot2::ggtitle('Finite jump activity') +
  ggplot2::geom_violin() +
  ggplot2::geom_boxplot(width=0.1, fill='white') +
  ggplot2::geom_hline(yintercept=theta_2, linetype="dashed", color = "black", size=1.1) +
  ggplot2::scale_fill_manual(values=colors[4:1]) + 
  ggplot2::ylab(expression(paste(theta[2]))) +
  ggplot2::xlab('') +
  ggplot2::ylim(c(0.98, 1.02) * theta_2) +
  ggplot2::theme(legend.position = "none", 
                 panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                 panel.background = ggplot2::element_blank(), 
                 strip.background = ggplot2::element_blank(),
                 axis.title.y = ggplot2::element_text(size = ggplot2::rel(multiplicator*1.5)),
                 axis.title.x = ggplot2::element_text(size = ggplot2::rel(multiplicator*1.2)),
                 strip.text.x = ggplot2::element_text(size = ggplot2::rel(multiplicator*.8), color = "black"),
                 # strip.text.x = ggplot2::element_text(size = ggplot2::rel(1.3), color = "black"),
                 axis.text.x = ggplot2::element_text(size = ggplot2::rel(multiplicator*1.4), color = "black"),
                 # legend.title = ggplot2::element_text(size = ggplot2::rel(multiplicator*1.2), color = "black"),
                 # legend.text = ggplot2::element_text(size = ggplot2::rel(multiplicator*1.2), color = "black"),
                 strip.text =  ggplot2::element_text(size = ggplot2::rel(multiplicator*1.0), color = "black"),
                 axis.text.y = ggplot2::element_text(size = ggplot2::rel(multiplicator*1.5), color = "black"),
                 plot.title = ggplot2::element_text(size = ggplot2::rel(multiplicator*1.5), color = "black", hjust = -.7))
ggplot_theta_1 <- ggplot2::ggplot(data=dt_layout, ggplot2::aes(x=factor_topo, y=theta_1, fill=layout)) +
  ggplot2::geom_violin() +
  ggplot2::ggtitle(' ') +
  ggplot2::geom_boxplot(width=0.1, fill='white') +
  ggplot2::geom_hline(yintercept=theta_1, linetype="dashed", color = "black", size=1.1) +
  ggplot2::scale_fill_manual(values=colors[4:1]) + 
  ggplot2::ylab(expression(paste(theta[1]))) +
  ggplot2::xlab('') +
  ggplot2::theme(legend.position = "none", 
                 panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                 panel.background = ggplot2::element_blank(), 
                 strip.background = ggplot2::element_blank(),
                 axis.title.y = ggplot2::element_text(size = ggplot2::rel(multiplicator*1.5)),
                 axis.title.x = ggplot2::element_text(size = ggplot2::rel(multiplicator*1.2)),
                 strip.text.x = ggplot2::element_text(size = ggplot2::rel(multiplicator*.8), color = "black"),
                 # strip.text.x = ggplot2::element_text(size = ggplot2::rel(1.3), color = "black"),
                 axis.text.x = ggplot2::element_text(size = ggplot2::rel(multiplicator*1.4), color = "black"),
                 # legend.title = ggplot2::element_text(size = ggplot2::rel(multiplicator*1.2), color = "black"),
                 # legend.text = ggplot2::element_text(size = ggplot2::rel(multiplicator*1.2), color = "black"),
                 strip.text =  ggplot2::element_text(size = ggplot2::rel(multiplicator*1.0), color = "black"),
                 axis.text.y = ggplot2::element_text(size = ggplot2::rel(multiplicator*1.5), color = "black"),
                 plot.title = ggplot2::element_text(size = ggplot2::rel(multiplicator*1.5), color = "black", hjust = -.7))
ggplot_theta_2 <- ggplot2::ggplot(data=dt_layout, ggplot2::aes(x=factor_topo, y=theta_2, fill=layout)) +
  ggplot2::ggtitle('Infinite jump activity') +
  ggplot2::geom_violin() +
  ggplot2::geom_boxplot(width=0.1, fill='white') +
  ggplot2::geom_hline(yintercept=theta_2, linetype="dashed", color = "black", size=1.1) +
  ggplot2::scale_fill_manual(values=colors[4:1]) + 
  ggplot2::ylab(expression(paste(theta[2]))) +
  ggplot2::xlab('') +
  ggplot2::ylim(c(0.98, 1.02) * theta_2) +
  ggplot2::theme(legend.position = "none", 
                 panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                 panel.background = ggplot2::element_blank(), 
                 strip.background = ggplot2::element_blank(),
                 axis.title.y = ggplot2::element_text(size = ggplot2::rel(multiplicator*1.5)),
                 axis.title.x = ggplot2::element_text(size = ggplot2::rel(multiplicator*1.2)),
                 strip.text.x = ggplot2::element_text(size = ggplot2::rel(multiplicator*.8), color = "black"),
                 # strip.text.x = ggplot2::element_text(size = ggplot2::rel(1.3), color = "black"),
                 axis.text.x = ggplot2::element_text(size = ggplot2::rel(multiplicator*1.4), color = "black"),
                 # legend.title = ggplot2::element_text(size = ggplot2::rel(multiplicator*1.2), color = "black"),
                 # legend.text = ggplot2::element_text(size = ggplot2::rel(multiplicator*1.2), color = "black"),
                 strip.text =  ggplot2::element_text(size = ggplot2::rel(multiplicator*1.0), color = "black"),
                 axis.text.y = ggplot2::element_text(size = ggplot2::rel(multiplicator*1.5), color = "black"),
                 plot.title = ggplot2::element_text(size = ggplot2::rel(multiplicator*1.5), color = "black", hjust = -.8))

setEPS()
postscript("../data/pictures/layout.eps")
# png("../data/pictures/layout.png", width = 1600, height = 800)
gridExtra::grid.arrange(
  ggplot_theta_1_finite, ggplot_theta_2_finite,
  ggplot_theta_1, ggplot_theta_2,
  ncol = 2, nrow = 2)
dev.off()
 
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


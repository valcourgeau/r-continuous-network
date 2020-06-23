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

clean_wind_data <- CleanData(df_load, frequency = 24, s.window = 24, t.window = 24)
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


set.seed(42)
n_paths <- 2
N <- 10000
levy_increment_sims <- list()
for(i in 1:n_paths){
  levy_increment_sims[[i]]<- matrix(ghyp::rghyp(n = N, object = ghyp_levy_recovery_fit$FULL), nrow=N)
}

DO_PARALLEL <- F

# choosing network type
network_types <- list(PolymerNetwork, LatticeNetwork, FullyConnectedNetwork, adj_grid)
network_types_name <- c('polymer_mles', 'lattice_mles', 'fc_mles', 're_europe_mles')

index_network <- 1
horizon_study <- list()

theta_1 <- mle_theta_vector[1]
theta_2 <- mle_theta_vector[2]

mesh_division <- 16
mesh_min <- 1e-3
mesh_max <- 0.5

for(f_network in network_types){
  warning('TODO: implement file saves')
  print(index_network)
  if(index_network < 4){
    network_topo <- f_network(d = n_nodes, theta_1 = theta_1, theta_2 = theta_2)
  }else{
    network_topo <- network_types[index_network][[1]] %>% as.matrix
  }
  network_topo <- RowNormalised(network_topo) * theta_1
  diag(network_topo) <- theta_2
  network_topo_raw <- RowNormalised(network_topo)
  diag(network_topo_raw) <- 1.0
  
  first_point <- head(core_wind, 1)
  time_init <- Sys.time() 
  
  beta <- 0.001
  mesh_sequence <- seq(from=log(mesh_min), to=log(mesh_max), length.out = mesh_division)
  
  if(DO_PARALLEL){
    num_cores <- detectCores()-1
    cl <- makeCluster(num_cores)
    clusterExport(cl, varlist=c("ConstructPath", "network_topo", "mesh_size", "first_point", "levy_increment_sims",
                                "network_topo_raw", "mesh_size", "n_nodes", "GrouMLE", "mesh_sequence",
                                "NodeMLELong", "CoreNodeMLE", "RowNormalised"))
    generated_paths <- parLapply(
      cl,
      levy_increment_sims,
      function(x){
        lapply(
          mesh_sequence,
          function(mesh){
            list(
              path=ConstructPath(nw_topo = network_topo, noise = x, delta_time = exp(mesh), first_point),
              mesh=mesh
            )
          })
      }
    )
  }else{
    generated_paths <- lapply(
      X = levy_increment_sims,
      FUN =
        function(x){
          lapply(
            mesh_sequence,
            function(mesh){
              list(
                path=ConstructPath(nw_topo = network_topo, noise = x, delta_time = exp(mesh), first_point),
                mesh=mesh
              )
            })
        }
    )
  }
  # generated_paths<- rlist::list.save(generated_paths, paste(network_types_name[index_network], '.RData', sep = ''))
  # generated_paths <- rlist::list.load(paste(network_types_name[index_network], '.RData', sep = ''))
  print('paths generated in')
  print(Sys.time()-time_init)
  
  # starting from topology
  #network_topo[which(abs(network_topo) > 1e-16)] <- 1
  n_row_generated <- N #nrow(generated_paths[[1]][[1]])
  
  time_init <- Sys.time()
  if(DO_PARALLEL){
    clusterExport(cl, varlist=c("n_row_generated"))
    generated_fit <- parLapply(
      cl,
      generated_paths,
      fun = function(x_with_differen_mesh){
        lapply(
          x_with_differen_mesh,
          function(x){
            gen_fit[['mle']] <- GrouMLE(times = seq(0, length.out = n_row_generated, by = exp(x$mesh)),
                                        adj = as.matrix(network_topo_raw), # TODO network_topo_raw???
                                        data = x$path,
                                        thresholds = rep(exp(x$mesh)^beta, n_nodes), # mesh_size^beta
                                        mode = 'node',
                                        output = 'matrix')
            gen_fit[['fasen']] <- FasenRegression(x$path)
            gen_fit[['truth']] <- network_topo
            gen_fit[['mesh']] <- x$mesh
            return(gen_fit)
          }
        )
      }
    )
    stopCluster(cl)
  }else{
    generated_fit <- lapply(
      X = generated_paths,
      FUN =
        function(x_with_differen_mesh){
          lapply(
            x_with_differen_mesh,
            function(x){
              gen_fit <- list()
              # print('path')
              # print(x$path)
              # print('mesh')
              # print(x$mesh)
              gen_fit[['mle']] <- GrouMLE(times = seq(0, length.out = n_row_generated, by = exp(x$mesh)),
                                          adj = as.matrix(network_topo_raw),
                                          data = as.matrix(x$path),
                                          thresholds = rep(exp(x$mesh)^beta, n_nodes), # mesh_size^beta
                                          mode = 'node',
                                          output = 'matrix')
              gen_fit[['fasen']] <- FasenRegression(x$path)
              gen_fit[['truth']] <- network_topo
              gen_fit[['mesh']] <- x$mesh
              return(gen_fit)
            }
          )
        }
    )
  }
  print('fit generated')
  print(Sys.time()-time_init)
  horizon_study[[index_network]] <- generated_fit
  index_network <- index_network + 1
}

horizon_array_mle <- array(dim=c(4, n_paths, mesh_division-2))
horizon_array_ls <- array(dim=c(4, n_paths, mesh_division-2))
expm_truth <- as.matrix(Matrix::expm(-horizon_study[[1]][[1]][[1]]$truth*mesh_size))
for(i in 1:4){
  for(j in 1:n_paths){
    for(k in 1:(mesh_division-2)){
      mle_val <- horizon_study[[i]][[j]][[k]]$mle
      mesh_val <- exp(horizon_study[[i]][[j]][[k]]$mesh)
      expm_mle <- as.matrix(Matrix::expm(-mle_val*mesh_val))
      expm_ls <- horizon_study[[i]][[j]][[k]]$fasen

      horizon_array_mle[i,j,k] <- sqrt(sum((expm_mle-expm_truth)^2))/sqrt(sum(expm_truth^2))
      horizon_array_ls[i,j,k] <- sqrt(sum((expm_ls-expm_truth)^2))/sqrt(sum(expm_truth^2))
    }
  }
}

colors <- c('#08605F', '#598381', '#8E936D', '#A2AD59')

par(mfrow=c(1,1), mar = c(5,5,2.5,1))
plot_unique_names <- c('RE-Europe 50', 'Polymer', 'Lattice', 'Complete')

plot_names <- as.vector(t(vapply(
  c('MLE', 'LS'),
  function(type_of_estimator){
    vapply(plot_unique_names, function(x){paste(x, type_of_estimator)}, 'w')
  }, rep('w', length(plot_unique_names)))))
plot_names

cut_tail <- -2

# setEPS()
# postscript("../data/pictures/horizon.eps")
png("../data/pictures/horizon.png", width = 800, height = 300)
x_matplot <- matrix(rep(exp(head(mesh_sequence, cut_tail)), 4), nrow=4, byrow = T)
matplot(x=t(x_matplot), t(apply(horizon_array_mle, c(1,3), mean))[,c(4,1:3)],
        type = 'b', ylim=c(0.025, 1), log='xy', pch=rep(22:25, each=1),   bty = "n",
        lwd=rep(2, 4),  lty=rep(1, 4), col=colors, cex.axis=1.4, cex.lab=1.5,
        xlab=expression(paste('Mesh size ', Delta[N], ' (log)')),
        ylab='Relative Error Norm', main='Relative L2-Norm Error vs Mesh size with N=10000')
legend(.001, .39, plot_names, lty=rep(c(1, 2), 4), bty = "n",
       pch=rep(22:25, each=2), lwd=rep(2, 8), cex=1.1, ncol = 1,
       col = rep(colors, each=2))
matplot(x=t(x_matplot), t(apply(horizon_array_ls, c(1,3), mean))[,c(4,1:3)],
        pch=23, type='b', add=T, lty=rep(2, 4), lwd=rep(2, 4), col=colors)
abline(v=2/24, lty=2, lwd=2)
dev.off()
# abline(v=0.095333310, lwd=2)


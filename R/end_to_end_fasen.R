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
finite_levy_recovery_fit <- FitBrownianMotionCompoundPoisson(data = levy_increments_recovery$increments, mesh_size = mesh_size, thresholds = rep(mesh_size^{.499}, n_nodes))

median_n_jumps <- round(finite_levy_recovery_fit$n_jumps %>% unlist %>% median)

set.seed(42)
n_paths <- 5
N <- 10000
levy_increment_sims <- list()
finite_increment_sims <- list()

for(i in 1:n_paths){
  levy_increment_sims[[i]]<- matrix(ghyp::rghyp(n = N, object = ghyp_levy_recovery_fit$FULL), nrow=N)
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

DO_PARALLEL <- F

# choosing network type
network_types <- list(PolymerNetwork, LatticeNetwork, FullyConnectedNetwork, adj_grid)
network_types_name <- c('polymer_mles', 'lattice_mles', 'fc_mles', 're_europe_mles')

index_network <- 1
fasen_study <- list()
fasen_study_finite <- list()

theta_1 <- mle_theta_vector[1]
theta_2 <- mle_theta_vector[2]

sig_division <- 10
sig_min <- 0.5
sig_max <- 1e4

sig_sequence <- seq(from=log(sig_min), to=log(sig_max), length.out = sig_division)
sig_sequence <- log(c(0.5, 1, 10, 100, 1000))

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
  
  beta <- 0.499
  
  if(DO_PARALLEL){
    num_cores <- detectCores()-1
    cl <- makeCluster(num_cores)
    clusterExport(cl, varlist=c("ConstructPath", "network_topo", "mesh_size", "first_point", "levy_increment_sims",
                                "network_topo_raw", "mesh_size", "n_nodes", "GrouMLE", "sig_sequence",
                                "NodeMLELong", "CoreNodeMLE", "RowNormalised"))
    generated_paths <- parLapply(
      cl,
      levy_increment_sims,
      function(x){
        lapply(
          sig_sequence,
          function(sig){
            ConstructPath(nw_topo = network_topo, noise = exp(sig)*x, delta_time = mesh_size, first_point)
          })
      }
    )
    generated_paths_finite <- parLapply(
      cl,
      finite_increment_sims,
      function(x){
        lapply(
          sig_sequence,
          function(sig){
            ConstructPath(nw_topo = network_topo, noise = exp(sig)*x, delta_time = mesh_size, first_point)
          })
      }
    )
  }else{
    generated_paths <- lapply(
      X = levy_increment_sims,
      FUN =
        function(x){
          lapply(
            sig_sequence,
            function(sig){
              print(sig)
              ConstructPath(nw_topo = network_topo, noise = exp(sig)*x, delta_time = mesh_size, first_point)
            })
        }
    )
    generated_paths_finite <- lapply(
      X = finite_increment_sims,
      FUN =
        function(x){
          lapply(
            sig_sequence,
            function(sig){
              print(sig)
              ConstructPath(nw_topo = network_topo, noise = exp(sig)*x, delta_time = mesh_size, first_point)
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
  n_row_generated <- nrow(generated_paths[[1]][[1]])
  
  time_init <- Sys.time()
  if(DO_PARALLEL){
    clusterExport(cl, varlist=c("n_row_generated"))
    generated_fit <- parLapply(
      cl,
      generated_paths,
      fun = function(x_with_differen_sigmas){
        lapply(
          x_with_differen_sigmas,
          function(x){
            gen_fit <- list()
            gen_fit[['mle']] <- GrouMLE(times = seq(0, length.out = n_row_generated, by = mesh_size),
                                        adj = as.matrix(network_topo_raw), # TODO network_topo_raw???
                                        data = x,
                                        thresholds = rep(mesh_size^beta, n_nodes), # mesh_size^beta
                                        mode = 'node',
                                        output = 'matrix')
            gen_fit[['fasen']] <- FasenRegression(x)
            gen_fit[['truth']] <- network_topo
            
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
        function(x_with_differen_sigmas){
          lapply(
            x_with_differen_sigmas,
            function(x){
              gen_fit <- list()
              gen_fit[['mle']] <- GrouMLE(times = seq(0, length.out = n_row_generated, by = mesh_size),
                                          adj = as.matrix(network_topo_raw),
                                          data = x,
                                          thresholds = rep(mesh_size^beta, n_nodes), # mesh_size^beta
                                          mode = 'node',
                                          output = 'matrix')
              # print('mle')
              # print(gen_fit[['mle']])
              # print('network_topo')
              # print(network_topo)
              gen_fit[['fasen']] <- FasenRegression(x)
              gen_fit[['truth']] <- network_topo
              return(gen_fit)
            }
          )
        }
    )
  }
  
  # FINITE
  if(DO_PARALLEL){
    clusterExport(cl, varlist=c("n_row_generated"))
    generated_fit_finite <- parLapply(
      cl,
      generated_paths_finite,
      fun = function(x_with_differen_sigmas){
        lapply(
          x_with_differen_sigmas,
          function(x){
            gen_fit <- list()
            gen_fit[['mle']] <- GrouMLE(times = seq(0, length.out = n_row_generated, by = mesh_size),
                                        adj = as.matrix(network_topo_raw), # TODO network_topo_raw???
                                        data = x,
                                        thresholds = rep(mesh_size^beta, n_nodes), # mesh_size^beta
                                        mode = 'node',
                                        output = 'matrix')
            gen_fit[['fasen']] <- FasenRegression(x)
            gen_fit[['truth']] <- network_topo
            
            return(gen_fit)
          }
        )
      }
    )
    stopCluster(cl)
  }else{
    generated_fit_finite <- lapply(
      X = generated_paths_finite,
      FUN =
        function(x_with_differen_sigmas){
          lapply(
            x_with_differen_sigmas,
            function(x){
              gen_fit <- list()
              gen_fit[['mle']] <- GrouMLE(times = seq(0, length.out = n_row_generated, by = mesh_size),
                                          adj = as.matrix(network_topo_raw),
                                          data = x,
                                          thresholds = rep(mesh_size^beta, n_nodes), # mesh_size^beta
                                          mode = 'node',
                                          output = 'matrix')
              # print('mle')
              # print(gen_fit[['mle']])
              # print('network_topo')
              # print(network_topo)
              gen_fit[['fasen']] <- FasenRegression(x)
              gen_fit[['truth']] <- network_topo
              return(gen_fit)
            }
          )
        }
    )
    print('generated_fit_finite')
  }
  print('fit generated')
  print(Sys.time()-time_init)
  fasen_study[[index_network]] <- generated_fit
  fasen_study_finite[[index_network]] <- generated_fit_finite
  index_network <- index_network + 1
}

fasen_array_mle <- array(dim=c(4, n_paths, length(sig_sequence)))
fasen_array_ls <- array(dim=c(4, n_paths, length(sig_sequence)))
expm_truth <- as.matrix(Matrix::expm(-fasen_study[[1]][[1]][[1]]$truth*mesh_size))
for(i in 1:4){
  for(j in 1:n_paths){
    for(k in 1:length(sig_sequence)){
      expm_mle <- as.matrix(Matrix::expm(-fasen_study[[i]][[j]][[k]]$mle*mesh_size))
      expm_ls <- fasen_study[[i]][[j]][[k]]$fasen
      fasen_array_mle[i,j,k] <- sqrt(sum((expm_mle-expm_truth)^2))/sqrt(sum(expm_truth^2))
      fasen_array_ls[i,j,k] <- sqrt(sum((expm_ls-expm_truth)^2))/sqrt(sum(expm_truth^2))
    }
  }
}

fasen_array_finite_mle <- array(dim=c(4, n_paths, length(sig_sequence)))
fasen_array_finite_ls <- array(dim=c(4, n_paths, length(sig_sequence)))
expm_truth <- as.matrix(Matrix::expm(-fasen_study_finite[[1]][[1]][[1]]$truth*mesh_size))
for(i in 1:4){
  for(j in 1:n_paths){
    for(k in 1:length(sig_sequence)){
      expm_mle <- as.matrix(Matrix::expm(-fasen_study_finite[[i]][[j]][[k]]$mle*mesh_size))
      expm_ls <- fasen_study_finite[[i]][[j]][[k]]$fasen
      fasen_array_finite_mle[i,j,k] <- sqrt(sum((expm_mle-expm_truth)^2))/sqrt(sum(expm_truth^2))
      fasen_array_finite_ls[i,j,k] <- sqrt(sum((expm_ls-expm_truth)^2))/sqrt(sum(expm_truth^2))
    }
  }
}

colors <- c('#08605F', '#598381', '#8E936D', '#A2AD59')
plot_unique_names <- c('RE-Europe 50', 'Polymer', 'Lattice', 'Complete')

# diverse plots
{
  par(mfrow=c(1,2), mar = c(5,5,2.5,1))
  plot_names <- as.vector(t(vapply(
    c('MLE', 'LS'),
    function(type_of_estimator){
      vapply(plot_unique_names, function(x){paste(x, type_of_estimator)}, 'w')
    }, rep('w', length(plot_unique_names)))))
  plot_names
  x_matplot <- matrix(rep(exp(sig_sequence), 4), nrow=4, byrow = T)
  matplot(x=t(x_matplot), t(apply(fasen_array_mle, c(1,3), mean))[,c(4,1:3)],
          type = 'b', ylim=c(0.005, 0.5), log='xy', pch=rep(22:25, each=1),  bty = "n",
          lwd=rep(2, 4),  lty=rep(1, 4), col=colors, cex.axis=1.4, cex.lab=1.5,
          xlab=expression(paste('Noise multiplicator ', sigma, ' (log)')), ylab='Relative Error Norm')
  legend(.5, 0.7, plot_names, lty=rep(c(1, 2), 4), bty = "n",
         pch=rep(22:25, each=2), lwd=rep(2, 8), cex=1.1, ncol = 2,
         col = rep(colors, each=2))
  matplot(x=t(x_matplot), t(apply(fasen_array_ls, c(1,3), mean))[,c(4,1:3)],
          pch=rep(22:25, each=1), type='b', add=T, lty=rep(2, 4), lwd=rep(2, 4), col=colors)
  # legend(1, 1.5, network_types_name,
  #        pch=22,
  #        col = 1:4)
  
  # subplot with zoom
  index_sub_plot <- x_matplot[1,] <= 150
  x_sub_matplot <- matrix(rep(x_matplot[1,][index_sub_plot], ncol(x_matplot)), nrow = ncol(x_matplot), byrow = T)
  matplot(x=t(x_sub_matplot), t(apply(fasen_array_mle, c(1,3), mean))[index_sub_plot,c(4,1:3)],
          type = 'b', ylim=c(0.005, .5), log='xy', pch=22:25,  bty = "n",
          lwd=rep(2, 4),  lty=rep(1, 4), col=colors, cex.axis=1.4, cex.lab=1.5,
          xlab=expression(paste('Noise multiplicator ', sigma, ' (log)')), ylab='Relative Error Norm (log)')
  legend(.5, .7, plot_names[(1:4-1)*2+1], lty=rep(c(1, 2), 4)[(1:4-1)*2+1],  bty = "n",
         pch=22:25, lwd=rep(2, 8)[(1:4-1)*2+1], cex=1.1,
         col = rep(colors, each=2)[(1:4-1)*2+1])
}

setEPS()
postscript("../data/pictures/fasen.eps")
# png("../data/pictures/fasen.png", width = 1600, height = 800)
index_bar_plot <- 1:length(sig_sequence)
bar_plot_mle <- t(apply(fasen_array_finite_mle, c(1,3), mean))[index_bar_plot, c(4,1:3)]
bar_plot_ls <- t(apply(fasen_array_finite_ls, c(1,3), mean))[index_bar_plot, c(4,1:3)]

bar_plot_mle_sd <- t(apply(fasen_array_finite_mle, c(1,3), sd))[index_bar_plot, c(4,1:3)]
bar_plot_ls_sd <- t(apply(fasen_array_finite_ls, c(1,3), sd))[index_bar_plot, c(4,1:3)]

dt_fasen <- data.frame(type=c(rep('LS', ncol(bar_plot_ls)*length(index_bar_plot)), rep('MLE', ncol(bar_plot_ls)*length(index_bar_plot))),
                       topo=rep(rep(plot_unique_names, each=length(index_bar_plot)), 2),
                       sigma=as.factor(rep(rep(exp(sig_sequence[index_bar_plot]), ncol(bar_plot_ls)), 2)),
                       ren=c(c(bar_plot_ls), c(bar_plot_mle)),
                       ren_sd=c(c(bar_plot_ls_sd, bar_plot_mle_sd)))

# reorder
multiplicator <- 2
factor_topo <- factor(dt_fasen$topo, level = c('RE-Europe 50', 'Polymer', 'Lattice', 'Complete'))
ggplot_finite <- ggplot2::ggplot(data=dt_fasen, ggplot2::aes(x=sigma, y=ren, fill=type)) +
  ggplot2::ggtitle('Finite jump activity') +
  ggplot2::facet_wrap(factor_topo, scales='free_x', ncol = 4) +
  ggplot2::geom_bar(stat="identity", position=ggplot2::position_dodge()) +
  ggplot2::scale_fill_manual(values=colors[c(3,1)]) + 
  ggplot2::geom_errorbar(ggplot2::aes(ymin=ren-1.96*ren_sd, ymax=ren+1.96*ren_sd), width=.5, size=1.2,
                         position=ggplot2::position_dodge(.9)) +
  ggplot2::ylab(expression(paste('Relative  ', L^2, '-Error Norm'))) +
  ggplot2::xlab(expression(paste('Noise multiplicator ', sigma))) +
  ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                 panel.background = ggplot2::element_blank(), 
                 strip.background = ggplot2::element_blank(),
                 axis.title.y = ggplot2::element_text(size = ggplot2::rel(multiplicator*1.3)),
                 axis.title.x = ggplot2::element_text(size = ggplot2::rel(multiplicator*1.2)),
                 strip.text.x = ggplot2::element_text(size = ggplot2::rel(multiplicator*.8), color = "black"),
                 # strip.text.x = ggplot2::element_text(size = ggplot2::rel(1.3), color = "black"),
                 axis.text.x = ggplot2::element_text(size = ggplot2::rel(multiplicator*1.4), color = "black"),
                 legend.title = ggplot2::element_text(size = ggplot2::rel(multiplicator*1.2), color = "black"),
                 legend.text = ggplot2::element_text(size = ggplot2::rel(multiplicator*1.2), color = "black"),
                 strip.text =  ggplot2::element_text(size = ggplot2::rel(multiplicator*1.0), color = "black"),
                 axis.text.y = ggplot2::element_text(size = ggplot2::rel(multiplicator*1.5), color = "black"),
                 plot.title = ggplot2::element_text(size = ggplot2::rel(multiplicator*1.5), color = "black", hjust = 0.5)) +
  ggplot2::labs(fill="Estimator")



par(mfrow=c(1,2), mar = c(5,5,2.5,1))
index_bar_plot <- 1:length(sig_sequence)
bar_plot_mle <- t(apply(fasen_array_mle, c(1,3), mean))[index_bar_plot, c(4,1:3)]
bar_plot_ls <- t(apply(fasen_array_ls, c(1,3), mean))[index_bar_plot, c(4,1:3)]

bar_plot_mle_sd <- t(apply(fasen_array_mle, c(1,3), sd))[index_bar_plot, c(4,1:3)]
bar_plot_ls_sd <- t(apply(fasen_array_ls, c(1,3), sd))[index_bar_plot, c(4,1:3)]

dt_fasen <- data.frame(type=c(rep('LS', ncol(bar_plot_ls)*length(index_bar_plot)), rep('MLE', ncol(bar_plot_ls)*length(index_bar_plot))),
           topo=rep(rep(plot_unique_names, each=length(index_bar_plot)), 2),
           sigma=as.factor(rep(rep(exp(sig_sequence[index_bar_plot]), ncol(bar_plot_ls)), 2)),
           ren=c(c(bar_plot_ls), c(bar_plot_mle)),
           ren_sd=c(c(bar_plot_ls_sd, bar_plot_mle_sd)))

# reorder
factor_topo <- factor(dt_fasen$topo, level = c('RE-Europe 50', 'Polymer', 'Lattice', 'Complete'))
ggplot_infinite <- ggplot2::ggplot(data=dt_fasen, ggplot2::aes(x=sigma, y=ren, fill=type)) +
  ggplot2::ggtitle('Infinite jump activity') +
  ggplot2::facet_wrap(factor_topo, scales='free_x', ncol = 4) +
  ggplot2::geom_bar(stat="identity", position=ggplot2::position_dodge()) +
  ggplot2::scale_fill_manual(values=colors[c(3,1)]) + 
  ggplot2::geom_errorbar(ggplot2::aes(ymin=ren-1.96*ren_sd, ymax=ren+1.96*ren_sd), width=.5, size=1.2,
                position=ggplot2::position_dodge(.9)) +
  ggplot2::ylab(expression(paste('Relative  ', L^2, '-Error Norm'))) +
  ggplot2::xlab(expression(paste('Noise multiplicator ', sigma))) +
  ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank(), 
        strip.background = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_text(size = ggplot2::rel(multiplicator*1.3)),
        axis.title.x = ggplot2::element_text(size = ggplot2::rel(multiplicator*1.2)),
        strip.text.x = ggplot2::element_text(size = ggplot2::rel(multiplicator*.8), color = "black"),
        axis.text.x = ggplot2::element_text(size = ggplot2::rel(multiplicator*1.5), color = "black"),
        legend.title = ggplot2::element_text(size = ggplot2::rel(multiplicator*1.2), color = "black"),
        legend.text = ggplot2::element_text(size = ggplot2::rel(multiplicator*1.2), color = "black"),
        strip.text =  ggplot2::element_text(size = ggplot2::rel(multiplicator*1.0), color = "black"),
        axis.text.y = ggplot2::element_text(size = ggplot2::rel(multiplicator*1.5), color = "black"),
        plot.title = ggplot2::element_text(size = ggplot2::rel(multiplicator*1.5), color = "black", hjust = 0.5)) +
  ggplot2::labs(fill="Estimator")

gridExtra::grid.arrange(ggplot_finite, ggplot_infinite, ncol = 1)
dev.off()

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
std_dev_core_wind <- sqrt(diag(cov(core_wind)))

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
levy_increments_recovery <- LevyRecovery(
  fitted_adj = mle_theta_matrix,
  data = core_wind,
  times = recovery_times,
  look_ahead = 1)
ghyp_levy_recovery_fit <- FitLevyRecoveryDiffusion(levy_increments_recovery$increments)
finite_levy_recovery_fit <- FitBrownianMotionCompoundPoisson(
  data = levy_increments_recovery$increments,
  mesh_size = mesh_size,
  thresholds = rep(mesh_size^{0.4999}, n_nodes))
finite_levy_recovery_fit$poisson
finite_levy_recovery_fit$n_jumps
finite_levy_recovery_fit$jump_sigma
xlim_plus <-.2
alpha_plot <- .45
cex_lab <- 1.6
cex_axis <- 1.7
cex_main <- 1.7
clrs <- viridis::viridis(n_nodes, alpha = alpha_plot)

# setEPS()
# postscript("../data/pictures/exploring.eps")
png("../data/pictures/exploring.png", width = 1200, height = 450)
par(mfrow=c(1, 2), mar=c(4.5,4.5,2,2))
# Plot of densities
plot(
  density(levy_increments_recovery$increments[,1]),
  ylim=c(0, 35), xlim=c(-xlim_plus, xlim_plus),
  col=clrs[1], main=paste('Densities of noise marginals'),
  xlab='Increments', ylab='Density', lwd=2, bty = 'n',
  cex.main=cex_main, cex.axis=cex_axis, cex.lab=cex_lab
)
apply_idx <- 2
for(apply_idx in 2:n_nodes){
  levy_increments_recovery$increments[, apply_idx] %>% (function(x){density(x)}) %>% (function(x){lines(x, col=clrs[apply_idx], lwd=2)})
}

# Plot of QQ-Plots
alpha_plot <- .40
clrs <- viridis::viridis(n_nodes, alpha = alpha_plot)
quantile_values <- seq(from=0.005, to=.995, length.out = 100)
q_norm <- qnorm(quantile_values, mean = 0, sd=1)
q_norm <- matrix(rep(q_norm, n_nodes), ncol=n_nodes)
increments_quantiles <- apply(levy_increments_recovery$increments, 2, function(x){quantile(x/sd(x), quantile_values)})
matplot(x = q_norm, y = increments_quantiles,
        main=paste('Gaussian QQ plots on noise marginals'),
        xlab='Gaussian quantiles', ylab='Data quantiles', lwd=2, bty = 'n', pch=22, type='l',
        cex.main=cex_main, cex.axis=cex_axis, cex.lab=cex_lab, col = clrs)
lines(q_norm, q_norm, type='l', lwd=2, lty=1)
dev.off()


png("../data/pictures/finite.png", width = 1200, height = 450)
par(mfrow=c(1, 2), mar=c(4.5,4.5,2,0))
colnames(finite_levy_recovery_fit$jump_sigma) <- lapply(1:50, function(i){paste('N', i, sep='')})
rownames(finite_levy_recovery_fit$jump_sigma) <- lapply(1:50, function(i){paste('N', i, sep='')})

colnames(finite_levy_recovery_fit$sigma) <- lapply(1:50, function(i){paste('N', i, sep='')})
rownames(finite_levy_recovery_fit$sigma) <- lapply(1:50, function(i){paste('N', i, sep='')})

corrplot::corrplot(
  cov2cor(finite_levy_recovery_fit$jump_sigma),
  type="upper", col = viridis::viridis(100), # method="shade",
  sig.level = 0.01, insig = "blank", tl.cex=.75, tl.col='black',
  title=expression(paste('Weiner covariance matrix ', Sigma^J)), mar=c(0,0,2,0), cex.main=1.7)
# corrplot::corrplot(
#   cov2cor(finite_levy_recovery_fit$sigma),
#   type="upper", col = viridis::viridis(100), # method="shade",
#   sig.level = 0.01, insig = "blank", tl.cex=.75, tl.col='black',
#   title=expression(paste('Jump covariance matrix ', Sigma)), mar=c(0,0,2,0), cex.main=1.7)
# diag(finite_levy_recovery_fit$jump_sigma)/diag(finite_levy_recovery_fit$sigma)
# dev.off()

N_plot <- 5000
finite_samples <- do.call(cbind, lapply(seq_len(length(finite_levy_recovery_fit$n_jumps)),
       function(j){
         set.seed(43)
          finite_samples <- BrownianMotionCompoundPoisson(
            n = N_plot,
            n_jumps = finite_levy_recovery_fit$n_jumps[[j]],
            sigma = finite_levy_recovery_fit$sigma,
            jump_sigma = finite_levy_recovery_fit$jump_sigma,
            delta_time = mesh_size)[,j]
          return(finite_samples)
    }
))

finite_samples <- BrownianMotionCompoundPoisson(
  n = N_plot,
  n_jumps = round(median(unlist(finite_levy_recovery_fit$n_jumps))),
  sigma = finite_levy_recovery_fit$sigma,
  jump_sigma = finite_levy_recovery_fit$jump_sigma,
  delta_time = mesh_size)

finite_quantiles <- apply(finite_samples, 2, function(x){quantile(x/sd(x), quantile_values)})
matplot(x =finite_quantiles , y = increments_quantiles,
        main=paste('Gaussian QQ plots on noise marginals'),
        xlab='Generated quantiles', ylab='Data quantiles', lwd=2, bty = 'n', pch=22, type='l',
        cex.main=cex_main, cex.axis=cex_axis, cex.lab=cex_lab, col = clrs)
lines((-10):10, (-10):10)
dev.off()


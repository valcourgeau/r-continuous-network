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


set.seed(42)
n_paths <- 10
N <- 2000
levy_increment_sims <- list()
for(i in 1:n_paths){
  levy_increment_sims[[i]]<- matrix(ghyp::rghyp(n = N, object = ghyp_levy_recovery_fit$FULL), nrow=N)
}

DO_PARALLEL <- TRUE

# choosing network type
network_types <- list(PolymerNetwork, LatticeNetwork, FullyConnectedNetwork, adj_grid)
network_types_name <- c('polymer_mles', 'lattice_mles', 'fc_mles', 're_europe_mles')

index_network <- 1
beta_study <- list()

theta_1 <- mle_theta_vector[1]
theta_2 <- mle_theta_vector[2]

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
  
  if(DO_PARALLEL){
    num_cores <- detectCores()-1
    cl <- makeCluster(num_cores)
    clusterExport(cl, varlist=c("ConstructPath", "network_topo", "mesh_size", "first_point", "levy_increment_sims",
                                "network_topo_raw", "beta", "mesh_size", "n_nodes", "GrouMLE",
                                "NodeMLELong", "CoreNodeMLE", "RowNormalised"))
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
  beta_seq <- seq(0.001, to = 0.499, length.out = 10)
  n_row_generated <- nrow(generated_paths[[1]])
  
  if(DO_PARALLEL){
    clusterExport(cl, varlist=c("n_row_generated", "beta_seq"))
    generated_fit <- parLapply(
      cl,
      generated_paths,
      function(x){
        lapply(
          beta_seq,
          function(beta_value){
            GrouMLE(times = seq(0, length.out = n_row_generated, by = mesh_size),
                    adj = as.matrix(network_topo),
                    data = x,
                    thresholds = rep(mesh_size^beta_value, n_nodes), # mesh_size^beta
                    mode = 'network',
                    output = 'vector')
          }
        )
      }
    )
    stopCluster(cl)
  }else{
    generated_fit <- lapply(
      X = generated_paths,
      FUN =
        function(x){
          lapply(
            beta_seq,
            function(beta_value){
              GrouMLE(times = seq(0, length.out = n_row_generated, by = mesh_size),
                      adj = as.matrix(network_topo),
                      data = x,
                      thresholds = rep(mesh_size^beta_value, n_nodes), # mesh_size^beta
                      mode = 'network',
                      output = 'vector')
            }
          )
        }
    )
  }
  print('fit generated')
  beta_study[[index_network]] <- generated_fit
  index_network <- index_network + 1
}

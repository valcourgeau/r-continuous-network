# loading the packages
source("~/GitHub/r-continuous-network/R/package_to_load.R")

# We do not use the load data but solar since load is too correlated and not very random
AS_SPARSE <- FALSE

# Functions and procedures to clean the data
data_path <- "~/GitHub/r-continuous-network/data/re-europe/"
n_df_load <- 1000
n_nodes <- 10
df_load <- read.csv(paste(data_path, "Nodal_TS/wind_signal_COSMO.csv", sep=""), nrows = n_df_load+10)[,2:(n_nodes+1)]
df_load <- df_load[-c(1:10),]

clean_wind_data <- CleanData(df_load)
core_wind <- clean_wind_data$remainders

# Network topology
load_topo <- read.csv(file=paste(data_path, "Static_data/network_edges.csv", sep=""))
load_nodes <- read.csv(file=paste(data_path, "Static_data/network_nodes.csv", sep=""))
load_topo$fromName <- load_nodes$name[load_topo$fromNode]
load_topo$toName <- load_nodes$name[load_topo$toNode]

forgotten_nodes <- (1:length(load_nodes$ID)) [which(! 1:length(load_nodes$ID) %in% c(load_topo$fromNode, load_topo$toNode))]
edges_list <- data.frame("from"=load_topo$fromNode, "to"=load_topo$toNode)
graph_grid <- graph.data.frame(
  edges_list
)
adj_grid <- igraph::as_adjacency_matrix(graph = graph_grid, sparse = AS_SPARSE)
for (i in 1:ncol(adj_grid)){ # add self-loop
  adj_grid[i,i] <- 1
}
adj_grid

# EXPLORATION PLOT
par(mfrow=c(1,1))
clrs <- list(color = colorRampPalette(brewer.pal(11,"Spectral"))(n_nodes))
plot(core_wind[,1], type="l", ylim=c(-0.4,0.53), ylab="Hourly load in MWh",
     main="Load across nodes", col = clrs$color[1])
for(i in 2:n_nodes){
  lines(core_wind[,i], type="l", col=clrs$color[i])
}

FasenRegression(core_wind) * adj_grid[1:n_nodes, 1:n_nodes]
norm(FasenRegression(core_wind) * adj_grid[1:n_nodes, 1:n_nodes], 'F')
vapply(seq(0.001,10,length.out=10),
       function(x){norm(
          FasenRegression(
            core_wind * matrix(rnorm(n_df_load*n_nodes, mean = 1.0, sd = x), nrow = n_df_load, ncol = n_nodes)
          ) * adj_grid[1:n_nodes, 1:n_nodes], 'F')},
       0.1)

# TODO
warning('get NOU + compare with Fasen')



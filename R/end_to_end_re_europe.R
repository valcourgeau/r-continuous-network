# loading the packages
source("~/GitHub/r-continuous-network/R/package_to_load.R")

# We do not use the load data but wind since load is too correlated and not very random
AS_SPARSE <- FALSE

#######################################################################
##################### DATA PREPARATION ################################
#######################################################################

# Functions and procedures to clean the data
data_path <- "~/GitHub/r-continuous-network/data/re-europe/"
n_df_load <- 10000
n_nodes <- 10
df_load <- read.csv(paste(data_path, "Nodal_TS/wind_signal_COSMO.csv", sep=""), nrows = n_df_load+10)[,2:(n_nodes+1)]
df_load <- df_load[-c(1:10),]

clean_wind_data <- CleanData(df_load)
core_wind <- clean_wind_data$remainders

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
topo_graph <- igraph::graph.data.frame(d = topo_edges, directed = TRUE, vertices = topo_nodes)

adj_grid <- igraph::as_adjacency_matrix(topo_graph, sparse = AS_SPARSE)

# EXPLORATION PLOT
par(mfrow=c(1,1))
clrs <- list(color = colorRampPalette(brewer.pal(11,"Spectral"))(n_nodes))
plot(core_wind[,1], type="l", ylim=c(-0.4,0.53), ylab="Hourly load in MWh",
     main="Load across nodes", col = clrs$color[1])
for(i in 2:n_nodes){
  lines(core_wind[,i], type="l", col=clrs$color[i])
}

#######################################################################
###################### FASEN COMPARISON ###############################
#######################################################################

FasenRegression(core_wind) * (adj_grid + diag(n_nodes))
GrouMLE(times=seq(0, by=1, length.out = n_sample), data=core_wind, adj = adj_grid, div = 1e3, mode="node", output = "matrix")

norm(FasenRegression(core_wind) * adj_grid[1:n_nodes, 1:n_nodes], 'F')
vapply(seq(0.001,10,length.out=10),
       function(x){norm(
          FasenRegression(
            core_wind * matrix(rnorm(n_df_load*n_nodes, mean = 1.0, sd = x), nrow = n_df_load, ncol = n_nodes)
          ) * adj_grid[1:n_nodes, 1:n_nodes], 'F')},
       0.1)




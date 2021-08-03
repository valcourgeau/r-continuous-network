# loading the packages
# source("~/GitHub/r-continuous-network/R/package_to_load.R")
# source("~/GitHub/r-continuous-network/R/utils.R")
# source("~/go/src/r-continuous-network/R/adjacency_generation.R")
# 
# source("~/GitHub/r-continuous-network/R/levy_recovery_v2.R")
# source("~/GitHub/r-continuous-network/R/utils.R")
# source("~/GitHub/r-continuous-network/R/path_generation.R")
# source("~/GitHub/r-continuous-network/R/grou_mle.R")

library(ntwk)
# We do not use the load data but wind since load is too correlated and not very random
AS_SPARSE <- F

#######################################################################
##################### DATA PREPARATION ################################
#######################################################################

# Functions and procedures to clean the data
data_path <- "~/Downloads/RE-Europe_dataset_package/"
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
load_nodes <- read.csv(file=paste(data_path, "Metadata/network_nodes.csv", sep=""))
load_nodes <- load_nodes[1:n_nodes,]
topo_nodes <- data.frame("name"=load_nodes$ID, 
                         "lon"= load_nodes$longitude,
                         "lat"=load_nodes$latitude)
load_edges <- read.csv(file=paste(data_path, "Metadata/network_edges.csv", sep=""))
load_edges <- load_edges[which(load_edges$fromNode %in% 1:n_nodes & 
                                 load_edges$toNode %in% 1:n_nodes),]
topo_edges <- data.frame("from" = load_edges$fromNode, 
                         "to" = load_edges$toNode)
topo_graph <- igraph::graph.data.frame(d = topo_edges, directed = FALSE, vertices = topo_nodes)

adj_grid <- igraph::as_adjacency_matrix(topo_graph, sparse = AS_SPARSE)
adj_grid <- as.matrix(adj_grid)
mesh_size <- 2/24

europeanUnion <- c("Spain", "Portugal")
maps::map(region=europeanUnion, col="grey80", fill=TRUE, bg="white", lwd=0.1)
igraph::plot.igraph(x = topo_graph, add=T, rescale=F,
            layout=topo_nodes[,2:3], vertex.label=NA, arrow.size=0.4, edge.color='black')

visNetwork::visNetwork(
  nodes=topo_nodes, edges=topo_edges
)

plot(igraph::graph_from_adjacency_matrix(adjmatrix = adj_grid))
igraph::layout_on_grid(igraph::graph_from_adjacency_matrix(PolymerNetwork(10, 2, 4)))
plot(igraph::layout_on_grid(igraph::graph_from_adjacency_matrix(PolymerNetwork(50, 2, 4)), width = 7, height = 8))
plot(igraph::graph_from_adjacency_matrix(PolymerNetwork(50, 2, 4)))

pdf("data/pictures/pdf_network_configurations_v2.pdf", width = 7, height = 2)
par(mfrow=c(1,4))
node.w <- 1
node.h <- 1

tmp <- PolymerNetwork(50, 2, 4)
tmp[lower.tri(tmp, diag = F)] <- 0
qgraph::qgraph(tmp, directed=F, parallelEdge=T, weighted=F, vTrans=255, 
               labels=F, node.width=node.w, node.height=node.h, layout='circular', edge.color='black',
               edge.width=1.5, esize=2, title='Polymer', title.cex=1.5)

tmp <- LatticeNetwork(50, 2, 4)
tmp[lower.tri(tmp, diag = F)] <- 0
qgraph::qgraph(tmp, directed=F, parallelEdge=T, weighted=F, vTrans=255, 
               labels=F, node.width=node.w, node.height=node.h, layout='spring', edge.color='black',
               edge.width=1.5, esize=2, title='Lattice', title.cex=1.5)

tmp <- FullyConnectedNetwork(50, 2, 4)
tmp[lower.tri(tmp, diag = F)] <- 0
qgraph::qgraph(tmp, directed=F, parallelEdge=T, weighted=F, vTrans=255, 
               labels=F, node.width=node.w, node.height=node.h, layout='spring', trans=0.2, 
               edge.width=1.5, esize=2, title='Fully-Connected', title.cex=1.5)

qgraph::qgraph(adj_grid, 
               directed=F, parallelEdge=T, weighted=F, vTrans=255, 
               labels=F, node.width=node.w, node.height=node.h, title.cex=1.5, edge.color='black',
               edge.width=1.5, esize=2, title='RE-Europe 50')
dev.off()


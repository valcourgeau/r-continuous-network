# loading the packages
setwd("~/GitHub/r-continuous-network/R/")
source("package_to_load.R")

# Functions and procedures to clean the data
setwd("~/GitHub/r-continuous-network/data/RE-Europe_dataset_package/")

#### Network topography
n_edges <- 1200
set.seed(42)

setwd("Metadata/")
load_topo <- read.csv(file="network_edges.csv")
load_nodes <- read.csv(file="network_nodes.csv")
load_topo$fromName <- load_nodes$name[load_topo$fromNode]
load_topo$toName <- load_nodes$name[load_topo$toNode]

forgotten_nodes <- (1:length(load_nodes$ID)) [which(! 1:length(load_nodes$ID) %in% c(load_topo$fromNode, load_topo$toNode))]
edges_list <- data.frame("from"=load_topo$fromNode, "to"=load_topo$toNode)
graph_grid <- graph.data.frame(
  edges_list
  )

# TODO Change labels

# clean_edges <- rep("o", length(V(graph_grid))) %>% as.character
# clean_edges[as.integer(load_nodes$ID)] <- as.character(load_nodes$name)
# missing_labels <- (1:length(V(graph_grid)))[-load_nodes$ID]
# clean_edges[-load_nodes$ID] <- missing_labels
# clean_edges <- factor(clean_edges, levels = clean_edges, labels=clean_edges)
# clean_edges %>% length
# clean_edges
# ggs <- set_vertex_attr(graph_grid, "label",  value=clean_edges)
# V(ggs)$label <- clean_edges
# plot(ggs)

adj_grid <- as_adjacency_matrix(graph = graph_grid)

setwd("../Nodal_TS/")
#n_cols <- 1200
load_data <- read.csv(file = "load_signal.csv", nrows = 10000)
colnames(load_nodes)[1] <- "ds"
head(load_nodes)

load_network <- network(as_adj(simplify(graph_grid)))

ideg <- degree(load_network, cmode="indegree")
odeg <- degree(load_network, cmode="outdegree")

load_network %>% plot
set.vertex.attribute(x = load_network, 
                     attrname = "load",
                     value = load_data[1,-1])
load_network %>% plot
get.vertex.attribute(load_network, attrname = "load")
plot(load_network, vertex.col=rgb(odeg/max(odeg),0,ideg/max(ideg)), 
     displayisolates = T, vertex.cex=3*pmax(odeg/max(odeg), ideg/max(ideg)),
     main="Red for out-degree, blue for in-degree")

mean_load <- colMeans(load_data[,-1])
mean_load_log <- mean_load %>% log
mean_load_log[which(mean_load_log == -Inf)] <- -1
mean_load_log <- mean_load_log - min(mean_load_log)
mean_load_log <- mean_load_log / sd(mean_load_log)
plot(mean_load_log)

par(mfrow=c(1,3))
layout(matrix(c(1,1,2,2,3,3, 1,1,2,2,4,4, 1,1,2,2,5,5), 3, 6, byrow = TRUE))
plot(load_network, vertex.col=rgb(0,mean_load/max(mean_load),0), 
     displayisolates = T, vertex.cex=3*pmax(odeg/max(odeg), ideg/max(ideg)),
     main="Greener = higher mean load")
plot(load_network, vertex.col=rgb(0,mean_load_log/max(mean_load_log),0), 
     displayisolates = T, vertex.cex=3*pmax(odeg/max(odeg), ideg/max(ideg)),
     main="Greener = higher mean load [LOG]")
plot(mean_load_log, ylab="rescaled log of mean load", xlab="Node",
     pch = 20, main="Rescaled log of mean load")
plot(sort(mean_load_log, decreasing = F), 
     ylab="Log of mean load",
     main="Sorted log of mean load",
     xlab="Sorted nodes by log of mean load")
plot(density(mean_load_log), 
     ylab="",
     main="PDF of log of mean-load",
     xlab="Log of mean-load")




# Prophet part
df <- load_nodes[,1:2]
colnames(df) <- c("ds", "y")
df_prophet <- prophet(df) # PROPHET

par(xpd=FALSE)
par(mfrow=c(1,1))

plot(df$y, type="l")
abline(v=which(df$ds %in% as.character(df_prophet$changepoints)))
diff(which(df$ds %in% as.character(df_prophet$changepoints)))/24
# changes every months

prophet:::plot_weekly(df_prophet)
prophet:::plot(df_prophet)


# STL
library("stats")
ts.x1 <- ts(frequency = 180*24, load_data$X1)
ts.x1.stl <- stl(ts.x1, s.window = "per")
plot(ts.x1.stl$time.series[1:1000,3], type="l")
abline(v=(1:300)*7*24)

NOUfit1D(times = 1:length(load_data$X1), 
         data = ts.x1.stl$time.series[,3], threshold = 0.6)

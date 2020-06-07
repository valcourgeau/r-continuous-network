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

mle_theta <- GrouMLE(times=seq(0, by=1, length.out = n_sample), data=core_wind, adj = adj_grid, div = 1e3, mode="network", output = "matrix")
mle_theta



par(mfrow=c(2,5), mar=c(4,4.5,2,4))
cex_value <- 1.8
for(i in c(10,40)){
  i <- round(i)
  print(i)
  
  quantil_vals <- 1:4999/5000
  
  val_x_qq <- quantile(x = lr$increments[,i], quantil_vals)
  set.seed(42)
  val_y_qq <- ghyp::rghyp(object = ghyp_fit_re_europe_50$NIG, n = length(lr$increments[,i]))[, i]
  val_y_qq <- quantile(x = val_y_qq, quantil_vals)
  qqplot(x = val_x_qq, y = val_y_qq, 
         ylab = 'NIG', xlab= paste('Levy increments', i), 
         cex.lab = cex_value, cex.axis=cex_value)
  
  if(i == 10){title('NIG', cex.main=cex_value, line = 0.8)}
  
  set.seed(42)
  val_y_qq <- ghyp::rghyp(object = ghyp_fit_re_europe_50$GAUSS, n = length(lr$increments[,i]))[, i]
  val_y_qq <- quantile(x = val_y_qq, quantil_vals)
  qqplot(x = val_x_qq, y = ghyp::rghyp(object = ghyp_fit_re_europe_50$GAUSS, n = length(lr$increments[,i]))[, i], 
         ylab = 'Gaussian', xlab= paste('Levy increments', i), 
         cex.lab = cex_value, cex.axis=cex_value)
  if(i == 10){title('Gaussian', cex.main=cex_value, line = 0.8)}
  
  set.seed(42)
  val_y_qq <- ghyp::rghyp(object = ghyp_fit_re_europe_50$VG, n = length(lr$increments[,i]))[, i]
  val_y_qq <- quantile(x = val_y_qq, quantil_vals)
  qqplot(x = val_x_qq, y = ghyp::rghyp(object = ghyp_fit_re_europe_50$VG, n = length(lr$increments[,i]))[, i], 
         ylab = 'VG', xlab= paste('Levy increments', i), 
         cex.lab = cex_value, cex.axis=cex_value)
  if(i == 10){title('VG', cex.main=cex_value, line = 0.8)}
  
  set.seed(42)
  val_y_qq <- ghyp::rghyp(object = ghyp_fit_re_europe_50$T, n = length(lr$increments[,i]))[, i]
  val_y_qq <- quantile(x = val_y_qq, quantil_vals)
  qqplot(x = val_x_qq, y = ghyp::rghyp(object = ghyp_fit_re_europe_50$T, n = length(lr$increments[,i]))[, i], 
         ylab = 'Student\'s t', xlab= paste('Levy increments', i),
         cex.lab = cex_value, cex.axis=cex_value)
  if(i == 10){title('Student\'s t', cex.main=cex_value, line = 0.8)}
  
  set.seed(42)
  val_y_qq <- ghyp::rghyp(object = ghyp_full, n = length(lr$increments[,i]))[, i]
  val_y_qq <- quantile(x = val_y_qq, quantil_vals)
  qqplot(x = val_x_qq, y = ghyp::rghyp(object = ghyp_full, n = length(lr$increments[,i]))[, i], 
         ylab = 'GHYP', xlab= paste('Levy increments', i), 
         cex.lab = cex_value, cex.axis=cex_value)
  if(i == 10){title('GHYP', cex.main=cex_value, line = 0.8)}
}

library("igraph")
library("raster")
setwd("~/GitHub/r-continuous-network/data/RE-Europe_dataset_package/")
setwd("Nodal_TS/")

#### Network topography
n_edges <- 1200
set.seed(42)

setwd("../Metadata/")
load_nodes <- read.csv(file = "network_nodes.csv")
#load_nodes <- load_nodes[-which(topo_nodes %in% c("Varna", "FromUkrainaWest")),]
colnames(load_nodes)
load_nodes <- load_nodes[1:n_edges,]
topo_nodes <- data.frame("name"=load_nodes$ID, 
                          "lon"= load_nodes$longitude,
                         "lat"=load_nodes$latitude
                         )
load_nodes

load_edges <- read.csv(file = "network_edges.csv")
colnames(load_edges)
load_edges <- load_edges[which(load_edges$fromNode %in% 1:n_edges & 
                                 load_edges$toNode %in% 1:n_edges),]
topo_edges <- data.frame("from" = load_edges$fromNode, 
                          "to" = load_edges$toNode)
topo_graph <- graph.data.frame(d = topo_edges, directed = TRUE, vertices = topo_nodes)

# European Union Map
# Member States of the European Union
europeanUnion <- c("Austria","Belgium","Bulgaria","Croatia","Cyprus",
                   "Czech Rep.","Denmark","Estonia","Finland","France",
                   "Germany","Greece","Hungary","Ireland","Italy","Latvia",
                   "Lithuania","Luxembourg","Malta","Netherlands","Poland",
                   "Portugal","Romania","Slovakia","Slovenia","Spain",
                   "Sweden","United Kingdom")

map(region=europeanUnion, col="grey80", fill=TRUE, bg="white", lwd=0.1)
plot.igraph(x = topo_graph, add=T, rescale=F,
            layout=topo_nodes[,2:3], vertex.label=NA, arrow.size=0.1)


### Numerics
setwd("../Nodal_TS/")
n_df_load <- 500
n_nodes <- 70
df_load <- read.csv("load_signal.csv", nrows = n_df_load)[,1:(n_nodes+1)]
#par(mfrow=c(round(sqrt(n_nodes)) ,round(sqrt(n_nodes+1))))
par(mfrow=c(1,1))
plot(df_load[,2], type="l", ylim=c(0,1800), ylab="Hourly load in MWh",
     main="Load across nodes")
for(i in 2:n_nodes){
  lines(df_load[,1+i], type="l")
}

n_df_load <- 72
n_nodes <- 70
df_load <- read.csv("load_signal.csv", nrows = n_df_load)[,1:(n_nodes+1)]
#par(mfrow=c(round(sqrt(n_nodes)) ,round(sqrt(n_nodes+1))))
par(mfrow=c(1,1))
library(colorspace)
library(RColorBrewer)
clrs <- list(color = colorRampPalette(brewer.pal(11,"Spectral"))(n_nodes))
plot(df_load[,2], type="l", ylim=c(0,1800), ylab="Hourly load in MWh",
     main="Load across nodes", col = clrs$color[1])
for(i in 2:n_nodes){
  lines(df_load[,1+i], type="l", col=clrs$color[i])
}

# study of one time series
load_nodes[load_nodes$ID == 50,]
library(colorspace)
par(mfrow=c(1,1))
n_df_load <- 1500
df_load <- read.csv("load_signal.csv", nrows = n_df_load)[,1:(n_nodes+1)]
plot(df_load[,50], type="l")
plot(df_load[1:(24*10),50], type="l")
abline(v = c(24*0:10+4))

library(zoo)
par(mfrow=c(2,5))
d <- sort(diff(df_load[10:1500,50]-rollapply(df_load[,50], 10, mean)))
hist(d,
     breaks=35, col=rgb(0.8*10/10,0.5,0.8*10/10,0.8), probability = T,
     main=paste("Incr of data-mean", 10))
lines(d, dnorm(x = d, mean = mean(d), sd = sd(d)), lwd=2)
for(i in 9:2){
  d <- sort(diff(df_load[i:1500,50]-rollapply(df_load[,50], i, mean)))
  hist(d, breaks=35, col=rgb(0.8*i/10,0.5,0.8*i/10,0.8), probability = T,
       main=paste("Incr of data-mean", i))
  lines(d, dnorm(x = d, mean = mean(d), sd = sd(d)), lwd=2)
}

  d <- sort(diff(df_load[i:1500,50]))
  hist(d, breaks=35, col=rgb(0.8*i/10,0.5,0.8*i/10,0.8), probability = T,
       main=paste("Incr of data-mean", 1))
  lines(d, dnorm(x = d, mean = mean(d), sd = sd(d)), lwd=2)
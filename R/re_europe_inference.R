library("igraph")
library("raster")
library(magrittr)
library("data.table")


setwd("~/GitHub/r-continuous-network/data/RE-Europe_dataset_package/")
setwd("Nodal_TS/")

#### Network topography
n_edges <- 50 # this can be easily increased to 1000s
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
df_load <- read.csv("solar_signal_COSMO.csv", nrows = n_df_load)[,1:(n_nodes+1)]
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
n_df_load <- 50000
n_nodes <- 50
# df_load <- read.csv("~/GitHub/r-continuous-network/data/RE-Europe_dataset_package/Nodal_TS/wind_signal_COSMO.csv", 
#                     nrows = n_df_load)[,1:(n_nodes+1)]
df_load <- data.table::fread("~/GitHub/r-continuous-network/data/RE-Europe_dataset_package/Nodal_TS/wind_signal_COSMO.csv", 
                             stringsAsFactors = F)
df_load <- df_load[,1:(n_nodes+1)]
df_load <- df_load[-1,]
dim(df_load)
par(mfrow=c(1,1))

plot(as.data.frame(df_load)$V50)

library(zoo)
par(mfrow=c(2,5))
d <- sort(diff(df_load[10:1500,50]-rollapply(df_load[,50], 10, mean)))
hist(d,
     breaks=100, col=rgb(0.8*10/10,0.5,0.8*10/10,0.8), probability = T,
     main=paste("Incr of data-mean", 10))
lines(d, dnorm(x = d, mean = mean(d), sd = sd(d)), lwd=2)
for(i in 9:2){
  d <- sort(diff(df_load[i:1500,50]-rollapply(df_load[,50], i, mean)))
  hist(d, breaks=35, col=rgb(0.8*i/10,0.5,0.8*i/10,0.8), probability = T,
       main=paste("Incr of data-mean", i))
  lines(d, dnorm(x = d, mean = mean(d), sd = sd(d)), lwd=2)
}


#JumpTest1D(data=df_load[1:5000,50], )

d_n <- 1/(24*3600)
plot(seq(0,1/2,length.out = 1000), vapply(seq(0,1/2,length.out = 1000), 
       function(x){length(which(absolute_increments < d_n^(x)))/length(absolute_increments)}, 1))
v_n <- d_n^0.3

whole_stream <- df_load[, 50]

absolute_increments <- diff(whole_stream) %>% abs
signed_increments <- diff(whole_stream)
filter_indices <- which(absolute_increments >= v_n) # finds 'jumps' with filter
cts_stream <- rep(0, length(whole_stream))
jump_stream <- rep(0, length(whole_stream))

cts_stream <- whole_stream
cts_stream[filter_indices+1] <- cts_stream[filter_indices+1] - signed_increments[filter_indices]
jump_stream[filter_indices+1] <- signed_increments[filter_indices]

plot(whole_stream[1:200], type='b')
lines(cts_stream[1:200], col='red')
lines(cts_stream[1:200]+jump_stream[1:200], col='green')
for(k in filter_indices){
  lines(c(k,k), c(0,1), col = 'blue')#rgb(1,1,0.5,alpha=0.1) )# rgb(0,0,0.5,alpha=0.9) )
}

rollapply(signed_increments[(absolute_increments >= v_n)] %>% as.numeric() %>% diff, FUN=sd, width=1000) %>% (function(x){plot(x,ylim=c(0,.5), main='Empirical SD: test for time-homogeneous Compound Poisson')})

nig_fit <- ghyp::fit.NIGmv(data = df_load[1:10000,2:10])
plot(density(ghyp::rghyp(n = 5000, nig_fit)[,1]))
hist(df_load[1:20000,2], breaks=100, add=T, probability=T)

ghyp::lik.ratio.test(df_load[1:1500,2:10], x.subclass = nig_fit)

library(NHPoisson)
library(cplm)

plot(signed_increments[(absolute_increments >= v_n)] %>% as.numeric() %>% diff %>% 
  (function(x){x^2}) %>% cumsum()/(2:length(which((absolute_increments >= v_n))) ))

(absolute_increments >= v_n) %>% as.numeric() %>% cumsum %>%  diff
signed_increments[(absolute_increments >= v_n)] %>% as.numeric() %>% diff %>% (function(x){hist(x[1000:2000], breaks=150)})
signed_increments[(absolute_increments >= v_n)] %>% as.numeric() %>% diff %>% mean
signed_increments[(absolute_increments >= v_n)] %>% as.numeric() %>% diff %>% sd
MASS::fitdistr(signed_increments[(absolute_increments >= v_n)] %>% as.numeric() %>% diff, densfun = 'normal')
lambda <- var(jump_stream) / mean(signed_increments[(absolute_increments >= v_n)] %>% as.numeric() %>% diff %>% (function(x){x^2})) / (length(jump_stream)*d_n)



d <- sort(diff(df_load[i:1500,50]))
hist(d, breaks=100, col=rgb(0.8*i/10,0.5,0.8*i/10,0.8), probability = T,
     main=paste("Incr of data-mean", 1))
lines(d, dnorm(x = d, mean = mean(d), sd = sd(d)), lwd=2)

library(ghyp)
ghyp::fit.NIGmv(data = df_load[1:1500,2:10])


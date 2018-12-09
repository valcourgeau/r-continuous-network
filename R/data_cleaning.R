# Functions and procedures to clean the data

library("prophet")
library("igraph")
library("raster")
setwd("~/GitHub/r-continuous-network/data/RE-Europe_dataset_package/")

#### Network topography
n_edges <- 1200
set.seed(42)

setwd("Nodal_TS/")
n_cols <- 20
load_nodes <- read.csv(file = "load_signal.csv", nrows = 10000)[,1:(n_cols+1)]
colnames(load_nodes)[1] <- "ds"
df <- load_nodes[,1:2]
colnames(df) <- c("ds", "y")
df_prophet <- prophet(df) # PROPHET

par(xpd=FALSE)
plot(df$y, type="l")
abline(v=which(df$ds %in% as.character(df_prophet$changepoints)))
diff(which(df$ds %in% as.character(df_prophet$changepoints)))/24
# changes every months

prophet:::plot_weekly(df_prophet)
prophet:::plot(df_prophet)


# STL
library("stats")
ts.x1 <- ts(frequency = 365*24, load_nodes$X1)
ts.x1.stl <- stl(ts.x1, s.window = "per")
plot(ts.x1.stl)
plot(ts.x1.stl$time.series[1:1000,3], type="l")
abline(v=(1:300)*7*24)

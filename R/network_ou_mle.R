### Functions to compute cts-time NAR-equivalent MLE and drift MLE

EnsureTopo <- function(nw_topo){
  # creates matrix with entries either 0 or 1.
  nw_topo[nw_topo != 0] <- 1
  return(nw_topo)
}

# Example
nw_topo <- c(2,1,1,
             0,2,3,
             0,0,-1)
nw_topo <- matrix(nw_topo, ncol = 3, nrow = 3)
nw_topo
EnsureTopo(nw_topo = nw_topo)

test_topo <- matrix(1, ncol=3, nrow=3)
EnsureTopo(nw_topo = test_topo)

StdTopo <- function(nw_topo){
  nw_topo <- EnsureTopo(nw_topo = nw_topo)
  diag(nw_topo) <- 0 
  std_coeff <- rowSums(nw_topo)
  std_coeff[std_coeff == 0] <- 1

    return(nw_topo/std_coeff)
}

# Example
nw_topo <- c(2,1,1,
             0,2,3,
             0,0,1)
nw_topo <- matrix(nw_topo, ncol = 3, nrow = 3)
StdTopo(nw_topo = nw_topo)

test_topo <- matrix(1, ncol=3, nrow=3)
StdTopo(nw_topo = test_topo)

DataFiltering <- function(data, thresholds, diff=F){
  # Shall we have the first or second point when diff=F?
  filtered_data <- diff(data)
  filtered_data <= thresholds
  
  under_thres <- t(apply(filtered_data, MARGIN = 1, '<=', thresholds))
  
  if(diff){
    return(filtered_data*under_thres)
  }else{
    #under_thres <- rbind(under_thres, rep(F, length(data[1,])))
    return(matrix(data[which(under_thres)], 
                  nrow=length(data[,1])-1))
  }
}

# Example
data <- c(1,2,5,
          1,3,6,
          1,4,7)
data <- t(matrix(data, ncol=3))
data
DataFiltering(data, c(2,3,6))
DataFiltering(data, c(2,3,6), diff = T)

TimeMatrix <- function(times, ncol, diff=F){
  if(diff){
    times <- diff(times)
  }
  
  return(t(matrix(rep(times, each=ncol), nrow=ncol)))
}

# Example
times <- c(1,2,4,7,11)
TimeMatrix(times, ncol = 2, diff=F)
TimeMatrix(times, ncol = 2, diff=T)

NOUfit <- function(nw_topo, times, data, thresholds){
  # fitting MLE to the network
  # TODO add support for jumps by using s- instead of s
  
  # Check 
  nw_topo <- StdTopo(nw_topo = nw_topo)
  nw_test <- nw_topo
  diag(nw_test) <- 1 
  if(max(Re(eigen(nw_test)$values)) <= 0){
    stop("Q matrix should be positive definite.")
  }
  
  d <- length(data[1,])
  N <- length(data[,1])

  if(length(thresholds) != d){
    stop('Wrong dimension between data and thresholds.')
  }
  if(length(times) != N){
    stop('Wrong dimensions between data and times.')
  }
  
  diff_filtered <- DataFiltering(data, thresholds, diff=T)
  diff_times <- TimeMatrix(times, ncol = d, diff=T)
  
  a_bar_y_t <- t(nw_topo %*% t(data))
  a_bar_y_t <- a_bar_y_t[-N,]
  y_t <- data[-N,]

  wide_mle_xi <- - matrix(c(
      sum(a_bar_y_t^2*diff_times), sum(a_bar_y_t*y_t*diff_times),
      sum(a_bar_y_t*y_t*diff_times), sum(y_t^2*diff_times)
    ), ncol=2)
  wide_mle_c_filtered <- c(
      sum(a_bar_y_t * diff_filtered),
      sum(y_t * diff_filtered)
    ) 
  
  #return(solve(a = wide_mle_xi, b = wide_mle_c_filtered))
  results <- list("MLE_wide"=
                    as.vector(solve(a = wide_mle_xi, b = wide_mle_c_filtered)))
  
  nw_topo <- results[["MLE_wide"]][1] * nw_topo
  diag(nw_topo) <- results[["MLE_wide"]][2]
  results[["Q"]] <- nw_topo
  
  return(results)
} 

# Example
#source("R/network_generation.R")

{
  d <- 5 #dims
  N <- 24*15 #numbers of points
  Y0 <- 1 #start point
  
  delta_t <- 1e-2
  
  sigma <- 1
  set.seed(42)
  nw_topo <- genRdmAssymetricGraphs(d = d, p.link = 0.25,
                         theta_1 = 1, theta_2 = 0)
  nw_topo <- StdTopo(nw_topo)
  
  times <- seq(from = 0, by = delta_t, length.out = N)
  nw_data_bm <- matrix(rnorm(n = d*N, mean = 0, sd = sigma*sqrt(delta_t)), ncol = d)
  nw_data <- matrix(0, ncol=d, nrow=N)
  nw_data[1,] <- 1:d/d
  
  nw_q <- 0.1*nw_topo
  diag(nw_q) <- 1
  
  for(index in 2:N){
    nw_data[index,] <-  -(nw_q %*% nw_data[index-1,]) * (times[index]-times[index-1]) + nw_data_bm[index,]
  }
  
  plot(nw_data[1:300,1], type="l")
  
  #nw_data <- Y0 + apply(X = nw_data, MARGIN = 2, FUN = function(x){cumsum(x)})
  thresholds <- rep(10, d)
  
  # Removing n.horizon points
  n.horizon <- 15
  nw_data_all <- nw_data
  nw_data <- nw_data[-((N-n.horizon+1):N),]
  times <- times[-((N-n.horizon+1):N)]
  N_all <- N
  N <- N-n.horizon
  
  res <- NOUfit(nw_topo = nw_topo, times = times,
                data = nw_data, thresholds = thresholds)
  res$MLE_wide
}

library("Matrix")
#set.seed(42)
n.subgrid <- 1000

levy_noise <- rnorm(n = n.horizon*n.subgrid*d, 
                    mean = 0, sd = sigma*sqrt(1/(24*n.subgrid)))
levy_noise <- matrix(levy_noise, ncol = d)
temp <- nw_data[N,]
pred <- matrix(0, ncol=d, nrow=n.horizon)
for(i in 1:n.horizon){
  subnoise <- levy_noise[(1+(i-1)*n.subgrid):(i*n.subgrid),]
  pred[i,] <- as.vector(t(expm(-res$Q*delta_t) %*% temp)) +
    colSums(t(expm(-sqrt(delta_t)*res$Q) %*% t(subnoise)))
  temp <- pred[i,]
}

library(viridis)
par(xpd=FALSE)
col.palette <- viridis(d)
rollback <- 50
plot((N-rollback):(N), nw_data[(N-rollback):N,1], 
     type="l", ylim=c(-5,20), xlim=c(N-rollback,N+1.25*n.horizon),
     col = col.palette[1], lwd=2,
     ylab="Value", xlab="Sample number (delta_t)")
lines((N+1):(N+n.horizon), pred[,1], 
      col = col.palette[1], lwd=2, lty=2)
lines((N+1):(N+n.horizon), nw_data_all[((N_all-n.horizon+1):N_all),1], 
      col = col.palette[1], lwd=2, lty=1)
abline(v=N)
for(j in 2:d){
  lines((N-rollback):(N), nw_data[(N-rollback):N,j], 
       type="l", ylim=c(-1,20), xlim=c(N-rollback,N),
       col = col.palette[j], lwd=2)
  lines((N+1):(N+n.horizon), nw_data_all[((N_all-n.horizon+1):N_all),j], 
        col = col.palette[j], lwd=2, lty=1)
  lines((N+1):(N+n.horizon), pred[,j], 
        col = col.palette[j], lwd=2, lty=2)
}

res$MLE_wide
L <- res$Q
library(igraph)
net<-graph_from_adjacency_matrix(L, mode = c("directed"), weighted = TRUE, add.colnames=NULL) 
#Compute number of vertecis
n_v<-gorder(net);

summary(net)

#Adding edge and vertex attributes:
s<-sample(5, size=d, replace=TRUE) 
V(net)$dealertype<-s

#Assign vertex colour corresponding to dealertype
V(net)$color <- ifelse(V(net)$dealertype %in% c(1,2,3, 4), "red", "lightblue")

#Assign vertex size corresponding to indegree
V(net)$size<-degree(net, mode="total")*6 #because 1 is a small size for a node, we multiply it by 3 to make it bigger
plot.igraph(net,layout=layout.circle, edge.width=E(net)$weight/20, edge.arrow.size=0.7, edge.label=L[L>0])
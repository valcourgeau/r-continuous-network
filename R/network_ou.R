library("Matrix")
library("foreach")
library("doParallel")

scaleMatrix <- function(x, time){
  # Allows to compute matrix times a collection of scalars
  # returns a list
  res <- list()
  n_max <- length(time)
  for(i in 1:n_max){
    res[[i]] <- x*time[i]
  }
  
  return(res)
}

collectionMatrix <- function(coll, x){
  # Allows to compute a collection of matrices by another matrix
  # returns a list
  res <- list()
  n_max <- length(time)
  for(i in 1:n_max){
    res[[i]] <- coll[[i]]*x
  }
  
  return(res)
}

nw_topo <- c(1,1,0,
             0,1,0,
             1,1,1)
nw_topo <- matrix(nw_topo, ncol = 3)
times <- c(0.1,0.4)
res <- scaleMatrix(x = nw_topo, time = times)

networkOUCov <- function(nw_topo, levy_cov, cor=FALSE, dtime=0.1){
  # This function creates NOU Covariance/Correlation matrices
  if(!identical(dim(nw_topo), dim(levy_cov)))
    stop("Wrong dimensions")

  dt <- seq(from = 0, to = 40, by = dtime)
  tmp <- Sys.time()
  seqMat <- scaleMatrix(-nw_topo, dt) 
  seqMat_exp <- mclapply(seqMat, FUN = Matrix::expm)
  seqMat_exp_t <- mclapply(t(seqMat), FUN = Matrix::expm)
  seqMat <- mclapply(1:length(dt), 
                     function(i){(seqMat_exp[[i]] %*% 
                         (levy_cov) %*% 
                         seqMat_exp_t[[i]])*dtime}) # TODO deltaT variable
  seqMat <- Reduce('+', seqMat)
  print(Sys.time() - tmp)
  return(if(cor){cov2cor(seqMat)}else{seqMat})
}

# Example networkOUCov
nw_topo <- c(1,1,0,
             1,1,0,
             0,0,1)
nw_topo <- matrix(nw_topo, ncol = 3)
nw_cov <- c(1.0,0.5,0.3,
            0.5,1.0,0.6,
            0.3,0.6,1.0)
nw_cov <- c(1.0,0.5,0.2,
            0.5,1.0,0.0,
            0.8,0.0,1.0)
nw_cov <- matrix(nw_cov, ncol = 3)
nw_cov

networkOUCov(nw_topo = nw_topo, levy_cov = nw_topo)
networkOUCov(nw_topo = nw_topo, levy_cov = nw_cov, cor = F, dtime = 0.05)
networkOUCov(nw_topo = nw_topo, levy_cov = nw_cov, cor = T, dtime = 0.05)
networkOUCov(nw_topo = matrix(rep(1, 9), ncol=3), levy_cov = nw_cov)


library(MASS)
nOUmleResiduals <- function(d, n, nw_topo, levy_cov, 
                            horizon, dtime = 0.1, empirical = F){
  # TODO Better second moment for \LL_1 including jumps
  nw_statio <- networkOUCov(nw_topo = nw_topo, 
                            levy_cov = levy_cov, dtime = dtime)
  nw_statio <- kronecker(levy_cov, solve(nw_statio))
  d.square <- dim(nw_statio)[1]
  return(mvrnorm(n = n, mu=rep(0, d.square), 
                 Sigma = nw_statio/horizon, empirical = empirical))
}

# Example
d <- 3
n <- 1000
nw_topo <- c(1,1,0,
             1,1,0,
             0,0,1)
nw_topo <- matrix(nw_topo, ncol = 3)
nw_cov <- c(1.0,0.5,0.3,
            0.5,1.0,0.6,
            0.3,0.6,1.0)
nw_cov <- matrix(nw_cov, ncol = 3)
horizon <- n * 1/24
sims <- nOUmleResiduals(d = d, n = n, nw_topo = nw_topo, levy_cov = nw_cov,
                horizon = horizon)
cor(sims)

# High-dim example
source("../../../R/network_generation.R")
d <- 10
p.link <- 0.2
n <- 1000
set.seed(42)
nw_topo <- genRdmAssymetricGraphs(d = d, p.link = p.link, theta_1=0.5/d)
eigen(nw_topo)$values
nw_cov <- genLevyCovMatrix(d = d)
horizon <- n * 1/24

sims <- nOUmleResiduals(d = d, n = n, nw_topo = nw_topo, levy_cov = nw_cov,
                        horizon = horizon)

library(mvtnorm)
nOUmleResidualsCI <- function(p, nw_topo, levy_cov, horizon, dtime = 0.1, tail = 'both'){
  # Computes the equicoordinate quantile function of the multivariate normal distribution
  # find x in R s.t. P(-x < X < x) = p
  # TODO Better second moment for \LL_1 including jumps
  nw_statio <- networkOUCov(nw_topo = nw_topo, 
                            levy_cov = levy_cov, dtime = dtime)
 
  nw_statio <- as.matrix(kronecker(solve(nw_statio), levy_cov))
  d.square <- dim(nw_statio)[1]
  return(qmvnorm(p=p, mean = rep(0, d.square), sigma = nw_statio/horizon, tail=tail))
}

# Example
d <- 3
n <- 1000
nw_topo <- c(1,1,0,
             1,1,0,
             0,0,1)
nw_topo <- matrix(nw_topo, ncol = 3)
nw_cov <- c(1.0,0.5,0.3,
            0.5,1.0,0.6,
            0.3,0.6,1.0)
nw_cov <- matrix(nw_cov, ncol = 3)
horizon <- n * 1/24
sims <- nOUmleResidualsCI(p = 0.95, nw_topo = nw_topo, levy_cov = nw_cov,
                        horizon = horizon)

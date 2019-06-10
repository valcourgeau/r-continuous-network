# This script attempts to provide construction functions for 4 main types of networks;
# They cover the following properties
#   - max(degree) = 0, 2, 4/8, d-1

#' @example IsolatedNetwork(d=10)
IsolatedNetwork <- function(d, theta_1, theta_2=1){
  # theta_2 is the diagonal element
  return(theta_2*diag(d))
}

IsolatedNetwork(d=10, 2)

#' @example PolymerNetwork(d=10)
PolymerNetwork <- function(d, theta_1, theta_2=1, directed=F){
  assertthat::assert_that(theta_2 > abs(theta_1))
  # theta_2 is the diagonal element
  mat_temp <-# IsolatedNetwork(d, theta_2 = theta_2)
  mat_temp <- AugmentedDiag(d = d, offset = 1) + if(directed){-AugmentedDiag(d = d, offset = 1)} else{AugmentedDiag(d=d,offset=-1)}
  mat_temp <- theta_1 * mat_temp / rowSums(abs(mat_temp))
  return(mat_temp + IsolatedNetwork(d = d, theta_2 = theta_2))
}

PolymerNetwork(d=10, theta_1 = 2, theta_2 = 5)


CheckThetas <- function(theta_1, theta_2){
  assertthat::assert_that(theta_2 > abs(theta_1), msg = paste('In', match.call()[1], 'must verify theta_2 > abs(theta_1)'))
}

#' @example LatticeNetwork(10, theta_1 = 1, theta_2 = 2)
LatticeNetwork <- function(d, theta_1, theta_2, diag.connections=F, directed=F){
  # theta_2 diagonal element
  # TODO add diag.connections and directed for it
  
  CheckThetas(theta_1, theta_2)
  
  if(diag.connections) stop('NotImplementedError')
  if(directed) stop('NotImplementedError')
  
  # net_matrix <- matrix(0, ncol = d, nrow = d)
  # net_temp <- diag(d-1)
  # net_temp <- cbind(rep(0, d-1), net_temp)
  # net_temp <- rbind(net_temp, rep(0, d))
  # net_matrix <- net_matrix + net_temp
  # if(directed){
  #   net_matrix <- net_matrix - t(net_temp)
  # }else{
  #   net_matrix <- net_matrix + t(net_temp)
  # }
  
  #net_matrix <- ngspatial::adjacency.matrix(d)
  d_real <- round(sqrt(d))
  net_matrix <- copCAR::adjacency.matrix(d_real)
  if(d_real^2 < d){
    net_matrix <- rbind(net_matrix, matrix(0, nrow=d-d_real^2, ncol=d_real*d_real))
    net_matrix <- cbind(net_matrix, rbind(matrix(0, ncol=d-d_real^2, nrow=d_real*d_real), diag(1, nrow = d - d_real^2)))
  }
 
  net_matrix <- net_matrix / rowSums(abs(net_matrix))
  net_matrix <- theta_1 * net_matrix
  diag(net_matrix) <- theta_2
 
  return(net_matrix)
}

LatticeNetwork(10, theta_1 = 1, theta_2 = 2) %>% rowSums

#' @example FullyConnectedNetwork(10, 1, 2)
FullyConnectedNetwork <- function(d, theta_1, theta_2, directed=F){
  net_temp <- matrix(1, d, d)
  diag(net_temp) <- 0
  net_temp <- theta_1 * net_temp / rowSums(net_temp)
  
  if(directed){
    net_temp[lower.tri(net_temp)] <- -net_temp[lower.tri(net_temp)]
  }
  
  diag(net_temp) <- theta_2
  return(net_temp)
}

FullyConnectedNetwork(10, 1, 2, F)


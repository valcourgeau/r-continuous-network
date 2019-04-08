# This script attempts to provide construction functions for 4 main types of networks;
# They cover the following properties
#   - max(degree) = 0, 2, 4/8, d-1

#' @example PolymerNetwork(d=10)
PolymerNetwork <- function(d, theta_2=1){
  # theta_2 is the diagonal element
  return(theta_2*diag(d))
}

PolymerNetwork(d=10, 2)


CheckThetas <- function(theta_1, theta_2){
  assertthat::assert_that(theta_2 > abs(theta_1), msg = paste('In', match.call()[1], 'must verify theta_2 > abs(theta_1)'))
}

#' @example LatticeNetwork(10, theta_1 = 1, theta_2 = 2)
LatticeNetwork <- function(d, theta_1, theta_2, diag.connections=F, directed=F){
  # theta_2 diagonal element
  # TODO add diag.connections and directed for it
  
  CheckThetas(theta_1, theta_2)
  
  if(diag.connections) stop('NotImplementedError')
  
  net_matrix <- matrix(0, ncol = d, nrow = d)
  net_temp <- diag(d-1)
  net_temp <- cbind(rep(0, d-1), net_temp)
  net_temp <- rbind(net_temp, rep(0, d))
  net_matrix <- net_matrix + net_temp
  if(directed){
    net_matrix <- net_matrix - t(net_temp)
  }else{
    net_matrix <- net_matrix + t(net_temp)
  }
  
  net_matrix <- net_matrix / rowSums(abs(net_matrix))
  net_matrix <- theta_1 * net_matrix
  diag(net_matrix) <- theta_2
 
  return(net_matrix)
}

LatticeNetwork(10, theta_1 = 1, theta_2 = 2)

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


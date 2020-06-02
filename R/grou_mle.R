
NodeMLE <- function(times, data, output="vector"){
  assertthat::assert_that(
    output %in% c("vector", "matrix"),
    msg=paste('output should be "node" or "network", given', output)
  )
  assertthat::are_equal(length(times), nrow(data))
  
  n_nodes <- ncol(data)
  n_data <- nrow(data)
  delta_t <- diff(times)
  delta_data <- apply(data, MARGIN = 2, diff)
  wo_last <- data[-n_data,]
  expanded_delta_data <- t(rep(1, n_nodes)) %x% delta_data
  expanded_wo_last <- wo_last %x% t(rep(1, n_nodes))
  
  numerator <- colSums(expanded_wo_last * expanded_delta_data)
  denominator <- t(wo_last * delta_t) %*% wo_last
  inv_denominator <- solve(denominator)
  mle <- -(inv_denominator %x% diag(n_nodes)) %*% numerator
  if(output == "vector"){
    return(mle)
  }else{
    return(matrix(mle, nrow = n_nodes, ncol = n_nodes, byrow = FALSE))
  }
}

CoreNodeMLE <- function(times, data){
  n_nodes <- ncol(data)
  n_data <- nrow(data)
  delta_t <- diff(times)
  delta_data <- apply(data, MARGIN = 2, diff)
  wo_last <- data[-n_data,]
  expanded_delta_data <- t(rep(1, n_nodes)) %x% delta_data
  expanded_wo_last <- wo_last %x% t(rep(1, n_nodes))
  
  numerator <- colSums(expanded_wo_last * expanded_delta_data)
  denominator <- t(wo_last * delta_t) %*% wo_last
  return(list(numerator=numerator, denominator=denominator))
}

NodeMLELong <- function(times, data, div=1e5, output="vector"){
  assertthat::assert_that(
    output %in% c("vector", "matrix"),
    msg=paste('output should be "node" or "network", given', output)
  )
  assertthat::are_equal(length(times), nrow(data))
  
  n_nodes <- ncol(data)
  n_data <- nrow(data)
  
  idx <- (0:(n_data %/% div)) * div
  if(tail(idx, n=1) != n_data){ # check if last elem is in
    idx <- c(idx, n_data)
  } 
  
  collection_num_denom <- lapply(
    1:(length(idx)-1),
    FUN = function(i){
      idx_couple <- idx[i]:idx[i+1]
      CoreNodeMLE(times[idx_couple], data[idx_couple,])
    }
  )
  numerator <- Reduce('+', lapply(collection_num_denom, function(coll){coll$numerator}))
  inv_denominator <- solve(Reduce('+', lapply(collection_num_denom, function(coll){coll$denominator})))
  mle <- - (inv_denominator %x% diag(n_nodes)) %*% numerator
  
  if(output == "vector"){
    return(mle)
  }else{
    return(matrix(mle, nrow = n_nodes, ncol = n_nodes, byrow = FALSE))
  }
}


GrouMLE <- function(times, data, adj=NA, mode = "node", output = "vector"){
  assertthat::assert_that(
    mode %in% c("node", "network"),
    msg=paste('mode should be "node" or "network", given', mode)
  )
  
  # node MLE without adj
  
  if(nrow(data) < 1e6){
    node_mle <- NodeMLE(times, data, output = "vector")
  }else{
    node_mle <- NodeMLELong(times, data, output = "vector")
  }
  
  if(any(is.na(adj))){
    return(node_mle)
  }
  n_nodes <- ncol(adj)
  adj_norm <- RowNormalised(adj)
  
  if(mode == "node"){
    if(output == "vector"){
      d_a <- c(diag(n_nodes)+adj_norm)
      return(d_a*node_mle)
    }else{
      if(output == "matrix"){
        adj_full <- diag(n_nodes)+adj_norm
        return(adj_full * matrix(node_mle, n_nodes, n_nodes, byrow = F))
      }
    }
  }else{
    a_l2 <- sqrt(sum(adj_norm^2))
    if(a_l2 < .Machine$double.eps){
      theta_1 <- 0.0
    }else{
      theta_1 <- (c(adj_norm) %*% node_mle) / a_l2
    }
    theta_2 <- (c(diag(n_nodes)) %*% node_mle)[1,1] / n_nodes
    if(output == "vector"){
      return(c(theta_1, theta_2))
    }else{
      if(output == "matrix"){
        return(theta_1*adj_norm + theta_2*diag(n_nodes))
      }
    }
  }
}

n_nodes <- 5
n_sample <- 10000000
sample_path <- ConstructPath(diag(n_nodes), matrix(rnorm(n_sample*n_nodes, 0, 1), ncol=n_nodes), n_nodes, 0.01)
# NodeMLELong(times=seq(0, by=0.01, length.out = n_sample), data=sample_path, output = "matrix")
GrouMLE(times=seq(0, by=0.01, length.out = n_sample), data=sample_path, adj = diag(n_nodes), mode="node", output = "matrix")

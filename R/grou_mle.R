
NodeMLE <- function(times, data, output="vector"){
  assertthat::assert_that(
    output %in% c("vector", "matrix"),
    msg=paste('output should be "node" or "network", given', output)
  )
  assertthat::are_equal(length(times), nrow(data))
  
  n_nodes <- ncol(data)
  delta_t <- diff(times)
  delta_data <- apply(data, MARGIN = 2, diff)
  wo_last <- data[-nrow(data),]
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

GrouMLE <- function(times, data, adj=NA, mode = "node", output = "vector"){
  assertthat::assert_that(
    mode %in% c("node", "network"),
    msg=paste('mode should be "node" or "network", given', mode)
  )
  
  # node MLE without adj
  node_mle <- NodeMLE(times, data, output = "vector")
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
    if(a_l2 == 0.0){
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

n_nodes <- 10
n_sample <- 100000
sample_path <- ConstructPath(diag(n_nodes), matrix(rnorm(n_sample*n_nodes, 0, 1), ncol=n_nodes), n_nodes, 0.01)
NodeMLE(times=seq(0, by=0.01, length.out = n_sample), data=sample_path, output = "matrix")
GrouMLE(times=seq(0, by=0.01, length.out = n_sample), data=sample_path, adj = diag(n_nodes), mode="node", output = "matrix")

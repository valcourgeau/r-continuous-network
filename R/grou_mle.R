Rcpp::cppFunction(depends = "RcppArmadillo", code = 'arma::mat mat_inv(const arma::mat &Am) {return inv(Am);}')

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
  sd_denom <- sd(denominator)
  inv_denominator <- solve(denominator / sd_denom)
  inv_denominator <- inv_denominator / sd_denom
  mle <- -(inv_denominator %x% diag(n_nodes)) %*% numerator
  if(output == "vector"){
    return(mle)
  }else{
    return(matrix(mle, nrow = n_nodes, ncol = n_nodes, byrow = FALSE))
  }
}

CoreNodeMLE <- function(times, data, thresholds=NA){
  # Computes the numerator and denominator of a batch of data
  n_nodes <- ncol(data)
  n_data <- nrow(data)
  delta_t <- diff(times)
  delta_data <- apply(data, MARGIN = 2, diff)
  wo_last <- data[-n_data,]
 
  kron_delta_data <- lapply(
    1:ncol(wo_last), 
    function(i){
      tmp <- t(delta_data*as.vector(wo_last[,i])) # transpose for the filtering
      if(!any(is.na(thresholds))){
        filters <- apply(abs(delta_data), 1, '<=', rep(thresholds[i], n_nodes))
        filters <- t(filters)
        tmp <- tmp * as.vector(filters)
      }
      # tmp <- delta_data*as.vector(wo_last[,i])
      
      return(colSums(t(tmp)))
    })
  numerator <- t(do.call(cbind, kron_delta_data))
  denominator <- t(wo_last * as.vector(delta_t)) %*% wo_last
  
  return(list(numerator=numerator, denominator=denominator))
}

NodeMLELong <- function(times, data, thresholds, div=1e5, output="vector"){
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
  idx <- idx[which(idx < n_data)]
  collection_num_denom <- lapply(
    1:(length(idx)-1),
    FUN = function(i){
      idx_couple <- idx[i]:idx[i+1]
      CoreNodeMLE(times = times[idx_couple], data = data[idx_couple,], thresholds = thresholds)
    }
  )
  
  numerator <- Reduce('+', lapply(collection_num_denom, function(coll){coll$numerator}))
  sd_numerator <- sd(numerator)
  numerator <- numerator / sd_numerator
  denominator <- Reduce('+', lapply(collection_num_denom, function(coll){coll$denominator}))
  sd_denominator <- sd(denominator)
  
  inv_denominator <- solve(denominator / sd_denominator)
  # inv_denominator <- mat_inv(denominator / sd_denominator)
  inv_denominator <- inv_denominator / sd_denominator
  numerator <- matrix(numerator, nrow = n_nodes, byrow = FALSE)
  
  mle <- apply(numerator, MARGIN = 2, function(x){- inv_denominator %*% x}) * sd_numerator
  if(output == "vector"){
    return(c(mle))
  }else{
    return(mle)
  }
}

GrouMLE <- function(times, data, adj=NA, thresholds=NA, div = 1e3, mode = "node", output = "vector"){
  assertthat::assert_that(
    mode %in% c("node", "network"),
    msg=paste('mode should be "node" or "network", given', mode)
  )
  if(length(times)!=nrow(data)) stop('length(times)!=nrow(data)')
  if(all(is.na(adj))){return(node_mle)}
  
  n_nodes <- ncol(adj)
  adj_normalised <- RowNormalised(adj)
  diag(adj_normalised) <- 0.0 # to be sure
  
  node_mle <- NodeMLELong(times, data, thresholds = thresholds, div = div, output = "vector")
  
  if(mode == "node"){
    adj_normalised[adj_normalised!=0] <- 1
    if(output == "vector"){
      d_a <- c(diag(n_nodes)+adj_normalised)
      return(d_a*node_mle)
    }else{
      if(output == "matrix"){
        adj_full <- diag(n_nodes)+adj_normalised
        return(adj_full * matrix(node_mle, n_nodes, n_nodes, byrow = F))
      }
    }
  }else{
    a_l2 <- sum(adj_normalised)
    adj_ones <- adj_normalised
    adj_ones[adj_ones!=0] <- 1
    diag(adj_ones) <- 0
    
    if(a_l2 < .Machine$double.eps){
      theta_1 <- 0.0
    }else{
      theta_1 <- (as.vector(adj_ones) %*% node_mle) / a_l2
    }
    
    theta_2 <- (c(diag(n_nodes)) %*% node_mle) / n_nodes
    if(output == "vector"){
      return(c(theta_1, theta_2))
    }else{
      if(output == "matrix"){
        return(c(theta_1)*adj_normalised + c(theta_2)*diag(n_nodes))
      }
    }
  }
}

FasenRegression <- function(data){
  data <- as.matrix(data)
  data_without_first <- data[-1,]
  data_without_last <- data[-nrow(data),]
  sub <- t(data_without_last) %*% data_without_last
  return(t(data_without_first) %*% data_without_last %*% solve(sub))
}

n_nodes <- 50
n_sample <- 50000
set.seed(42)
adj_test <- diag(n_nodes)
adj_test[2,1] <- 0.5

mesh_size <- 0.01
sample_path <- ConstructPath(adj_test, matrix(rnorm(n_sample*n_nodes, 0, 1*mesh_size^(1/2)), ncol=n_nodes), rep(0, n_nodes), mesh_size)
GrouMLE(times=seq(0, by=mesh_size, length.out = n_sample), data=sample_path, adj = adj_test, div = 1e3, mode="node", output = "matrix")
adj_test
GrouMLE(times=seq(0, by=mesh_size, length.out = n_sample), data=sample_path, adj = adj_test, div = 1e3, mode="network", output = "vector")
GrouMLE(times=seq(0, by=mesh_size, length.out = n_sample), data=sample_path, adj = adj_test, div = 1e3, mode="network", output = "matrix")

adj_test <- diag(n_nodes)
adj_test[2,1] <- 0.4
adj_test[1,2] <- 0.2
adj_test[1,3] <- 0.2

mesh_size <- 0.01
sample_path <- ConstructPath(adj_test, matrix(rnorm(n_sample*n_nodes, 0, 1*mesh_size^(1/2)), ncol=n_nodes), rep(0, n_nodes), mesh_size)
GrouMLE(times=seq(0, by=mesh_size, length.out = n_sample), data=sample_path, adj = adj_test, div = 1e3, mode="node", output = "matrix")
GrouMLE(times=seq(0, by=mesh_size, length.out = n_sample), data=sample_path, adj = adj_test, div = 1e3, mode="network", output = "vector")



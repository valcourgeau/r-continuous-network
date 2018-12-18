### Function to assess goodness-of-fit of MLE estimations
# based on Renyi entropy and statistics (Morales et al. (2000))

NOUGof <- function(theta, theta_true, data, times, quiet=F){
  N <- nrow(data)
  d <- ncol(data)
  diff_times <- TimeMatrix(times, d)
  
  theta <- as.vector(theta)
  theta_true <- as.vector(theta_true)
  
  # trying to do kronecker product with d*d matrices
  k_n <- matrix(colSums(rep.col(data[-N,], n_copy = d) * 
                          rep.mat(data[-N,], n_copy = d) * 
                          rep.col(diff_times[-N,], d)), d)
  theta_mat <- matrix(theta-theta_true, nrow = d, byrow = T)
  rv_value <- as.vector(t(theta_mat %*% k_n)) %*%(theta-theta_true)
  
  eigen_vals <- eigen(k_n, only.values = T)$value
  first <- imhof(q=rv_value, lambda = eigen_vals)
  abserr <- first$abserr
  first <- 1-first$Qq
  last <- imhof(q=rv_value, lambda = eigen_vals)$Qq
  
  probability <- abs(last-first)

  if(probability > 0.975){
    probability <- '***'
  }else{
    if(probability > 0.95){
      probability <- '**'
    }else{
      if(probability > 0.90){
        probability <- '*'
      }
    }
  }
  test_values <- list("lower"=first,
                     "upper"=last,
                     "97.5%='***' / 95%='**' / 90%='*'"=probability,
                     "Abs.error"=abserr)
  if(!quiet){
    cat("Goodness-of-test based on Renyi statistics\n")
  }
  return(test_values)
}

# https://cran.r-project.org/web/packages/CompQuadForm/CompQuadForm.pdf
ans <- NOUGof(theta = mle_matrix_test, theta_true = nw_q, data=nw_data, times =times, quiet = T)
ans

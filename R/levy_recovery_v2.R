
LevyRecovery <- function(fitted_adj, data, times, look_ahead=1){
  # only for OU processes; t_{n} -> t_{n+look_ahead}
  assertthat::are_equal(length(times), nrow(data))
  
  diff_x <- apply(data, MARGIN = 2, FUN = function(x){diff(x, lag = look_ahead)})
  diff_times <- diff(times, lag =  1)
  integrated_x <- zoo::rollapply(cbind(data, times),
                                 width = look_ahead+1,
                                 by.column = FALSE,
                                 FUN = function(sub_integ){
                                   tt <- sub_integ[,ncol(sub_integ)]
                                   y_col_bind <- sub_integ[,-ncol(sub_integ)]
                                   return(apply(y_col_bind,
                                                MARGIN = 2,
                                                FUN = function(y){
                                                  pracma::trapz(x = tt,
                                                                y = y)}
                                   ))}
  ) # create integral component-wise
  
  # few dimension checks
  assertthat::are_equal(nrow(integrated_x), nrow(diff_x))
  assertthat::are_equal(nrow(data)-look_ahead, nrow(diff_x))
  assertthat::are_equal(dim(fitted_adj)[1], ncol(integrated_x))
  
  q_integrated_x <- apply(integrated_x, MARGIN = 1, function(x){return(fitted_adj %*% x)})
  if(any(vapply(class(fitted_adj), function(x){x == 'dgCMatrix'}, T))){
    q_integrated_x <- lapply(q_integrated_x, as.matrix) # getting rid of potential sparse matrix
    q_integrated_x <- t(as.matrix(as.data.frame(q_integrated_x)))
  }else{
    q_integrated_x <- t(q_integrated_x)
  }
  
  cat('q_integrated_x', dim(q_integrated_x), '\n')
  
  # TODO implement ghyp
  # TODO implement with m > 1
  assertthat::are_equal(dim(diff_x)[1], dim(q_integrated_x)[1])
  assertthat::are_equal(dim(diff_x)[2], dim(q_integrated_x)[2])
  
  recover <- diff_x + q_integrated_x
  
  recovery_results <- list()
  recovery_results[['increments']] <- recover
  recovery_results[['cumsum']] <- apply(recover, MARGIN = 2, FUN = cumsum)
  recovery_results[['fitted_adj']] <- fitted_adj
  return(recovery_results)
}

FitLevyRecoveryDiffusion <- function(data, look_ahead=1){
  # m how much we hop at every increments
  # data is increments
  # ghyp_models <- c(ghyp::fit.NIGmv, ghyp::fit.gaussmv, ghyp::fit.VGmv, ghyp::fit.tmv, ghyp::fit.ghypmv)
  ghyp_models <- c(ghyp::fit.ghypmv)
  ghyp_name <- c('NIG', 'GAUSS', 'VG', 'T', 'FULL')
  ghyp_name <- c('FULL')
  i <- 1
  res <- list()
  for(model in ghyp_models){
    # res[[ghyp_name[i]]] <- model(data = data[which(1:nrow(data) %% (look_ahead+1) == 1),])
    res[[ghyp_name[i]]] <- model(data = data)
    i <- i + 1
  }
  return(res)
}

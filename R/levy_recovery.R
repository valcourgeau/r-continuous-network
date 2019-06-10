# Script to implement levy increment recovery given network topology or
CHECK_EXAMPLE <- FALSE
install.packages('pracma')

LevyRecovery <- function(nw_topo, data, times, m=1, fitted=FALSE, on_matrix=FALSE){
   # only for OU processes
  # m = number of look-ahead points, default is one to do t_{n} -> t_{n+1}
  assertthat::are_equal(length(times), nrow(data))

  if(!fitted){
    if(on_matrix){
      Q_hat <- NOUmatrix(nw_topo = nw_topo, data = data, times = times, thresholds = rep(1000, ncol(data)))
    }else{
      Q_hat <- NOUfit(nw_topo = nw_topo, data = data, times = times, thresholds = rep(1000, ncol(data)))
      Q_hat <- Q_hat[['Q']]
    }
  }else{
    Q_hat <- nw_topo
    warning('Given fitted Q, ignoring flag on_matrix;')
  }
  
  diff_x <- apply(data, MARGIN = 2, FUN = function(x){diff(x, lag = m)})
  diff_times <- diff(times, lag =  1)
  # diff_times_2 <- diff(times, lag = 2)
  # diff_times <- c(diff_times[1], diff_times_2, diff_times[length(diff_times)])
  
  # integrated_x <- zoo::rollapply(t(t(data[-nrow(data),])*diff_times), width = m, FUN = sum) # create integral component-wise
  integrated_x <- zoo::rollapply(cbind(data, times),
                                 width = m+1,
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
  assertthat::are_equal(nrow(data)-m, nrow(diff_x))
  assertthat::are_equal(dim(Q_hat)[1], ncol(integrated_x))

  q_integrated_x <- apply(integrated_x, MARGIN = 1, function(x){return(Q_hat %*% x)})
  if(class(Q_hat) == 'dgCMatrix'){
    q_integrated_x <- lapply(q_integrated_x, as.matrix) # getting rid of potential sparse matrix
    q_integrated_x <- t(as.matrix(as.data.frame(q_integrated_x)))
  }else{
    q_integrated_x <- t(q_integrated_x)
  }

  # TODO implement ghyp
  # TODO implement with m > 1
  assertthat::are_equal(dim(diff_x)[1], dim(q_integrated_x)[1])
  assertthat::are_equal(dim(diff_x)[2], dim(q_integrated_x)[2])
  
  recover <- diff_x + q_integrated_x

  recovery_results <- list()
  recovery_results[['increments']] <- recover
  recovery_results[['cumsum']] <- apply(recover, MARGIN = 2, FUN = cumsum)
  recovery_results[['m']] <- m
  return(recovery_results)
}

FitLevyRecoveryDiffusion <- function(data, m=1){
  # m how much we hop at every increments
  # data is increments
  ghyp_models <- c(ghyp::fit.NIGmv, ghyp::fit.gaussmv, ghyp::fit.VGmv, ghyp::fit.tmv)
  ghyp_name <- c('NIG', 'GAUSS', 'VG', 'T')
  i <- 1
  res <- list()
  for(model in ghyp_models){
    res[[ghyp_name[i]]] <- model(data = data[which(1:nrow(data) %% (m+1) == 1),])
    i <- i + 1
  }
  return(res)
  
}

# after doing re_europe_inference.R
data_tmp <- as.matrix(df_load[,-1])
data_tmp <- (rowSums(data_tmp) * data_tmp)
adj_tmp <- igraph::as_adjacency_matrix(topo_graph)

par(mfrow=c(2,5))
layout(matrix(c(1,3,5,7,9,2,4,6,8,10), 2, 5, byrow = TRUE))
for(m in 1:5){
  lr <- LevyRecovery(nw_topo = adj_tmp, data = data_tmp, times = (1:nrow(data_tmp))/24, m = m, fitted = F, on_matrix = T)
  if(m == 1){
    acf(lr$increments[,15], main=paste('Node 15 ACF:', m,'-step'), 
        cex.axis=1.3, cex.lab = 1.3, cex.main = 1.3)
    acf(lr$increments[,35], main=paste('Node 35 ACF:', m,'-step'), 
        cex.axis=1.3, cex.lab = 1.3, cex.main = 1.3)
  }else{
    acf(lr$increments[which(1:nrow(lr$increments) %% (m) == 1),15], main=paste('Node 15 ACF:', m,'-step'), 
        cex.axis=1.3, cex.lab = 1.3, cex.main = 1.3)
    acf(lr$increments[which(1:nrow(lr$increments) %% (m) == 1),35], main=paste('Node 35 ACF:', m,'-step'), 
        cex.axis=1.3, cex.lab = 1.3, cex.main = 1.3)
  }
  
}


LevyRecovery(nw_topo = adj_tmp, data = data_tmp, times = (1:nrow(data_tmp))/24, m = 1, fitted = F, on_matrix = F)

lr <- LevyRecovery(nw_topo = adj_tmp, data = data_tmp, times = (1:nrow(data_tmp))/24, m = 1, fitted = F, on_matrix = F)

ind_seq <- seq(1, to=nrow(lr$increments), length.out = 5) %>% round

library(RColorBrewer)
par(mfrow=c(2,4))
# layout(matrix(c(1,3,5,7,9,2,4,6,8,10), 2, 5, byrow = TRUE))
n_numbers <- c(15, 35)
for(node_number in n_numbers){
  marker <- viridis(5+length(ind_seq))[5:(length(ind_seq)+5)]
  print(marker)
  hist(lr$increments[which(1:ind_seq[2] %% (lr$m) == 1), node_number], breaks = 20, xlim=c(-.8,.8), probability=T,
       col=marker[1],
       main=paste('Node', node_number, '2-step\nsemester', 1), cex.main = 1.3,
       xlab='Increment value', cex.lab = 1.3)
  for(j in 2:(length(ind_seq) - 1)){
    hist(lr$increments[ind_seq[j]-1+which(ind_seq[j]:ind_seq[j+1] %% (lr$m) == 1), node_number], breaks = 20, probability=T,
         col=marker[j], xlim=c(-.8,.8),
         main=paste('Node', node_number, '2-step\nsemester', j), cex.main = 1.3,
         xlab='Increment value', cex.lab = 1.3)
  }
  
  # plot(stats::density(lr$increments[which(1:ind_seq[2] %% (lr$m) == 1), node_number]), 
  #      col =rgb(red=col2rgb(marker$color[1])[1,1]/255,green=col2rgb(marker$color[1])[2,1]/255, 
  #               blue = col2rgb(marker$color[1])[3,1]/255, alpha=0.9), lwd=2, xlim=c(-1,1),
  #      main = paste('Node', node_number, '2-step increment density'), cex.main = 1.3,
  #      xlab='Increment value', cex.lab = 1.3)
  # for(j in 1:(length(ind_seq) - 1)){
  #   print(j)
  #   lines(stats::density(lr$increments[which(ind_seq[j]:ind_seq[j+1] %% lr$m == 1), node_number]), 
  #         col=marker$color[j], lwd=2)
  #   
  # }
}

# lr <- LevyRecovery(nw_topo = adj_tmp, data = data_tmp, times = (1:nrow(data_tmp))/24, m = 2, fitted = F, on_matrix = F)

fit_ghyp <- FitLevyRecoveryDiffusion(data = lr$increments, m=lr$m)

fit_ghyp['NIG']$NIG@sigma %>% heatmap
cov2cor(fit_ghyp['NIG']$NIG@sigma) %>% 
  (function(x){heatmap.2(x, Rowv = F, Colv = F, col = viridis, scale = 'none', trace='none', 
                         labRow = vapply(1:50, FUN = function(i){if(i<25){paste('P-',stringr::str_pad(string = i, width = 2, pad = '0'),sep = '')}else{
                           paste('E-',stringr::str_pad(string = i-24, width = 2, pad = '0'),sep = '')
                         }}, FUN.VALUE = ''),
                         labCol = vapply(1:50, FUN = function(i){if(i<25){paste('P-',stringr::str_pad(string = i, width = 2, pad = '0'),sep = '')}else{
                           paste('E-',stringr::str_pad(string = i-24, width = 2, pad = '0'),sep = '')
                         }}, FUN.VALUE = ''),
                         colsep = 24, rowsep = 24,
                         lhei = c(0.1,1,0.2),
                         lwid = c(0.1,1,0.1),
                         cexRow = 1.3, cexCol = 1.3,
                         lmat=matrix(c(0,0,0, 0,1,2,0, 4,3), byrow = T, ncol=3), 
                         key.title = '', key.par = list(cex.lab=1.5),
                         key.ylab = '')})

cov2cor(ghyp_full@sigma) %>% 
  (function(x){heatmap.2(x, Rowv = F, Colv = F, col = viridis, scale = 'none', trace='none', 
                         labRow = vapply(1:50, FUN = function(i){if(i<25){paste('P-',stringr::str_pad(string = i, width = 2, pad = '0'),sep = '')}else{
                           paste('E-',stringr::str_pad(string = i-24, width = 2, pad = '0'),sep = '')
                         }}, FUN.VALUE = ''),
                         labCol = vapply(1:50, FUN = function(i){if(i<25){paste('P-',stringr::str_pad(string = i, width = 2, pad = '0'),sep = '')}else{
                           paste('E-',stringr::str_pad(string = i-24, width = 2, pad = '0'),sep = '')
                         }}, FUN.VALUE = ''),
                         colsep = 24, rowsep = 24,
                         lhei = c(0.1,1,0.2),
                         lwid = c(0.1,1,0.1),
                         cexRow = 1.3, cexCol = 1.3,
                         lmat=matrix(c(0,0,0, 0,1,2,0, 4,3), byrow = T, ncol=3), 
                         key.title = '', key.par = list(cex.lab=1.5),
                         key.ylab = '')})

fit_ghyp['NIG']$NIG@aic
fit_ghyp['GAUSS']$GAUSS@aic
fit_ghyp['VG']$VG@aic
fit_ghyp['T']$T@aic

fit_ghyp['NIG']$NIG@llh
fit_ghyp['GAUSS']$GAUSS@llh
fit_ghyp['VG']$VG@llh
fit_ghyp['T']$T@llh

rlist::list.save(fit_ghyp, '50_nodes_ghyp_fit_v3.RData')
cov2cor(fit_ghyp['VG']$VG@sigma) %>% (function(x){heatmap.2(x, Rowv = F, Colv = F)})

Q_used <- NOUfit(nw_topo = adj_tmp, data = data_tmp, times = (1:nrow(data_tmp))/24, thresholds = rep(1000, ncol(data_tmp)))
ghyp::rghyp(1000, fit_ghyp['VG']$VG)[,1] %>% plot


if(CHECK_EXAMPLE){
  source("R/network_generation.R")
  d <- 5 #dims
  N <- 24*1500 #numbers of points
  M <- 10 # number of simulations
  Y0 <- 1 #start point
  delta_t <- 1/24
  draw_n <- seq(from=8, to=N, by=1000)
  mle_fit_example <- rep(0, length(draw_n)*M*2)
  mle_fit_example <- array(data = mle_fit_example,
                           dim=c(M,length(draw_n),2))
  
  sigma <- 0.1
  set.seed(42)
  nw_topo <- genRdmAssymetricGraphs(d = d, p.link = 0.25,
                                    theta_1 = 1, theta_2 = 0)
  nw_topo <- StdTopo(nw_topo)
  
  times <- seq(from = 0, by = delta_t, length.out = N)
  thresholds <- rep(0, d)
  for(i in 1:M){
    nw_data_bm <- matrix(rnorm(n = d*N, mean = 0, sd = sigma*sqrt(delta_t)), ncol = d)
    nw_data <- matrix(0, ncol=d, nrow=N)
    nw_data[1,] <- nw_data_bm[1,]*sqrt(delta_t)
    
    nw_q <- 0.1*nw_topo
    diag(nw_q) <- 1
    
    for(index in 2:N){
      nw_data[index,] <-  nw_data[index-1,] - (nw_q %*% nw_data[index-1,]) * (times[index]-times[index-1]) + nw_data_bm[index,]
    }
    
    #plot(nw_data[1:300,1], type="l")
    
    j <- 0
    for(n_temp in draw_n){
      res <- NOUfit(nw_topo = nw_topo, times = times[1:n_temp],
                    data = nw_data[1:n_temp,], thresholds = thresholds)
      mle_fit_example[i,j,] <- res$MLE_wide
      
      j <- j + 1
    }
  }
  mle_true <- c(0.1, 1.0)
  par(mfrow=(c(1,2)))
  index_mle_param <- 1
  
  plot(draw_n[-length(draw_n)], colMeans(mle_fit_example[,-length(draw_n),index_mle_param])
       - mle_true[index_mle_param], main="First para MLE / Error and Empirical var")
  abline(h=0)
  lines(draw_n[-length(draw_n)], colMeans((mle_fit_example[,-length(draw_n),index_mle_param]
                                           -colMeans(mle_fit_example[,-length(draw_n),index_mle_param]))^2),
        col="red", lty=2, lwd=2)
  
  index_mle_param <- 2
  
  plot(draw_n[-length(draw_n)], colMeans(mle_fit_example[,-length(draw_n),index_mle_param])
       - mle_true[index_mle_param], main="First para MLE / Error and Var", ylim=c(-0.01,0.05))
  abline(h=0)
  lines(draw_n[-length(draw_n)], colMeans((mle_fit_example[,-length(draw_n),index_mle_param]
                                           -colMeans(mle_fit_example[,-length(draw_n),index_mle_param]))^2),
        col="red", lty=2, lwd=2)
}

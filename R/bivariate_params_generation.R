# This script attmepts to provide the functions necessary to generate pairs of parameters 

BivariateGrid <- function(interval_1, interval_2, n_1, n_2){
  assertthat::assert_that(interval_2[1] >= 0.0, msg=paste('In', match.call()[1], 'interval2 (for diagonal element) must have non-negative elements'))
  assertthat::assert_that(length(interval_1) == 2)
  assertthat::assert_that(length(interval_2) == 2)
  assertthat::assert_that(diff(interval_2) >= 0)
  assertthat::assert_that(n_1 > 0)
  assertthat::assert_that(n_2 > 0)
  
  if(diff(interval_2) == 0){
    seq_1 <- seq(interval_1[1], interval_1[2], length.out = n_1)
      return(cbind(seq_1, rep(interval_2[1], n_1)))
  }else{
    #res <- matrix(0, ncol = 2, nrow = n_1 * n_2)
    # res <- , 
    #                     rep(interval_2[1], interval_2[2], length.out = n_1))
    # TODO HERE
    
    seq_1 <- seq(interval_1[1], interval_1[2], length.out = n_1)
    seq_2 <- seq(interval_2[1], interval_2[2], length.out = n_2)
    return(expand.grid(x=seq_1, y=seq_2))
  }
}


BivariateGrid(c(1,2), c(10,11), 3, 3) %>% is.data.frame()
plot(BivariateGrid(c(0,1),c(2,5), 10, 5))

CornerGrid <- function(top_x, top_y, n_x, n_y, scale='exp'){
  top_x <- abs(top_x)
  top_y <- abs(top_y)
  
  seq_y <- seq(0, top_y, length.out = n_y)
  grid_pts <- matrix(0, ncol = 2, nrow = n_x * n_y)
  
  for(i in 2:n_y){
    if(scale == 'exp'){
      tmp_x <- seq_y[i]-exp(seq(log(1e-5), log(seq_y[i]), length.out = round(n_x/2)))
      if(n_x %% 2 == 0){
        tmp_x <- c(-tmp_x, tmp_x)
      }else{
        tmp_x <- c(-tmp_x, 0, tmp_x)
      }
    }else{
      if(scale == 'linear'){
        tmp_x <- seq(-seq_y[i], seq_y[i], length.out = round(n_x))
      }
    }
    
    grid_pts[((i-1)*n_x+1):(i*n_x),] <- cbind(
        tmp_x, 
        rep(seq_y[i], n_x)
      )
  }
  return(data.frame(x=grid_pts[,1], y=grid_pts[,2]))
}

CornerGrid(top_x = 2, top_y = 3, n_x = 50, n_y = 50) %>% plot

# This script attmepts to provide the functions necessary to generate pairs of parameters 

BivariateGrid <- function(interval_1, interval_2, n_1, n_2){
  assertthat::assert_that(interval_2[1] >= 0.0, msg=paste('In', match.call()[1], 'interval2 (for diagonal element) must have non-negative elements'))
  assertthat::assert_that(length(interval_1) == 2)
  assertthat::assert_that(length(interval_2) == 2)
  assertthat::assert_that(diff(interval_2) >= 0)
  assertthat::assert_that(n_1 > 0)
  assertthat::assert_that(n_2 > 0)
  
  if(diff(interval_2) == 0){
      return(cbind(seq(interval_1[1], interval_1[2], length.out = n_1), rep(interval_2[1], n_1)))
  }else{
    #res <- matrix(0, ncol = 2, nrow = n_1 * n_2)
    # res <- , 
    #                     rep(interval_2[1], interval_2[2], length.out = n_1))
    # TODO HERE
    return(res)
  }
}


plot(BivariateGrid(c(0,1),c(2,5), 10, 5))

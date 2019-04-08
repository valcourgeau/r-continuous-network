setwd('GitHub/r-continuous-network/R/')
source('package_to_load.R')
source('adjacency_generation.R')

ConstructPath <- function(q_matrix, noise, y_init, delta_time){
  # generates path given noise increments and fixed delta_time.
  assertthat::assert_that(delta_time > 0)
  
  
}
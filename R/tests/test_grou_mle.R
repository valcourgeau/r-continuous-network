testthat::test_that('CoreNodeMLE_kron_numerator_unif_times', {
  N <- 5
  d <- 5
  times <- 1:N
  data <- matrix(rep(1:N, d), nrow = N, ncol = d, byrow = F)
  cn_mle <- CoreNodeMLE(times = times, data = data, thresholds = NA)
  perfect_denominator <- matrix(sum(1:(N-1)), d, d)
  testthat::expect_equal(cn_mle$numerator, perfect_denominator)
})

testthat::test_that('CoreNodeMLE_kron_numerator_doubling_unif_times', {
  N <- 5
  d <- 5
  times <- 1:N
  # component j have increments of size j
  data <- do.call(cbind, lapply(1:d, function(j){j*(1:N)}))
  cn_mle <- CoreNodeMLE(times = times, data = data, thresholds = NA)
  # i is for the increment, j for the component
  component_wise <- vapply(1:d, function(j){sum(j*(1:(N-1)))}, 1.0) #rows
  perfect_denominator <- do.call(cbind, lapply(1:d, function(i){i * component_wise})) #increments-wise, building cols
  testthat::expect_equal(cn_mle$numerator, perfect_denominator)
})
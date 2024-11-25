test_that("gibbs_bdp_dw runs on sample white noise data", {
  test_data <- rnorm(150)
  test_data <- test_data - mean(test_data)
  # just very few steps to make sure that the algorithm runs...
  mcmc <- gibbs_bdp_dw(data=test_data, 
                       m=5, 
                       likelihood_thinning=1, 
                       rescaled_time=(1:2)/3, 
                       freq=(1:2)/3 * pi,
                       Ntotal=2000, 
                       burnin=100, 
                       thin=1)
  # ... and produces results (i.e. tv-PSD estimate & log posterior trace)
  expect_equal(length(mcmc$lpost), 1900)
  expect_false(any(is.na(mcmc$lpost)))
  expect_true(length(mcmc$tvpsd.median)>0)
  expect_false(all(is.na(mcmc$tvpsd.median)))
})
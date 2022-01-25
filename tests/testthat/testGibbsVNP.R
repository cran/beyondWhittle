test_that("gibbs_vnp runs on sample white noise data", {
  test_data <- matrix(data=rnorm(200), nrow=100, ncol=2)
  test_data <- apply(test_data, 2, function(x) x-mean(x))
  # just very few steps to make sure that the algorithm runs...
  mcmc <- gibbs_vnp(data=test_data, Ntotal=2000, burnin=100, thin=1)
  # ... and produces results (i.e. PSD estimate & log posterior trace)
  expect_true(length(mcmc$lpost) > 0)
  expect_false(any(is.na(mcmc$lpost)))
  expect_true(length(mcmc$psd.median)>0)
  expect_false(all(is.na(mcmc$psd.median)))
})
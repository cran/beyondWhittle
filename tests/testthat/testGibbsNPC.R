test_that("gibbs_npc runs on sample white noise data", {
  test_data <- rnorm(100)
  test_data <- test_data - mean(test_data)
  # just very few steps to make sure that the algorithm runs...
  mcmc <- gibbs_npc(data=test_data, ar.order=2, Ntotal=2000, burnin=100, thin=1)
  # ... and produces results (i.e. PSD estimate & log posterior trace)
  expect_equal(length(mcmc$lpost), 1900)
  expect_false(any(is.na(mcmc$lpost)))
  expect_true(length(mcmc$psd.median)>0)
  expect_false(all(is.na(mcmc$psd.median)))
})
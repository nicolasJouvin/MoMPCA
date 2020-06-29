Q = 6
K = 4

N = 8

simu <- simulate_BBC(N = N, L = 100, epsilon = 0, theta_true = NULL)
dtm.full = simu$dtm.full
V = dim(dtm.full)[2]
Ytruth = simu$Ytruth

test_that("Yinit vector returns error", code = {
  Yinit = rep(1, N)
  expect_error(mmpca_clust_modelselect(dtm.full, Qs = c(2,3), Ks = 2, Yinit = Yinit))
})

test_that("beta.init matrix returns error", code = {
  V = ncol(dtm.full)
  beta.init = matrix(1/V, K, V)
  expect_error(mmpca_clust_modelselect(dtm.full, Qs = 2, Ks = c(2,3), init.beta = beta.init))
})

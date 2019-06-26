Q = 6
K = 4

theta_true = matrix(0, Q, K)
for (i in 1:K) {
  theta_true[i,i] = 1
}
theta_true[5,] = 1/2 * c(1, 0, 1, 0)
theta_true[6,] = 1/2 * c(0, 1, 0, 1)

N = 8

simu <- simulate_BBC(N = N, L = 100, epsilon = 0, theta_true = theta_true)
dtm.full = simu$dtm.full
V = dim(dtm.full)[2]
Ytruth = simu$Ytruth

res <- mmpca_clust(simu$dtm.full, Q = 6, K = 4,
                   Yinit = Ytruth,
                   max.epochs = 1,
                   keep = 1,
                   verbose = 0)

test_that("check slots", code = {
  expect_is(res, "mmpca_clust")
  expect_equal(res@K, 4)
  expect_equal(res@Q, 6)
  expect_length(unique(res@clustering), Q)
  expect_length(res@clustering, N)
  expect_true(res@icl < 0)
  expect_true(res@llhood < 0)
  # check increasing bounds
  expect_equal(res@logLikelihoods, sort(res@logLikelihoods, decreasing = F))
})

test_that("check dims", code = {
  expect_equal(dim(res@beta), c(K,V))
  expect_equal(dim(res@gamma), c(6, 4))
})

test_that("Check Beta and Gamma are proba matrix", code = {
  expect_equal(Matrix::rowSums(res@beta), rep(1, K))
  expect_equal(Matrix::rowSums(res@gamma), rep(1, Q))
})

test_that("Check wrong Yinit error message", code = {
  expect_error(mmpca_clust(simu$dtm.full, Q = 6, K = 4,
                           Yinit = c(1),
                           max.epochs = 1,
                           keep = 1,
                           verbose = 0))
  expect_error(mmpca_clust(simu$dtm.full, Q = 6, K = 4,
                           Yinit = rep(1, N),
                           max.epochs = 1,
                           keep = 1,
                           verbose = 0))
  expect_error(mmpca_clust(simu$dtm.full, Q = 6, K = 4,
              Yinit = 'wrong!',
              max.epochs = 1,
              keep = 1,
              verbose = 0))
})

test_that("Check wrong Yinit error message", code = {
  init.beta = matrix(1, nrow = 4, data = V)
  expect_error(mmpca_clust(simu$dtm.full, Q = 6, K = 4,
                           init.beta = init.beta,
                           max.epochs = 1,
                           keep = 1,
                           verbose = 0))
  expect_error(mmpca_clust(simu$dtm.full, Q = 6, K = 4,
                           init.beta = 'wrong !',
                           max.epochs = 1,
                           keep = 1,
                           verbose = 0))
})



skip_on_cran()






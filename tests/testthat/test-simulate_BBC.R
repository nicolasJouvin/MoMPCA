Q = 6
K = 4
N = 10
simu <- simulate_BBC(N = N, L = 100, epsilon = 0.1, lambda = 1, theta_true = NULL)


test_that("Check that simu dimensions are correct", code = {
  expect_equal(dim(simu$dtm.full)[1], N)
  expect_equal(length(simu$Ytruth), N)
  expect_equal(dim(simu$theta_true), c(Q,K))
})

test_that("All cluster are present", code = {
  expect_equal(length(unique(simu$Ytruth)), Q)
})


Q = 6
K = 4

theta_true = matrix(0, Q, K)
for (i in 1:K) {
  theta_true[i,i] = 1
}
theta_true[5,] = 1/2 * c(1, 0, 1, 0)
theta_true[6,] = 1/2 * c(0, 1, 0, 1)


test_that("Check that simu dimensions are correct", code = {

  N = 10
  simu <- simulate_BBC(N = N, L = 100, epsilon = 0.1, theta_true = theta_true)
  expect_equal(dim(simu$dtm.full)[1], N)
  expect_equal(dim(simu$Z), c(N, K))
  expect_equal(length(simu$Ytruth), N)
  expect_equal(simu$theta_true, theta_true)
})

test_that("All cluster are present", code = {
  N = 10
  simu <- simulate_BBC(N = N, L = 100, epsilon = 0.1, theta_true = theta_true)

  expect_equal(length(unique(simu$Ytruth)), Q)
})

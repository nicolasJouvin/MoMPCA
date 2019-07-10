
simu = simulate_BBC(N = 10, L = 10)
Q = 3
K = 2

test_that('Check init_Y w/ random', code = {
  Y = initialize_Y(simu$dtm.full, Q, K, init = 'random')
  expect_length(unique(Y), Q)
})

testthat::skip_on_cran()
test_that("Check init_Y w/ kmeans_lda", code = {
  Y = initialize_Y(simu$dtm.full, Q, K, init = 'kmeans_lda')
  expect_length(unique(Y), Q)
})

# test_that("Check init_Y w/ custom 'init' function HTSClust", code = {
#   Y = initialize_Y(simu$dtm.full, Q, K, init = benchmark.htsclust)
#   expect_length(unique(Y), Q)
# })
#
# test_that("Check init_Y w/ custom 'init' function multmix", code = {
#   Y = initialize_Y(simu$dtm.full, Q, K, init = benchmark.multmix)
#   expect_length(unique(Y), Q)
# })
#
# test_that("Check init_Y w/ custom 'init' function nmfem", code = {
#   Y = initialize_Y(simu$dtm.full, Q, K, init = benchmark.nmfem)
#   expect_length(unique(Y), Q)
# })

#
# test_that("Check init_Y w/ custom 'init' function NMF", code = {
#   Y = initialize_Y(simu$dtm.full, Q, K, init = benchmark.nmf)
#   expect_length(unique(Y), Q)
# })

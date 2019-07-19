
simu = simulate_BBC(N = 10, L = 10)
Q = 3
K = 2
dtm.full = simu$dtm.full
V = dim(dtm.full)[2]

test_that('Check init_Y w/ random', code = {
  Y = initialize_Y(dtm.full, Q, K, init = 'random')
  expect_length(unique(Y), Q)
})

test_that("Check init_Y w/ kmeans_lda", code = {
  testthat::skip_on_cran()
  Y = initialize_Y(dtm.full, Q, K, init = 'kmeans_lda')
  expect_length(unique(Y), Q)
})

test_that("Check initializeBeta w/ lda and given control", code = {
  testthat::skip_on_cran()
  control_test = methods::new("LDA_VEMcontrol",
                              var = methods::new('OPTcontrol',
                                                 iter.max = 1L, tol = 1e-1),
                              em = methods::new('OPTcontrol',
                                                iter.max = 1L, tol = 1e-1)
                              )
  beta = initializeBeta(dtm.full, 'lda', K,
                        verbose = 0,
                        control_lda_init = control_test)
  expect_equal( dim(beta), c(K, V))
  expect_equal(rowSums(beta), rep(1, K))
})

test_that("Check initializeBeta w/ lda and NULL control", code = {
  testthat::skip_on_cran()
  beta = initializeBeta(dtm.full, 'lda', K,
                        verbose = 0,
                        control_lda_init = NULL)
  expect_equal( dim(beta), c(K, V))
  expect_equal(rowSums(beta), rep(1, K))
})

# test_that("Check init_Y w/ custom 'init' function HTSClust", code = {
#   Y = initialize_Y(dtm.full, Q, K, init = benchmark.htsclust)
#   expect_length(unique(Y), Q)
# })
#
# test_that("Check init_Y w/ custom 'init' function multmix", code = {
#   Y = initialize_Y(dtm.full, Q, K, init = benchmark.multmix)
#   expect_length(unique(Y), Q)
# })
#
# test_that("Check init_Y w/ custom 'init' function nmfem", code = {
#   Y = initialize_Y(dtm.full, Q, K, init = benchmark.nmfem)
#   expect_length(unique(Y), Q)
# })

#
# test_that("Check init_Y w/ custom 'init' function NMF", code = {
#   Y = initialize_Y(dtm.full, Q, K, init = benchmark.nmf)
#   expect_length(unique(Y), Q)
# })

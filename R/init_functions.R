
#' @title Beta initialization
#' @description Used in the \code{\link{mmpca_clust}}() function to initialize
#'   beta. It can be either "random" or "lda". Please not that the
#'   \code{\link{mmpca_clust}}() function also allow for a user given beta
#'   matrix. In this case, this function is not used.
#' @param dtm An object of class \code{\link{DocumentTermMatrix}}
#' @param init.beta A string specifying the method, either \itemize{\item
#'   'random': Initialization a la Blei et. al. with 1/V coefficient everywhere
#'   + a small uniform noise U[0, 1e-10] on every coefficients. \item 'lda':
#'   Recommended. Uses the beta of LDA algorithm via a VEM algorithm, with an
#'   initialization of 5 repeats of the gibbs sampling algorithm with 1000
#'   burning iterations and 1000 iterations.  }
#' @param K The number of topics (dimension of the latent space).
#' @param verbose The verbosity level. Only prints a message at function
#'   activation.
#' @param control_lda_init The control for \code{\link{LDA}}().
#'   Only used when \code{init.beta} == 'lda' and initilialized to the default
#'   \code{"LDA_VEMcontrol"} of the \code{\linkS4class{TopicModelcontrol}}
#'   class.
#'
#' @return A KxV matrix with each row summing to 1.
#' @export
#' @examples
#' \donttest{
#' simu = simulate_BBC(N = 100, L = 100)
#' K = 4
#' beta = initializeBeta(simu$dtm.full, 'lda', K, verbose = 1)
#' }
initializeBeta = function(dtm, init.beta, K,
                          verbose = 0,
                          control_lda_init = NULL) {

  V = dim(dtm)[2]

  if (verbose > 0) message('Init beta with the ', init.beta, ' method...')
  if (init.beta == 'random') {

    coefs = rep(1/V, K*V) + stats::runif(n = K*V, min = 0, max = 1e-10)
    coefs = matrix(coefs, nrow = K, ncol = V)
    beta.init = coefs / rowSums(coefs)

  } else if (init.beta == 'lda') {
    ## Combination of Gibbs sampling + VEM to find the best beta
    if (is.null(control_lda_init))
      control_lda_init = methods::new("LDA_VEMcontrol")

    # 1. Gibbs
    nstart.gibbs = 5
    seed.gibbs  = rep(NA, nstart.gibbs) # See ?TopicModelcontrol-class --> seed field explanation
    control_gibbs = list(alpha = control_lda_init@alpha,
                              estimate.beta = TRUE,
                              nstart = nstart.gibbs,
                              seed = seed.gibbs,
                              burnin = 1000,
                              iter = 1000,  # nbre iter after burnin
                              thin = 10,
                              best = TRUE, # keep sample according to best posterior
                              verbose = 0
                              )
    lda.init.gibbs = topicmodels::LDA(dtm, k = K, method = 'gibbs',
                                      control = control_gibbs)
    # 2. VEM
    lda.init = topicmodels::LDA(dtm, k = K, method = 'VEM',
                                control = control_lda_init,
                                model = lda.init.gibbs)

    ## exponentiate and renormalize beta (for underflows)
    beta.init = exp(lda.init@beta) / rowSums(exp(lda.init@beta))

  } else {
    stop(paste0(init.beta, ' method is not currently implemented (check for typos).'))
  }

  if (verbose > 0) message('... Finished. ')
  beta.init
}

##****************************
## Initialize clustering

#' @title Clustering initialization
#' @description Perform a \code{DocumentTermMatrix} clustering via default
#'   routines or allow for user specified function
#'
#' @param dtm An object of class \code{\link{DocumentTermMatrix}}
#' @param Q The number of cluster
#' @param K The dimension of the latent space. It is mandatory, for
#'   compatibility reasons but not always used (e.g. random do not use it).
#' @param init Either: \itemize{\item \code{'random'}: Random initialization.
#'   \item \code{'kmeans_lda'}: A Q-kmeans on the latent space (theta matrix) of
#'   a K-topic LDA.  \item A user defined function which MUST take the following
#'   structure for compatibility \code{init <- function(dtm, Q, K, nruns, ...)}}
#'
#' @return A vector of size equal to the number of row of \code{dtm}, containing
#'   a Q-clustering
#' @details For more details see \code{\link{benchmarks-functions}}
#' @export
#'
#' @examples
#' \donttest{
#' simu = simulate_BBC(N = 100, L = 100)
#' Q = 6
#' K = 4
#' Y = initialize_Y(simu$dtm.full, Q, K, init = 'kmeans_lda')
#' }
initialize_Y <- function(dtm, Q, K, init='random'){

  INIT_registry <- list(benchmark.kmeans_lda = c("kmeans_lda",
                                                 "kmeans+lda",
                                                 "kmeans.lda"),
                        benchmark.random = c('random',
                                             'rand'))

  if (!is.function(init)) {
    MATCH <- which(sapply(INIT_registry, function(x) length(grep(tolower(init), tolower(x)))) > 0)
    if (!length(MATCH) == 1)
      stop("'init'", init, "not specified correctly")
    method <- get(names(INIT_registry)[MATCH])
  } else {
    method <- init
  }

  method(dtm = dtm, Q = Q, K = K, nruns = 1)
}

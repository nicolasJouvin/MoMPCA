

##****************************
## Initialize matrix beta

initializeBeta = function(dtm, init.beta, K, verbose, control_lda_init) {
  #' Function to calculate Beta on a given corpus
  #' Several choices are available
  #' @return a KxV matrix that which rows sums to 1

  V = dim(dtm)[2]

  if (init.beta == 'random') {

    if (verbose > 0) cat('\n Init beta randomly...')
    coefs = rep(1/V, K*V) + stats::runif(n = K*V, min = 0, max = 1e-10)
    coefs = matrix(coefs, nrow = K, ncol = V)
    beta.init = coefs / rowSums(coefs)

  } else if (init.beta == 'lda') {
    ## Combination of Gibbs sampling + VEM to find the best beta
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

    ## renormalize beta for underflows
    beta.init = exp(lda.init@beta) / rowSums(exp(lda.init@beta))

  } else if (init.beta == 'nmf') {
    if (!requireNamespace("NMF", quietly = TRUE)) {
      stop('Package NMF needed for this initialization function to work. Please install it.',
           call. = FALSE)
    }
    ## Use NMF alogorithm of Lee and Seung. In parrallel if possible.
    X.normalized = as.matrix(dtm) / apply(dtm, 1, function(x) sqrt(sum(x^2)))
    res.nmf = NMF::nmf(x = t(X.normalized), method = 'lee', rank = K, nrun = 5, .opt = 'p')
    beta.init = t(NMF::basis(res.nmf))

  } else {
    stop(paste0(init.beta, ' method is not currently implemented (check for typos).'))
  }

  return(beta.init)
}

##****************************
## Initialize clustering

initialize_Y <- function(dtm, Q, K, init='random'){
  if (init == 'random') {
    N = dim(dtm)[1]
    Y = apply(rmultinom(N, 1, rep(1/Q, Q)), 2, which.max)
    # check if #meta-docs < Q
    while (length(unique(Y)) != Q) Y = apply(rmultinom(N, 1, rep(1/Q, Q)), 2, which.max)
  } else if (init == 'HTSCluster') {
    conds = 1:V
    norm = rep(1,V)
    run =  HTSCluster::PoisMixClus(as.matrix(dtm), g = Q, conds = conds, norm = norm, init.type = 'small-em')
    Y = run$labels
  } else if (init == 'nmf') {
    Y = benchmark.nmf(dtm.full = dtm, Q = Q)
  } else if (init == 'nmf_kmeans') {
    Y = benchmark.nmf_kmeans(dtm.full = dtm, Q = Q, K = K)
  } else if (init == 'kmeans_lda') {
    Y = benchmark.kmeans_lda(dtm.full = dtm, Q = Q, K = K)
  } else if (init == 'gmm_lda') {
    Y = benchmark.gmm_lda(dtm = dtm, Q = Q, K = K)
  } else {
    stop(paste0('Method ', init, ' is not implemented yet. Try supplying the clustering vector
                to Yinit directly.'))
  }

  return(Y)
}


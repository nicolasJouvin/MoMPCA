
initializeBeta = function(dtm, init.beta, K, verbose, control_gibbs_init, control_lda_init) {
  #' Function to calculate Beta on a given corpus
  #' Several choices are available
  #' @return a KxV matrix that which rows sums to 1

  V = dim(dtm)[2]

  if (init.beta == 'random') {

    if (verbose > 0) cat('\n Init beta randomly...')
    coefs = rep(1/V, K*V) + stats::runif(n = K*V, min = 0, max = 1e-10)
    coefs = matrix(coefs, nrow = K, ncol = V)
    coefs = coefs / rowSums(coefs)
    lbeta.init = log(matrix(coefs, nrow = K, ncol = V))
    beta.init = exp(lbeta.init)

  } else if (init.beta == 'lda') {

    lda.init.gibbs = topicmodels::LDA(dtm, k = K, method = 'gibbs', control = control_gibbs_init)
    lda.init = topicmodels::LDA(dtm, k = K, method = 'VEM', control = control_lda_init,
                                model = lda.init.gibbs)
    beta.init = exp(lda.init@beta) / rowSums(exp(lda.init@beta))

  } else if (init.beta == 'nmf') {

    X.normalized = as.matrix(dtm) / apply(dtm, 1, function(x) sqrt(sum(x^2)))

    res.nmf = NMF::nmf(x = t(X.normalized), method = 'lee', rank = K, nrun = 10, .opt = 'p')
    beta.init = t(basis(res.nmf))

  } else {
    stop(paste0(init.beta, ' method is not currently implemented (check for typos).'))
  }

  return(beta.init)
}


initialize_Y <- function(dtm, Q, K, init='random'){

  if(init=='random'){

    Y = apply(rmultinom(D, 1, rep(1/Q, Q)), 2, which.max) # random init
    while(length(unique(Y)) != Q) Y = apply(rmultinom(D, 1, rep(1/Q, Q)), 2, which.max) # check if #meta-docs < Q

  }else if(init == 'HTSCluster'){

    conds = 1:V
    norm = rep(1,V)
    run =  PoisMixClus(as.matrix(dtm), g = Q, conds=conds, norm=norm, init.type = 'small-em')
    Y = run$labels

  }else if(init == 'nmf'){

    Y = benchmark.nmf(dtm.full = dtm, Q = Q)

  } else if(init == 'nmf_kmeans'){

    Y = benchmark.nmf_kmeans(dtm.full = dtm, Q = Q, K=K)

  }
  else if (init == 'kmeans_lda'){

    Y = benchmark.kmeans_lda(dtm.full = dtm, Q = Q, K = K)

  }else if(init == 'truth'){

    Y = Ytruth
    cat(' / !! \ On init avec Ytruth. / !! \ \n ')

  }

  return(Y)
}


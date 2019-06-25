## TODO add NAMESPACE checks to see if packages are installed


######################
# Benchmarks
benchmark.htsclust = function(dtm.full, Q, ...){
  #' Clustering of count data via the Poisson mixture model of Rau et. al.
  #' @param ... : params to pass onto the PoisMixClus function of the HTSCluster package

  if (!requireNamespace("HTSCluster", quietly = T)) {
    stop('Package HTSCluster needed for this initialization function to work. Please install it.',
         call. = FALSE)
  }

  V = dim(dtm.full)[2]
  conds = 1:V
  norm = rep(1,V) # optional, avoid numerical error due to normalization by 0 when a word doesn't appear in the corpus (rare).
  run <- HTSCluster::PoisMixClus(as.matrix(dtm.full), g = Q , conds = conds, norm = norm, alg.type = 'CEM', ...)

  return(run$labels)


}

benchmark.nmf = function(dtm.full, Q, ...){

  if (!requireNamespace("NMF", quietly = TRUE)) {
    stop('Package NMF needed for this initialization function to work. Please install it.',
         call. = FALSE)
  }

  Wtilde = as.matrix(dtm.full)
  dtm.tfidf.normalized =  suppressWarnings(tm::as.DocumentTermMatrix(Wtilde, weighting = tm::weightTfIdf))
  stopifnot(sum(is.na(dtm.tfidf.normalized$v)) == 0)
  #Sanity check : remove words that never appear in the corpus
  voc.zero =  which(slam::col_sums(dtm.tfidf.normalized) == 0)
  if (length(voc.zero) != 0) {
    dtm.tfidf.normalized = dtm.tfidf.normalized[,-voc.zero]
  }
  remove(Wtilde)
  X.normalized = as.matrix(dtm.tfidf.normalized) / apply(dtm.tfidf.normalized, 1, function(x) sqrt(sum(x^2)))

  res.nmf = NMF::nmf(x = t(X.normalized), method = 'lee', rank = Q, ...)

  H = coef(res.nmf)

  Y.est.nmf = apply(H, 2 , which.max)

  return(Y.est.nmf)

}


benchmark.LDA = function(dtm.full, Q, ...){
  # Cluster with MAP on theta, theta is estimated with a LDA topicmodels with Q topics
  baseline.lda = topicmodels::LDA(dtm.full,
                                  control = list(estimate.alpha = FALSE,
                                                 estimate.beta = TRUE,
                                                 alpha = 1,
                                                 verbose = 0,
                                                 nstart = nruns,
                                                 var = list(iter.max = 5000)
                                  ),
                                  k = Q)
  Y.baseline.lda = apply(baseline.lda@gamma, 1, which.max)

  return(Y.baseline.lda)
}

benchmark.kmeans_lda = function(dtm.full, Q, K, ...){
  #' Cluster the matrix theta obtained by a topicmodels LDA with K.est topics

  baseline.lda = topicmodels::LDA(dtm.full,
                                  k = K,
                                  control = list(estimate.alpha = FALSE,
                                                 estimate.beta = TRUE,
                                                 alpha = 1,
                                                 verbose = 0,
                                                 nstart = 4,
                                                 var = list(iter.max = 5000)
                                  )
  )
  Y.baseline.kmeans_theta = stats::kmeans(baseline.lda@gamma, centers = Q, ...)$cluster

  return(Y.baseline.kmeans_theta)
}


benchmark.kmeans_dtm = function(dtm.full, Q, ...){
  #' Naive Clustering on the dtm directly

  Y.baseline.kmeans_theta = stats::kmeans(dtm.full, centers = Q, ...)$cluster
  return(Y.baseline.kmeans_theta)
}

benchmark.multmix = function(dtm.full, Q, nruns) {
  #' Clustering with the mixture of multinomials models

  if (!requireNamespace("mixtools", quietly = T)) {
    stop('Package mixmult needed for this initialization function to work. Please install it.',
         call. = FALSE)
  }

  X = as.matrix(DTMtoSparse(dtm.full))

  results = list()

  if (requireNamespace("foreach", quietly = TRUE)) {
    results <- foreach(i = 1:nruns, .inorder = F) %dopar% {
      invisible(capture.output( res <- mixtools::multmixEM( X, k = Q, verb = F)))
      return(res)
    }
  } else {
    for (i in 1:nruns)
      invisible(capture.output( res[[i]] <- mixtools::multmixEM( X, k = Q, verb = F)))
  }

  best_res = unlist(sapply(1:nruns, function(i) results[[i]]$loglik))
  best_res = which.max(best_res)
  Y.est.mixtmult = apply(results[[best_res]]$posterior, 1, which.max)

  return(Y.est.mixtmult)

}

benchmark.nmfem = function(dtm, Q, K, nruns) {
  #' NMFEM algo from carel et. al.
  #' MLE in a frequentist setting of MMPCA model

  if (!requireNamespace("nmfem", quietly = TRUE)) {
    stop('Package nmfem needed for this initialization function to work. Please install it.',
         call. = FALSE)
  }

  X = as.data.frame(as.matrix(dtm))

  results = list()

  if (requireNamespace("foreach", quietly = TRUE)) {
    results <- foreach(i = 1:nruns, .inorder = F) %dopar% {
      invisible(capture.output(res <- tryCatch( expr = {
        # try runnning nmfem
        nmfem::nmfem_mult(X = X, H = K, K = Q)
      },
      error = function(e) {
        # Return list with $llh = -in if it fails
        list(llh = -Inf)
      })
      ))
      return(res)
    }
  } else {
    for (i in 1:nruns)
      invisible(capture.output( res[[i]] <- nmfem::nmfem_mult(X = X, H = K, K = Q)))
  }

  best_res = unlist(sapply(1:nruns, function(i) results[[i]]$llh))
  best_res = which.max(best_res)
  Y.nmfem =  apply(results[[best_res]]$posterior, 1, which.max)
  return(Y.nmfem)


}

benchmark.gmm_lda = function(dtm, Q, K , nruns=1, seed=NULL) {
  #' Fit a Q-GMM in a K-LDA's latent space (Rmixmod)

  if (!requireNamespace("Rmixmod", quietly = TRUE)) {
    stop('Package Rmixmod needed for this initialization function to work. Please install it.',
         call. = FALSE)
  }

  baseline.lda = topicmodels::LDA(dtm,
                                  control = list(estimate.alpha = FALSE,
                                                 estimate.beta = TRUE,
                                                 alpha = 1,
                                                 verbose = 0,
                                                 nstart = 4,
                                                 var = list(iter.max = 5000)),
                                  k = K)
  strat = Rmixmod::mixmodStrategy(algo = 'CEM', nbTry = nruns)
  tryModels = Rmixmod::mixmodGaussianModel(family = 'all')
  res = Rmixmod::mixmodCluster(data = as.data.frame(baseline.lda@gamma),
                               nbCluster = Q,
                               models = tryModels,
                               strategy = strat,
                               criterion = 'ICL',
                               seed = seed)

  Y.est.gmm_mpca = res@bestResult@partition
  return(Y.est.gmm_mpca)
}

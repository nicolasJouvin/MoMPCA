
#' convert a dtm from package tm to sparseMatrix from package Matrix
#' without converting it to full matrix.
#'@param dtm a document-term-matrix from package 'tm'
#'
#'@return a sparse dcgMatrix from package Matrix
DTMtoSparse <- function(dtm){

  return(Matrix::sparseMatrix(i = dtm$i, j = dtm$j, x = dtm$v, dimnames = dtm$dimnames, dims = dim(dtm)))
}


#' Wrapper on mmpca_clust for several restart
#'
#' Returns the best of all runs, i.e. the one wich achieves the greatest bound.
#' It works with the foreach package, if a parralel backend is registered it will run
#' the different initializations in parrallel.
#' @param nruns number of restart of the algo
#' @param parrallel True if a parrallel or False
#' @param Q.est The number of cluster
#' @param K.est The number of topics
#'
#' @return
#' @export
#'
#' @examples
wrapper.mmpca_clust = function(nruns = 5, Q.est, K.est){
  #' Wrapper function on mLDA to evaluate it on BBC (or other data set) with several restarts
  models = list()


  models <- foreach::foreach(i=1:nruns,
                    .inorder = F,
                    .export = ls(globalenv()),
                    .packages = c('topicmodels', 'Matrix', 'slam', 'aricode', 'NMF')
  ) %dopar% {
    model = mmpca_clust(dtm.full, Q = Q.est, K = K.est, Ytruth=Ytruth,
                             control_lda_loop = control_lda_loop,
                             max.iter = max.iter,
                             control_lda_init=control_lda_init,
                             init=init,
                             method='greedy',
                             verbose=verbose)
    return(model)
  }

  bounds = sapply(1:nruns, function(n) models[[n]]$final_bound)
  best_model = model[[which.max(bounds)]]

  return(best_model)
}


######################
# Benchmarks
benchmark.htsclust = function(dtm.full, Q, ...){
  #' Clustering of count data via the Poisson mixture model of Rau et. al.
  #' @param ... : params to pass onto the PoisMixClus function of the HTSCluster package

  V = dim(dtm.full)[2]
  conds = 1:V
  norm = rep(1,V) # optional, avoid numerical error due to normalization by 0 when a word doesn't appear in the corpus (rare).
  run <- HTSCluster::PoisMixClus(as.matrix(dtm.full), g = Q , conds = conds, norm = norm, alg.type = 'CEM', ...)

  return(run$labels)
}

benchmark.nmf = function(dtm.full, Q, ...){
  Wtilde = as.matrix(dtm.full)
  dtm.tfidf.normalized =  suppressWarnings(tm::as.DocumentTermMatrix(Wtilde, weighting = weightTfIdf))
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

  X = as.matrix(DTMtoSparse(dtm.full))

  results = list()
  results <- foreach(i = 1:nruns, .inorder = F) %dopar% {
    invisible(capture.output( res <- mixtools::multmixEM( X, k = Q, verb = F)))
    return(res)
  }

  best_res = unlist(sapply(1:nruns, function(i) results[[i]]$loglik))
  best_res = which.max(best_res)
  Y.est.mixtmult = apply(results[[best_res]]$posterior, 1, which.max)

  return(Y.est.mixtmult)
}

benchmark.nmfem = function(dtm, Q, K, nruns) {
  #' NMFEM algo from carel et. al.
  #' MLE in a frequentist setting of MMPCA model

  X = as.data.frame(as.matrix(dtm))

  results = list()
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

  best_res = unlist(sapply(1:nruns, function(i) results[[i]]$llh))
  best_res = which.max(best_res)
  Y.nmfem =  apply(results[[best_res]]$posterior, 1, which.max)
  return(Y.nmfem)
}

benchmark.mpca_gmm = function(dtm, Q, K , nruns=1, seed=NULL) {
  #' Fit a Q-GMM in a K-LDA's latent space (Rmixmod)

  baseline.lda = topicmodels::LDA(dtm,
                     control = list(estimate.alpha=FALSE, estimate.beta=TRUE, alpha = 1, verbose=0,
                                    nstart=4, var=list(iter.max=5000)),
                     k = K)
  strat = Rmixmod::mixmodStrategy(algo = 'CEM', nbTry = nruns)
  tryModels = Rmixmod::mixmodGaussianModel(family = 'all')
  res = Rmixmod::mixmodCluster(data=as.data.frame(baseline.lda@gamma),
                      nbCluster = Q,
                      models = tryModels, strategy = strat, criterion = 'ICL', seed = seed)

  Y.est.gmm_mpca = res@bestResult@partition

  return(Y.est.gmm_mpca)
}


##############################################
# Text analysis

preprocess_corpus = function(corpus, stop_words='english', language='english'){
  corpus <- tm_map(corpus,content_transformer(tolower))
  corpus <- tm_map(corpus, content_transformer(removeNumbers))
  corpus <- tm_map(corpus, removePunctuation)
  corpus <- tm_map(corpus, stemDocument, language = language)
  corpus <- tm_map(corpus, stripWhitespace)
  return(corpus)
}



#############################################
# Init functions for mmpca_clust()
#############################################

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


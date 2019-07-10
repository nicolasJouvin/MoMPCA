#' @title Benchmarks functions for clustering
#' @description These are wrapper to other methods for the clustering of count
#'   data. They can be used to initialize the clustering. It is also possible to
#'   implement your own benchmark function depending on other packages.
#' @param dtm an S4 object of class \code{\linkS4class{mmpcaClust}}
#' @param Q The number of clusters
#' @param K Number of topics (dimension of the latent space).
#' @param nruns Number of restart of the \code{\link[stats]{kmeans}}()
#'   algorithm.
#' @param ... Some argument to be consistent with the function's skeleton :
#'   \code{K} and \code{nruns} are optional arguments for some of them.
#' @return A vector of size equal to the number of row of \code{dtm}, containing
#'   a Q-clustering
#' @rdname benchmarks-functions
#' @name benchmarks-functions
NULL

#' @section \code{benchmark.random}:
#' Random initialisation of the clustering. Arguments K and nruns are unused
#' @rdname benchmarks-functions
#' @export
benchmark.random <- function(dtm, Q, ...) {
  N = dim(dtm)[1]
  Y = apply(stats::rmultinom(N, 1, rep(1/Q, Q)), 2, which.max)
  # check if #meta-docs < Q
  while (length(unique(Y)) != Q)
    Y = apply(stats::rmultinom(N, 1, rep(1/Q, Q)), 2, which.max)
  Y
}



#' @section \code{benchmarks.kmeans_lda}: Cluster the matrix theta obtained by a
#'   topicmodels LDA with K topics
#' @rdname benchmarks-functions
#' @export
benchmark.kmeans_lda = function(dtm, Q, K, nruns = 1, ...){

  baseline.lda = topicmodels::LDA(dtm,
                                  k = K,
                                  control = list(estimate.alpha = FALSE,
                                                 estimate.beta = TRUE,
                                                 alpha = 1,
                                                 verbose = 0,
                                                 nstart = 4,
                                                 var = list(iter.max = 5000)
                                  )
  )
  Y.baseline.kmeans_theta = stats::kmeans(baseline.lda@gamma,
                                          centers = Q, nstart = nruns)$cluster

  return(Y.baseline.kmeans_theta)
}



#' #' @section benchmark.htsclust:
#' #' Clustering of count data via the Poisson mixture model of Rau et. al.
#' #' @rdname benchmarks-functions
#' #' @param dtm
#' #' @param Q
#' #' @param ... : params to pass onto the PoisMixClus function of the HTSCluster package
#' #' @return
#' #' @export
#' #' @
#' #' @examples
#' benchmark.htsclust = function(dtm, Q, nruns=1, ...){
#'
#'   if (!requireNamespace("HTSCluster", quietly = T)) {
#'     stop('Package HTSCluster needed for this initialization function to work. Please install it.',
#'          call. = FALSE)
#'   }
#'
#'   loadNamespace('HTSCluster')
#'   V = dim(dtm)[2]
#'   conds = 1:V
#'   norm = rep(1,V) # optional, avoid numerical error due to normalization by 0 when a word doesn't appear in the corpus (rare).
#'   run <- HTSCluster::PoisMixClus(as.matrix(dtm),
#'                                  g = Q ,
#'                                  conds = conds,
#'                                  norm = norm,
#'                                  alg.type = 'CEM',
#'                                  init.runs = nruns)
#'   unloadNamespace('HTSCluster')
#'   return(run$labels)
#'
#'
#' }
#'
#'
#'
#' #' @section \code{benchmarks.nmf}:
#' #' Perform NMF document clustering of Lee and Seung .
#' #' @rdname benchmarks-functions
#' #' @param dtm
#' #' @param Q
#' #' @param ...
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' benchmark.nmf = function(dtm, Q, nruns = 1, ...){
#'
#'   if (!requireNamespace("NMF", quietly = TRUE)) {
#'     stop('Package NMF needed for this initialization function to work. Please install it.',
#'          call. = FALSE)
#'   }
#'
#'   Wtilde = as.matrix(dtm)
#'   dtm.tfidf.normalized =  suppressWarnings(tm::as.DocumentTermMatrix(Wtilde, weighting = tm::weightTfIdf))
#'   stopifnot(sum(is.na(dtm.tfidf.normalized$v)) == 0)
#'   #Sanity check : remove words that never appear in the corpus
#'   voc.zero =  which(slam::col_sums(dtm.tfidf.normalized) == 0)
#'   if (length(voc.zero) != 0) {
#'     dtm.tfidf.normalized = dtm.tfidf.normalized[,-voc.zero]
#'   }
#'   remove(Wtilde)
#'   X.normalized = as.matrix(dtm.tfidf.normalized) / apply(dtm.tfidf.normalized, 1, function(x) sqrt(sum(x^2)))
#'
#'   res.nmf = NMF::nmf(x = t(X.normalized), method = 'lee', rank = Q,
#'                      nrun = nruns, .options = 'p')
#'
#'   H = NMF::coef(res.nmf)
#'
#'   Y.est.nmf = apply(H, 2 , which.max)
#'
#'   return(Y.est.nmf)
#'
#' }
#'
#'
#' #' @title benchmark.LDA Cluster with MAP on theta, theta is estimated with a
#' #'   LDA topicmodels with Q topics
#' #' @rdname benchmarks-functions
#' #'
#' #' @param dtm
#' #' @param Q
#' #' @param nruns
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' benchmark.LDA = function(dtm, Q, nruns = 1, ...){
#'   baseline.lda = topicmodels::LDA(dtm,
#'                                   control = list(estimate.alpha = FALSE,
#'                                                  estimate.beta = TRUE,
#'                                                  alpha = 1,
#'                                                  verbose = 0,
#'                                                  nstart = nruns,
#'                                                  var = list(iter.max = 5000)
#'                                   ),
#'                                   k = Q)
#'   Y.baseline.lda = apply(baseline.lda@gamma, 1, which.max)
#'
#'   return(Y.baseline.lda)
#' }
#'
#'
#' #' @title benchmark.kmeans_dtm
#' #' @section \code{benchmarks.kmeans_dtm}:
#' #' Naive Clustering on the dtm directly
#' #' @rdname benchmarks-functions
#' #' @param dtm
#' #' @param Q
#' #' @param ...
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' benchmark.kmeans_dtm = function(dtm, Q, ...){
#'
#'   Y.baseline.kmeans_theta = stats::kmeans(dtm, centers = Q, ...)$cluster
#'   return(Y.baseline.kmeans_theta)
#' }
#'
#'
#'
#' #' @section \code{benchmark.multmix} Cluster dtm through the mixture of multinomials models.
#' #'   Essentially a wrapper around \code{\link{mixtools}{multmix}}().
#' #' @rdname benchmarks-functions
#' #' @param dtm
#' #'
#' #' @param Q
#' #' @param nruns
#' #'
#' #' @importFrom utils capture.output
#' #' @importFrom foreach %dopar%
#' benchmark.multmix = function(dtm, Q, nruns = 1, ...) {
#'
#'   if (!requireNamespace("mixtools", quietly = T)) {
#'     stop('Package mixmult needed for this initialization function to work. Please install it.',
#'          call. = FALSE)
#'   }
#'
#'   loadNamespace("mixtools")
#'   X = as.matrix(DTMtoSparse(dtm))
#'
#'   results = list()
#'
#'   if (requireNamespace("foreach", quietly = TRUE)) {
#'     loadNamespace("foreach")
#'     results <- foreach::foreach(i = 1:nruns, .inorder = F) %dopar% {
#'       invisible(capture.output( res <- mixtools::multmixEM( X, k = Q, verb = F)))
#'       return(res)
#'     }
#'     unloadNamespace("foreach")
#'   } else {
#'     for (i in 1:nruns)
#'       invisible(capture.output( res[[i]] <- mixtools::multmixEM( X, k = Q, verb = F)))
#'   }
#'
#'   best_res = unlist(sapply(1:nruns, function(i) results[[i]]$loglik))
#'   best_res = which.max(best_res)
#'   Y.est.mixtmult = apply(results[[best_res]]$posterior, 1, which.max)
#'   unloadNamespace("mixtools")
#'
#'   Y.est.mixtmult
#' }
#'
#'
#' #' @section \code{benchmark.nmfem} NMFEM algo from carel et. al. Amounts to do MLE in a frequentist
#' #'   setting of MMPCA model.
#' #' @rdname benchmarks-functions
#' #'
#' #' @param dtm
#' #' @param Q
#' #' @param K
#' #' @param nruns
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' #' @importFrom foreach %dopar%
#' benchmark.nmfem = function(dtm, Q, K, nruns, ...) {
#'
#'
#'   if (!requireNamespace("nmfem", quietly = TRUE)) {
#'     stop('Package nmfem needed for this initialization function to work. Please install it.',
#'          call. = FALSE)
#'   }
#'
#'   loadNamespace("nmfem")
#'   X = as.data.frame(as.matrix(dtm))
#'
#'   results = list()
#'
#'   if (requireNamespace("foreach", quietly = TRUE)) {
#'     loadNamespace("foreach")
#'     results <- foreach::foreach(i = 1:nruns, .inorder = F) %dopar% {
#'       invisible(capture.output(res <- tryCatch( expr = {
#'         # try runnning nmfem
#'         nmfem::nmfem_mult(X = X, H = K, K = Q)
#'       },
#'       error = function(e) {
#'         # Return list with $llh = -in if it fails
#'         list(llh = -Inf)
#'       })
#'       ))
#'       return(res)
#'     }
#'     unloadNamespace("foreach")
#'   } else {
#'     for (i in 1:nruns)
#'       invisible(capture.output( res[[i]] <- nmfem::nmfem_mult(X = X, H = K, K = Q)))
#'   }
#'
#'   best_res = unlist(sapply(1:nruns, function(i) results[[i]]$llh))
#'   best_res = which.max(best_res)
#'   Y.nmfem =  apply(results[[best_res]]$posterior, 1, which.max)
#'
#'   unloadNamespace("nmfem")
#'   Y.nmfem
#' }
#'
#'
#' #' @title \code{benchmark.gmm_lda}() Fit a Q-GMM in a K-LDA's latent space (Rmixmod)
#' #' @rdname benchmarks-functions
#' #' @section benchmarks.gmm_lda:
#' #' dsncs
#' #' sdcknsqc
#' #'
#' #' sdlkvnsd
#' #'
#' #' @param dtm
#' #' @param Q
#' #' @param K
#' #' @param nruns
#' #' @param seed
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' benchmark.gmm_lda = function(dtm, Q, K , nruns=1, seed=NULL, ...) {
#'
#'   if (!requireNamespace("Rmixmod", quietly = TRUE)) {
#'     stop('Package Rmixmod needed for this initialization function to work. Please install it.',
#'          call. = FALSE)
#'   }
#'
#'   loadNamespace("Rmixmod")
#'   baseline.lda = topicmodels::LDA(dtm,
#'                                   control = list(estimate.alpha = FALSE,
#'                                                  estimate.beta = TRUE,
#'                                                  alpha = 1,
#'                                                  verbose = 0,
#'                                                  nstart = 4,
#'                                                  var = list(iter.max = 5000)),
#'                                   k = K)
#'   strat = Rmixmod::mixmodStrategy(algo = 'CEM', nbTry = nruns)
#'   tryModels = Rmixmod::mixmodGaussianModel(family = 'all')
#'   res = Rmixmod::mixmodCluster(data = as.data.frame(baseline.lda@gamma),
#'                                nbCluster = Q,
#'                                models = tryModels,
#'                                strategy = strat,
#'                                criterion = 'ICL',
#'                                seed = seed)
#'
#'   Y.est.gmm_mpca = res@bestResult@partition
#'   unloadNamespace("Rmixmod")
#'
#'   Y.est.gmm_mpca
#' }

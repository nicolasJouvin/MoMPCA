#' @title Greedy procedures for joint inference and clustering in MMPCA
#' @description Perform clustering of count data using the MMPCA model.
#'
#' @param dtm an NxV \code{\link{DocumentTermMatrix}} with term-frequency
#'   weighting.
#' @param Q The number of clusters
#' @param K The number of topics (latent space dimension)
#' @param model A given model in which to take the controls for the VE-steps in
#'   the greedy procedure. If NULL, a model of class
#'   \code{\linkS4class{mmpcaClust}} is created with default controls (see
#'   \code{\linkS4class{mmpcaClustcontrol}} class for more details).
#' @param Yinit Parameter for the initialization of Y. It can be either:
#'   \itemize{ \item a string or a function specifying the initialization
#'   procedure. It should be one of ('random', 'kmeans_lda'). See
#'   \code{\link{benchmarks-functions}} functions for more details. \item A
#'   vector of length N with Q modalities, specifying the initialization
#'   clustering. }
#' @param method The clustering algorithm to be used. Only "BBCVEM" is available
#'   : it corresponds to the branch and bound C-VEM of the original article.
#' @param init.beta Parameter for the initialization of the matrix beta. It can
#'   be either: \itemize{ \item a string specifying the initialization
#'   procedure. It should be one of ('random', 'lda'). See
#'   \code{\link{initializeBeta}}() for more details. \item A KxV matrix with
#'   each row summing to 1.}
#' @param keep The evolution of the bound is tracked every \code{keep} iteration
#' @param max.epochs Specifies the maximum number of pass allowed on the whole
#'   dataset.
#' @param verbose verbosity level
#' @param nruns number of runs of the algorithm (default to 1) : the run
#'   achieving the best evidence lower bound is selected.
#' @param mc.cores The number of CPUs to use when fitting in parallel the different
#'   models (only for non-Windows platforms). Default is the number of available
#'   cores minus 1.
#'
#' @return An object of class \code{"\linkS4class{mmpcaClust}"} containing the
#'   fitted model.
#' @export
#' @importFrom parallel detectCores
#' @importFrom foreach %dopar%
mmpca_clust = function(dtm,
                       Q,
                       K,
                       model = NULL,
                       Yinit = 'random',
                       method = 'BBCVEM',
                       init.beta = 'lda',
                       keep = 1L,
                       max.epochs = 10L,
                       verbose = 1L,
                       nruns = 1L,
                       mc.cores = max(1, (detectCores() - 1))) {


  cl <- parallel::makeCluster(mc.cores)
  doParallel::registerDoParallel(cl)


  RES <- foreach::foreach(i = 1:nruns,
                        .inorder = F,
                        .export = c("bbcvem")
                        # .packages = c('topicmodels', 'Matrix', 'slam', 'aricode', 'NMF')
                        ) %dopar% {
                          res = bbcvem(dtm,
                                       Q,
                                       K,
                                       model = model,
                                       Yinit = Yinit,
                                       method = method,
                                       init.beta = init.beta,
                                       keep = keep,
                                       max.epochs = max.epochs,
                                       verbose = verbose)
                          return(res)
                }

  parallel::stopCluster(cl)
  rm(cl)
  # return the run with the best evidence lower bound
  ELBOs = sapply(1:nruns, function(n) RES[[n]]@llhood)
  best = which.max(ELBOs)
  return(RES[[best]])
}

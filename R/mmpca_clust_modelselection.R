#' @title Model selection for MMPCA
#' @description A wrapper on \code{\link{mmpca_clust}}() to perform model
#'   selection with an Integrated Classification Likelihood (ICL) criterion.
#'
#' @param dtm an NxV \code{\link{DocumentTermMatrix}} with term-frequency
#'   weighting.
#' @param Qs The vector of clusters to be tested.
#' @param Ks The number of topics to be tested.
#' @param Yinit Parameter for the initialization of Y. It can be either:
#'   \itemize{ \item a string or a function specifying the initialization
#'   procedure. It should be one of ('random', 'kmeans_lda'). See
#'   \code{\link{benchmarks-functions}} functions for more details. \item (Only
#'   when Qs is a singleton) A vector of length N with Q modalities, specifying
#'   the initialization clustering. }
#' @param method The clustering algorithm to be used. Only "BBCVEM" is available
#'   : it corresponds to the branch and bound C-VEM of the original article.
#' @param init.beta Parameter for the initialization of the matrix beta. It can
#'   be either: \itemize{ \item a string specifying the initialization
#'   procedure. It should be one of ('random', 'lda'). See
#'   \code{\link{initializeBeta}}() for more details. \item (Only when Ks is a
#'   singleton) A KxV matrix with each row summing to 1.}
#' @param keep The evolution of the bound is tracked every \code{keep}
#'   iteration.
#' @param max.epochs Specifies the maximum number of pass allowed on the whole
#'   dataset.
#' @param verbose verbosity level.
#' @param nruns number of runs of the algorithm for each (K,Q) pair (default to
#'   1) : the run achieving the best evidence lower bound is selected.
#' @param mc.cores The number of CPUs to use when fitting in parallel the
#'   different models. Default is the number of available cores minus 1.
#'
#' @return \itemize{\item An object of class \code{"\linkS4class{mmpcaClust}"}
#'   containing the best selected model. \item A matrix containing the value of
#'   the ICL for each pair (K,Q).}
#'
#' @export
#' @examples
#' \donttest{
#' ## generate data with the BBCmsg
#' simu = simulate_BBC(N = 100, L = 250)
#' ## Define a grid
#' Qs = 5:6
#' Ks = 3:4
#' ## Run model selection with MoMPCA
#' res <- mmpca_clust_modelselect(simu$dtm.full, Qs = Qs, Ks = Ks,
#'                                Yinit = 'kmeans_lda',
#'                                init.beta = 'lda',
#'                                method = 'BBCVEM',
#'                                max.epochs = 7,
#'                                nruns = 2,
#'                                verbose = 1,
#'                                mc.cores = 2)}
mmpca_clust_modelselect <- function(dtm,
                                    Qs,
                                    Ks,
                                    Yinit = 'random',
                                    method = 'BBCVEM',
                                    init.beta = 'lda',
                                    keep = 1L,
                                    max.epochs = 10L,
                                    verbose = 1L,
                                    nruns = 5L,
                                    mc.cores = (detectCores() - 1)) {

  nQ = length(Qs)
  nK = length(Ks)
  if (nQ != 1 & is.vector(Yinit) & (length(Yinit) > 1))
    stop('Yinit must be a function returning a Q-clustering (since Q is allowed to change)')
  if (nK != 1 & is.matrix(init.beta)) stop('init.beta must be a function return a KxV topic matrix (since K is allowed to change).')
  ICL = matrix(-Inf, nK, nQ)
  rownames(ICL) = Ks
  colnames(ICL) = Qs
  res = NULL
  icl = -Inf
  for (k in 1:nK) {
    K = Ks[k]
    if (verbose > 0) message('-- K = ', K, ' initialize Beta...')
    if (!is.matrix(init.beta)) {
      control_lda_init <- methods::new("LDA_VEMcontrol",
                                       estimate.alpha = FALSE,
                                       estimate.beta = TRUE,
                                       alpha = 1,
                                       verbose = 0L,
                                       nstart = 1L
      )
      Beta = initializeBeta(dtm = dtm, init.beta = init.beta, K = K, verbose = 0, control_lda_init = control_lda_init )
      # Beta = Beta / rowSums(Beta)
    } else {
      Beta = init.beta
    }
    for (q in 1:nQ) {
      Q = Qs[q]
      if (verbose > 0) message('----- Q = ', Q)
      temp = mmpca_clust(dtm = dtm,
                         Q = Q,
                         K = K,
                         Yinit = Yinit,
                         method = method,
                         init.beta = Beta,
                         keep = keep,
                         max.epochs = max.epochs,
                         verbose = 0L,
                         nruns = nruns,
                         mc.cores = mc.cores)
      ICL[k, q] = temp@icl
      if (temp@icl > icl ) {
        icl = temp@icl
        res = temp
      }
    }
  }

  return(list(model = res, ICL = ICL))
}

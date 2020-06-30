
#' @include mmpcaClustcontrol_class.R
NULL



#' @title mmpcaClust class
#' @description An S4 class representing a fitted mmpca model.
#' @name mmpcaClust-class
#'
#' @section Objects from the class:
#' Object of class "\code{mmpcaClust}" are returned by
#'   \code{\link{mmpca_clust}}()
#'
#' @slot call A \code{\link{call}} object specifying the call
#' @slot method The method used in the call
#' @slot clustering The final partition found by the algorithm
#' @slot controls An object of class \code{\linkS4class{mmpcaClustcontrol}}
#'   containing the controls used in the VEM algorithm on the aggregated DTM
#'   during the loop. The slots \code{controls@@control_lda_init} where only use
#'   when init.beta == 'lda'.
#' @slot K An integer specifying the number of topics.
#' @slot Q An integer specifying he number of clusters.
#' @slot N An integer specifying the number of observations.
#' @slot V An integer specifying the number of variables.
#' @slot beta The (KxV) topic matrix.
#' @slot gamma A (QxK) matrix containing the variational paramaters of the
#'   variational distribution of each $\\theta_q$ in its rows.
#' @slot lda_algo An object of class "\code{LDA}" (cf.
#'   \code{\linkS4class{TopicModel}}) containing the results of the
#'   \code{\link{LDA}}() function applied to the aggregated DTM,
#'   with control \code{controls@@control_lda_loop}
#' @slot max.epochs The maximum number of pass through the whole dataset in the
#'   algorithm.
#' @slot logLikelihoods A numeric vector containing the evolution of the
#'   variational bound every \code{keep} iteration.
#' @slot keep An integer specifying the . Mostly useful for the plot function.
#' @slot n_epochs The number of pass through the datasets before convergence.
#'   see details
#' @slot llhood The final value of the variational lower bound.
#' @slot Yinit The value of the initial partition.
#' @slot icl The Integrated Classification Likelihood value.
#'
#' @details The BB-CVEM method is the branch & bound greedy procedure proposed in
#'   the original paper of Jouvin et. al. \url{https://arxiv.org/abs/1909.00721}. The number of epochs in the \code{n_epochs} slot
#'   is actually the true number of pass minus 1 (unless \code{max.epochs} was
#'   reached). Indeed, the last pass before convergence does not change either
#'   the bound or the \code{clustering}, hence it is removed of the counter.
#'
#' @importClassesFrom topicmodels LDA
#' @export
setClass("mmpcaClust",
         slots = list(call = "call",
                      method = "character",
                      clustering = "ANY",
                      K = "integer",
                      Q = "integer",
                      N = "integer",
                      V = "integer",
                      beta = "matrix",
                      gamma = "matrix",
                      lda_algo = "LDA_VEM",
                      max.epochs = "integer",
                      logLikelihoods = "vector",
                      keep = "integer",
                      n_epochs = "integer",
                      llhood = "numeric",
                      Yinit = "vector",
                      icl = "numeric",
                      controls = "mmpcaClustcontrol"
                      )
         )

## ****************************
## plot methods

setMethod("show",
          signature = signature(object = "mmpcaClust"),
          definition = function(object) {
            cat("A", class(object), " model with", object@Q, "clusters and", object@K, "topics.\n")
            }
          )




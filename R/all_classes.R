#' @importClassesFrom topicmodels LDA_VEMcontrol OPTcontrol
#'
#' @slot alpha numeric.
#' @slot control_lda_init list
#' @slot control_lda_loop list
setClass(Class = "mmpca",
                  slots = list(control_lda_init = "LDA_VEMcontrol",
                               control_lda_loop = "LDA_VEMcontrol")
                  )

init_mmpca = function(.Object, control_lda_init, control_lda_loop, ...) {

  if (missing(control_lda_init)) {
    control_lda_init <- new("LDA_VEMcontrol",
                            estimate.alpha = FALSE,
                             estimate.beta = TRUE,
                             alpha = 1,
                             verbose = 0L,
                             nstart = 1L
                             )
  }

  if (missing(control_lda_loop)) {
    control_lda_loop <- new("LDA_VEMcontrol",
                            estimate.alpha = FALSE,
                             alpha = 1,
                             estimate.beta = FALSE,
                             initialize = "model",
                             nstart = 1L,
                             verbose = 0L,
                             var = new('OPTcontrol', iter.max = 1e8L, tol = 1e-9),
                             em = new('OPTcontrol', iter.max = 1L)
    )
  }

  return(list(.Object = .Object,
              control_lda_init = control_lda_init,
              control_lda_loop = control_lda_loop,
              ... = ...
              )
         )
}

setMethod(f = "initialize", signature = "mmpca",
          definition = function(.Object, control_lda_init, control_lda_loop, ...) {
            args <- init_mmpca(.Object, control_lda_init, control_lda_loop, ...)
            .Object <- do.call("callNextMethod", args)
            invisible(.Object)
          }
)


#' @importClassesFrom topicmodels LDA
setClass("mmpca_clust",
         contains = c("mmpca"),
         slots = list(
           # call = "call",
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
                      n_epochs = "integer",
                      llhood = "numeric",
                      Yinit = "vector",
                      icl = "numeric"
                      )
         )

setMethod("show",
          signature = signature(object = "mmpca_clust"),
          definition = function(object) {
            cat("A", class(object), "mmpca model with", object@Q, "clusters and", object@K, "topics.\n")
            }
          )

# new('mmpca_clust', call = call('sqrt', 1))

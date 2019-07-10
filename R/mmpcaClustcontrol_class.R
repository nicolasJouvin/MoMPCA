

#' @importClassesFrom topicmodels LDA_VEMcontrol OPTcontrol
#' @title mmpcaClustcontrol
#' @description  An S4 class for \code{\link{mmpca_clust}}(). It is mainly a
#'   wrapper around the class \code{\linkS4class{TopicModelcontrol}}
#'   (specifically: LDA_VEMcontrol).
#' @name mmpcaClustcontrol-class
#'
#' @slot control_lda_init Object of class \code{"LDA_VEMcontrol"}; specifies the
#'   controls of the VEM algorithm used for the initialization of beta.
#' @slot control_lda_loop Object of class \code{"LDA_VEMcontrol"}; specifies the
#'   controls for the VEM algorithm used after a swap in the branch & bound.
#'
#' @export
setClass(Class = "mmpcaClustcontrol",
         slots = list(control_lda_init = "LDA_VEMcontrol",
                      control_lda_loop = "LDA_VEMcontrol")
)

init_mmpca_control = function(.Object, control_lda_init, control_lda_loop, ...) {

  if (missing(control_lda_init)) {
    control_lda_init <- methods::new("LDA_VEMcontrol",
                            estimate.alpha = FALSE,
                            estimate.beta = TRUE,
                            alpha = 1,
                            verbose = 0L,
                            nstart = 1L
    )
  }

  if (missing(control_lda_loop)) {
    control_lda_loop <- methods::new("LDA_VEMcontrol",
                            estimate.alpha = FALSE,
                            alpha = 1,
                            estimate.beta = FALSE,
                            initialize = "model",
                            nstart = 1L,
                            verbose = 0L,
                            var = methods::new('OPTcontrol', iter.max = 1e8L, tol = 1e-9),
                            em = methods::new('OPTcontrol', iter.max = 100L, tol = 1e-4)
    )
  }

  return(list(.Object = .Object,
              control_lda_init = control_lda_init,
              control_lda_loop = control_lda_loop,
              ... = ...
  )
  )
}

setMethod(f = "initialize", signature = "mmpcaClustcontrol",
          definition = function(.Object, control_lda_init, control_lda_loop, ...) {
            args <- init_mmpca_control(.Object, control_lda_init, control_lda_loop, ...)
            .Object <- do.call("callNextMethod", args)
            invisible(.Object)
          }
)

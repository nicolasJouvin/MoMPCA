
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





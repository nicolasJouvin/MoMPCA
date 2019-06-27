
#' convert a dtm from package tm to sparseMatrix from package Matrix
#' without converting it to full matrix.
#'@param dtm a document-term-matrix from package 'tm'
#'
#'@return a sparse dcgMatrix from package Matrix
DTMtoSparse <- function(dtm){

  return(Matrix::sparseMatrix(i = dtm$i, j = dtm$j, x = dtm$v, dimnames = dtm$dimnames, dims = dim(dtm)))
}







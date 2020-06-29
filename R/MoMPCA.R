#' MoMPCA: Greedy clustering of count data through a mixture of multinomial PCA
#'
#' The MMPCA package implements the branch & bound classification-Variational
#' Expectation Maximisation algorithm described in Jouvin et. al.
#' \url{https://arxiv.org/abs/1909.00721}. It enables the clustering of counts
#' data such as document/term matrix modeled as a mixture of multinomial PCA
#' model. Model selection is performed via an approximated form of the
#' Integrated Classification Likelihood. .
#'
#' The main entry point is the \code{\link{mmpca_clust}}() function to perfom
#' the clustering.
#' @docType package
#' @name MoMPCA
NULL

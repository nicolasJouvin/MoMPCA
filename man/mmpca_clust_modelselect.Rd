% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mmpca_clust_modelselection.R
\name{mmpca_clust_modelselect}
\alias{mmpca_clust_modelselect}
\title{Model selection for MMPCA}
\usage{
mmpca_clust_modelselect(
  dtm,
  Qs,
  Ks,
  Yinit = "random",
  method = "BBCVEM",
  init.beta = "lda",
  keep = 1L,
  max.epochs = 10L,
  verbose = 1L,
  nruns = 5L,
  mc.cores = (detectCores() - 1)
)
}
\arguments{
\item{dtm}{an NxV \code{\link{DocumentTermMatrix}} with term-frequency
weighting.}

\item{Qs}{The vector of clusters to be tested.}

\item{Ks}{The number of topics to be tested.}

\item{Yinit}{Parameter for the initialization of Y. It can be either:
\itemize{ \item a string or a function specifying the initialization
procedure. It should be one of ('random', 'kmeans_lda'). See
\code{\link{benchmarks-functions}} functions for more details. \item (Only
when Qs is a singleton) A vector of length N with Q modalities, specifying
the initialization clustering. }}

\item{method}{The clustering algorithm to be used. Only "BBCVEM" is available
: it corresponds to the branch and bound C-VEM of the original article.}

\item{init.beta}{Parameter for the initialization of the matrix beta. It can
be either: \itemize{ \item a string specifying the initialization
procedure. It should be one of ('random', 'lda'). See
\code{\link{initializeBeta}}() for more details. \item (Only when Ks is a
singleton) A KxV matrix with each row summing to 1.}}

\item{keep}{The evolution of the bound is tracked every \code{keep}
iteration.}

\item{max.epochs}{Specifies the maximum number of pass allowed on the whole
dataset.}

\item{verbose}{verbosity level.}

\item{nruns}{number of runs of the algorithm for each (K,Q) pair (default to
1) : the run achieving the best evidence lower bound is selected.}

\item{mc.cores}{The number of CPUs to use when fitting in parallel the
different models. Default is the number of available cores minus 1.}
}
\value{
\itemize{\item An object of class \code{"\linkS4class{mmpcaClust}"}
  containing the best selected model. \item A matrix containing the value of
  the ICL for each pair (K,Q).}
}
\description{
A wrapper on \code{\link{mmpca_clust}}() to perform model
  selection with an Integrated Classification Likelihood (ICL) criterion.
}
\examples{
\donttest{
## generate data with the BBCmsg
simu = simulate_BBC(N = 100, L = 250)
## Define a grid
Qs = 5:6
Ks = 3:4
## Run model selection with MoMPCA
res <- mmpca_clust_modelselect(simu$dtm.full, Qs = Qs, Ks = Ks,
                               Yinit = 'kmeans_lda',
                               init.beta = 'lda',
                               method = 'BBCVEM',
                               max.epochs = 7,
                               nruns = 2,
                               verbose = 1,
                               mc.cores = 2)}
}

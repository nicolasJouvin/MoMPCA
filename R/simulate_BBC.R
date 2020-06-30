#' @title simulate_BBC
#'
#' @description This function simulate from the MMPCA model with an additional
#'   noise parameter epsilon. The number of cluster is Q=6 for K=4 topics. The
#'   parameter beta is taken to be the row normalized document-term matrix of 4
#'   BBC messages contained in BBCmsg.
#'
#' @param N number of observations.
#' @param L vector of length N containing the total count per observations.
#'   Duplicated if integer.
#' @param epsilon The noise level in the latent space. Quantify how far the
#'   distribution is from theta_true
#' @param lambda A parameter quantifying the class proportion. lambda=1 means
#'   balanced cluster sizes, lower means that the last clusters are bigger, with
#'   an geometric decay in cluster size for the first ones.
#' @param theta_true The true parameter theta for the simulation. If \code{NULL}
#'   (default) then it is initialized to the default value of the experimental
#'   section of the paper.
#'
#' @return A list with names \itemize{ \item \code{dtm.full}: A
#'   \code{\link{DocumentTermMatrix}} object containing the simulated
#'   document-term matrix \item \code{Ytruth}: the simulated partition \item
#'   theta_true The parameter of the simulation }
#' @export
#'
#' @examples
#' simu <- simulate_BBC(N = 100, L = 200, epsilon = 0, lambda = 1)
#' dtm <- simu$dtm.full
#' Ytruth <- simu$Ytruth
#'
simulate_BBC = function(N, L, epsilon = 0, lambda=1, theta_true = NULL){
  Q = 6
  K = 4
  BBCmsg <- MoMPCA::BBCmsg

  if (is.null(theta_true)) {
    theta_true = matrix(0.17, Q, K)

    for (i in 1:K) {
      theta_true[i,i] = 0.5
    }
    theta_true[5,c(1,3)] = 0.33
    theta_true[6,c(2,4)] = 0.33
  } else if (dim(theta_true) != c(6L, 4L)) {
    stop("'theta_true should be a matrix of dimension c(6, 4), currently of dim ",
                dim(theta_true)
         )
  }


  if (length(L) == 1) L <- rep(L, N)

  Pi = lambda^(Q - 0:(Q - 1))
  Pi = Pi / sum(Pi)

  all.clust.present = F
  while (!all.clust.present) {
    Ytruth = apply(stats::rmultinom(N, 1, Pi), 2, which.max)
    if (length(unique(Ytruth)) == Q) all.clust.present = T
  }

  # simulate dtm :
  doc = list()
  Z = matrix(0, N, K)
  for (d in 1:N) {

    doc[[d]] = ""
    for (n in 1:L[d]) {

      test = 1 - stats::rbinom(1, 1, epsilon)
      if (test == 1) {
        Zd = which.max(stats::rmultinom(1, 1, theta_true[Ytruth[d], ]))
      }
      else {
        Zd = sample(1:K, 1)
      }
      Z[d, Zd] = Z[d, Zd] + 1
      msg = switch(Zd,
                  "1" = BBCmsg$msg1,
                  "2" = BBCmsg$msg2,
                  "3" = BBCmsg$msg3,
                  "4" = BBCmsg$msg4
      )

      doc[[d]] = c(doc[[d]], msg[sample(1:length(msg), 1 , replace = TRUE)])
    }
    doc[[d]] = paste(doc[[d]], collapse = " ")
  }

  x = tm::Corpus(tm::VectorSource(doc))

  dtm.full <- tm::DocumentTermMatrix(x, control = list())
  voc.zero =  which(slam::col_sums(dtm.full) == 0)
  if (length(voc.zero) != 0) {
    dtm.full = dtm.full[,-voc.zero]
  }
  #rearange columns in alphabetical order
  dtm.full = dtm.full[,sort(dtm.full$dimnames$Terms)]

  return(list(dtm.full = dtm.full, Ytruth = Ytruth, theta_true = theta_true))
}

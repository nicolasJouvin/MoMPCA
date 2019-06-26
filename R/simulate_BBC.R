#' @title simulate_BBC
#'
#'   This function simulate from the MMPCA model with an additional noise
#'   parameter epsilon.
#'
#' @param N number of observations
#' @param L vector of length N containing the total count per observations.
#'   Duplicated if integer.
#' @param epsilon The noise level in the latent space. Quantify how far the
#'   distribution is from theta_true
#' @param lambda A parameter quantifying the class proportion. lambda=1 means
#'   balanced cluster sizes, lower means that the first cluster are bigger.
#' @param theta_true The true parameter theta for the simulation.
#'
#' @return
#' @export
#'
#' @examples
simulate_BBC = function(N, L, epsilon = 0, lambda=1, theta_true = NULL){

  data("BBCmsg")
  Q = 6
  K = 4

  if (is.null(theta_true)) {
    theta_true = matrix(0.17, Q, K)

    for (i in 1:K) {
      theta_true[i,i] = 0.5
    }
    theta_true[5,c(1,3)] = 0.33
    theta_true[6,c(2,4)] = 0.33
  }


  if (length(L) == 1) L <- rep(L, N)

  Pi = lambda^(Q - 0:(Q - 1))
  Pi = Pi / sum(Pi)

  all.clust.present = F
  while (!all.clust.present) {
    Ytruth = apply(rmultinom(N, 1, Pi), 2, which.max)
    if (length(unique(Ytruth)) == Q) all.clust.present = T
  }

  # simulate dtm :
  doc = list()
  Z = matrix(0, N, K)
  for (d in 1:N) {

    doc[[d]] = ""
    for (n in 1:L[d]) {

      test = 1 - rbinom(1, 1, epsilon)
      if (test == 1) {
        Zd = which.max(rmultinom(1, 1, theta_true[Ytruth[d], ]))
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

  return(list(dtm.full = dtm.full, Ytruth = Ytruth, Z = Z))
}

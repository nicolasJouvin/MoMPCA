bbcvem <- function(dtm,
                   Q,
                   K,
                   model = NULL,
                   Yinit = 'random',
                   method = 'BBCVEM',
                   init.beta = 'lda',
                   keep = 1L,
                   max.epochs = 10L,
                   verbose = 1L) {

  if (as.integer(K) != K || as.integer(K) < 2)
    stop("'K' needs to be an integer of at least 2")
  if (as.integer(Q) != Q || as.integer(Q) < 2)
    stop("'Q' needs to be an integer of at least 2")
  if (as.integer(keep) != keep || keep < 0)
    stop("'keep' needs to be an non-negative integer")

  mycall = match.call()
  if (is.null(model)) {
    model <- methods::new("mmpcaClust")
  }

  control_lda_init = model@controls@control_lda_init
  control_lda_loop = model@controls@control_lda_loop

  if (control_lda_loop@alpha != control_lda_init@alpha)
    warning("Different hyper-parameters alpha for init and loop. Taking the loop one.")

  alpha = control_lda_loop@alpha
  N = dim(dtm)[1]
  V = dim(dtm)[2]

  ## ------------------------------------------------------------
  ## ------------------ Initialisation of Y ---------------------
  ##
  ## Several benchmarks possible.

  if (is.character(Yinit) && length(Yinit) == 1) {
    init_method <- Yinit
    Yinit <- initialize_Y(dtm, Q, K, init = init_method)
  } else if (!is.vector(Yinit) && !is.matrix(Yinit)) {
    stop('Yinit must be either a user-given clustering or a string
         specifying the initialization method.')
  } else if (length(Yinit) != N) {
    stop(paste0("Yinit must be a vector of length ", N,
                ' (length(Yinit) = ', length(Yinit), ').'))
  } else if (length(unique(Yinit)) != Q) {
    stop(paste0("Yinit must be a clustering with ", Q,
                ' groups, not ', length(unique(Yinit)), '.'))
  }

  Y <- Yinit

  ## ------------------------------------------------------------
  ## ----------- Different initialization for beta --------------
  ##
  ## * User defined : in this case init.beta is a user
  ##                  defined init matrix, of dim (KxV),
  ##                  which rows sums to 1
  ##
  ## * 'random' : beta_ij = 1/V + U([0,1e-10])
  ##
  ## * 'lda' : LDA on full dtm. Succession of
  ##           1. LDA with giibs : 5 restarts
  ##           2. LDA with C-VEM
  ##           hyper-parameter alpha is handeled carefully in the controls.
  ##
  ## * 'nmf' : take the basis of the NMF decomposition of the full dtm.

  ## dummy topicmodels::lda object to store beta init
  lda_init = topicmodels::LDA(dtm, k = K,
                              method = 'VEM',
                              control = list(estimate.alpha = FALSE,
                                             estimate.beta = FALSE,
                                             alpha = alpha,
                                             verbose = 0,
                                             nstart = 1,
                                             var = list(iter.max = 1),
                                             em = list(iter.max = 1)
                              )
  )

  ## initialize the lda_init@beta slot
  if (is.matrix(init.beta) ) {
    if (!all.equal(dim(init.beta), c(K, V)) ||
        !all.equal(Matrix::rowSums(init.beta), rep(1, K), tolerance = 1e-12)) {
      stop('init.beta must be a KxV matrix which rows sums to 1.')
      # stop(paste0("K =", K, " V=", V, '\n rowsums = ', Matrix::rowSums(init.beta), ' / diff = ', rep(1, K) - Matrix::rowSums(init.beta), '\n'))
    }
    if (verbose > 0) message('Beta initialisation with a user given beta.')
    lda_init@beta = log(init.beta)
  } else if (is.character(init.beta) && length(init.beta) == 1) {
    beta.init = initializeBeta(dtm, init.beta, K, verbose, control_lda_init)
    lda_init@beta = log(beta.init)
  } else {
    stop('init.beta argument must be a matrix or a string specifying
         the initialization method for beta.')
  }

  ## ------------------------------------------------------------
  ## Compute the bounds with topicmodels::
  ##
  ## Only the var step is done, control
  ## estimate.beta=F ensures the M-step is
  ## locked, hence not re-restimation of beta

  ## Construct meta observations
  P = Matrix::t(Matrix::sparseMatrix(i = 1:N, j = Y, x = rep(1,N)))
  dtm_init = slam::as.simple_triplet_matrix(P %*% DTMtoSparse(dtm))

  ## An LDA to get the individual meta obervations bounds
  ## /!\ We do not reestimate beta. /!\
  lda_aggr = topicmodels::LDA(dtm_init,
                              k = K,
                              control = list(estimate.alpha = FALSE,
                                             estimate.beta = F,
                                             var = list(iter.max = control_lda_loop@var@iter.max,
                                                        tol = control_lda_loop@var@tol),
                                             keep = 1),
                              model = lda_init)


  ## ------------------------------------------------------------
  ## -------------- Branch & Bound with topicmodels -------------
  ##
  ##
  ##   1 loop with beta fixed = B&B
  ##      1 loop over each observation = test each swaps
  ##         1 loop over q = VE-step to obtain the bound modification


  ## Setup quantities used in the B&B
  csize = table(Y)
  PI = csize/N
  lda_algo = lda_aggr
  estimate.beta.algo = control_lda_loop@estimate.beta
  ite_cpt = 0
  # individual meta-observations' contribution to bound
  meta_bounds = rep(0,Q)
  for (q in 1:Q) {
    meta_bounds[q] = lda_algo@loglikelihood[q] + csize[q]*log(PI[q])
  }
  current_bound = sum(meta_bounds)

  ## stock the evolution of the bound every keep iterations
  bounds = c()

  if (verbose > 1) message(' ---- Begin B&B on ', N, ' documents. ----')
  for (epoch in 1:max.epochs) {

    if (verbose > 1) {
      message(' ---------- Epoch number ', epoch, ' :  ')
    }

    # M-step on $\pi$.
    csize = table(Y)
    PI = csize/N

    # optimisation of Y : greedy
    Ycurrent <- Y
    check_bound = current_bound
    for (d in 1:N) {

      ## stock bound every keep iterations
      if (ite_cpt %% keep == 0)
        bounds = c(bounds, current_bound)

      delta_bound = rep(0, Q) # lower bounds change induced by each possible swap
      lda_temp = list() # list of LDA models to stock the bound change for each swap

      # only apply the swap if the cluster is not going to be empty (=> numerical issues)
      if (csize[Y[d]] > 1) {

        for (q in setdiff(1:Q, Y[d])) {
          Y_temp = Y
          Y_temp[d] = q # temporarily swap d to cluster q

          ## temporary M-step on $\pi$
          csize_temp = table(Y_temp)
          Pi_temp = csize_temp/N

          ## Temporary meta-observations. Only construct those impacted by the swap : Y[d] & q
          P_temp = Matrix::t(Matrix::sparseMatrix(i = 1:N, j = Y_temp, x = rep(1,N)) )
          dtm_temp = slam::as.simple_triplet_matrix(P_temp[c(Y[d], q),] %*% DTMtoSparse(dtm))

          ## Update meta-bounds : temporary VE-step
          lda_temp[[q]] <- topicmodels::LDA(dtm_temp,
                                            model = lda_algo,
                                            control = control_lda_loop,
                                            k = K)

          ## Compute the bound difference. Metabounds only change for q and Y[d].
          delta_bound[q] = (lda_temp[[q]]@loglikelihood[1] - meta_bounds[Y[d]]) +
            (lda_temp[[q]]@loglikelihood[2] - meta_bounds[q]) +
            sum(csize_temp[c(Y[d],q)] * log(Pi_temp[c(Y[d], q)]))
        }

        ## ------------------------------------------------------------
        # Select best swap

        best = which.max(delta_bound) # find the best swap

        if (best != Y[d]) { # apply the swap : document d goes into cluster best

          ## M-step on $\pi$
          csize[best] = csize[best] + 1
          csize[Y[d]] = csize[Y[d]] - 1
          PI = csize/N

          if (current_bound > 0)
          {
            stop('Computational error: current bound (log) > 0')
          }

          ## update the meta-observations individual bounds for next iteration
          meta_bounds[Y[d]] =  sum(lda_temp[[best]]@loglikelihood[1])  + csize[Y[d]]*log(PI[Y[d]]) # Y[d] is first meta docs in dtm_temp
          meta_bounds[best] = sum(lda_temp[[best]]@loglikelihood[2])  +  csize[best]*log(PI[best])  # best is second metadocs in dtm_temp
          current_bound = sum(meta_bounds)

          ## apply the swap
          Y[d] = best

          ## Update meta-observations
          P = Matrix::t(Matrix::sparseMatrix(i = 1:N, j = Y, x = rep(1,N)) )
          dtm_algo = slam::as.simple_triplet_matrix(P %*% DTMtoSparse(dtm))

          if (estimate.beta.algo == TRUE) {
            ## re-estimate Beta on aggreated dtm at each swap : unstable and long.
            lda_algo = topicmodels::LDA(dtm_algo,
                                        k = K,
                                        model = lda_algo,
                                        control = control_lda_loop
            )
          }
        }
      }else{
        if (verbose > 1) message(' Can\'t swap observation ', d, ' !')
      }

      ite_cpt = ite_cpt + 1
    }


    #break of outer for loop (of the B&B)
    if (identical(Ycurrent, Y) && check_bound == current_bound)
      break


  }

  ## ------------------------------------------------------------
  ## Re-estimate beta on the aggregated DTM at the end of B&B ---


  # lda_algo = topicmodels::LDA(dtm_algo,
  #                k = K,
  #                method = 'VEM',
  #                model = lda_algo,
  #                control = list(estimate.alpha = FALSE,
  #                               estimate.beta = TRUE,
  #                               alpha = alpha,
  #                               verbose = 0,
  #                               nstart = 1,
  #                               var = list(iter.max = 100)                                                                               )
  #                )

  beta <- exp(lda_algo@beta) / rowSums(exp(lda_algo@beta))


  ## ------------------------------------------------------------
  # Compute final bound

  final_bound = sum(meta_bounds)

  ## ------------------------------------------------------------
  # compute final $\gammas$

  P = Matrix::t(Matrix::sparseMatrix(i = 1:N, j = Y, x = rep(1,N)) )
  dtm_final = slam::as.simple_triplet_matrix(P %*% DTMtoSparse(dtm))

  Vgamma_final <- computeVgamma(lda = lda_algo, dtm.aggr = dtm_final)

  ## ------------------------------------------------------------
  # Compute ICL
  icl <- final_bound - ((K * (V - 1))/2) * log(Q) - ((Q - 1)/2) * log(N)

  res <- methods::new("mmpcaClust",
                      call = mycall,
                      method = method,
                      Yinit = Yinit,
                      clustering = Y,
                      K = as.integer(K),
                      Q = as.integer(Q),
                      N = N,
                      V = V,
                      beta = beta,
                      gamma = Vgamma_final,
                      lda_algo = lda_algo,
                      max.epochs = as.integer(max.epochs),
                      logLikelihoods = bounds[1:(ceiling((epoch - 1) * N)/keep)],
                      keep = as.integer(keep),
                      n_epochs = if (epoch != 1 || epoch < max.epochs) as.integer(epoch - 1) else epoch,
                      llhood = final_bound,
                      icl = icl,
                      controls = model@controls
  )
  res
}


computeVgamma = function(lda, dtm.aggr){

  alpha = lda@alpha
  Q = dim(dtm.aggr)[1]
  K = lda@k
  b <- topicmodels::dtm2ldaformat(dtm.aggr)
  Nq = c();
  for (q in 1:Q) {
    Nq[q] = sum(b$documents[[q]][2,]) # total word count in cluster q
  }

  Egamma = lda@gamma
  Vgamma = matrix(0, Q, K) # true gamma of Blei's paper
  for (q in 1:Q) {
    Vgamma[q,] = Egamma[q,]*(alpha*K + Nq[q])
  }

  Vgamma
}

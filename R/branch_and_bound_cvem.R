mmpca_clust <- function(dtm, Q, K,
                             Yinit='random',
                             init.beta=NULL,
                             keep=-1,
                             control_lda_loop=NULL,
                             control_lda_init=NULL,
                             max.iter=100,
                             verbose=1) {


  if (is.null(control_lda_init)) control_lda_init = default_control_init()
  if (is.null(control_lda_loop)) control_lda_loop = default_control_loop()

  N = dim(dtm)[1]
  # V = dim(dtm)[2]
  alpha = control_lda_loop$alpha

  ###########################################
  #' Initialisation of Y
  #'
  #' Several benchmarks possible.
  ###########################################

  init_method = 'user provided'
  if (isString(Yinit)) {
    init_method <- Yinit
    Yinit <- initialize_Y(dtm, Q, K, init)
  } else if (!is.vector(Yinit) || !is.matrix(Yinit)) {
    stop('Yinit must be either a user-given clustering or a string specifying the initialization method')
  }
  Y <- Yinit


  ###################################################
  #' ------ Different initialization for beta ---------
  #'
  #' * User defined : in this case init.beta is a user defined init matrix
  #'
  #' * 'random' : beta_ij = 1/V + U([0,1e-10])
  #'
  #' * 'lda' : LDA on full dtm. Succession of
  #'           1. LDA with giibs : 5 restarts
  #'           2. LDA with C-VEM
  #'           hyper-parameter alpha is handeled carefully in the controls.
  #'
  #' * 'nmf' : take the basis of the NMF decomposition of the full dtm.
  ###################################################

  ## dummy topicmodels::lda object to store beta init
  lda_init = topicmodel::LDA(dtm, k = K,
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

  ## initialize the la_init@beta slot
  if (is.matrix(init.beta) ) {
    if (verbose > 0) cat('Beta initialisation with custom user given beta.\n')
    lda_init@beta = log(init.beta)
  } else if (isString(init.beta)) {
    if (verbose > 0) cat('Beta initialisation with ',init.beta, '...\n')
    beta.init = initializeBeta(dtm, init.beta, K, verbose, control_gibbs_init, control_lda_init)
    lda_init@beta = log(beta.init)
    if (verbose > 0) cat('... Finished. \n')
  } else {
    stop('init.beta argument must be a matrix or a string specifying the initialization method for beta.')
  }

  ###########################################
  #' Compute the bounds with topicmodels
  #'
  #' Only the var step is done, control
  #' estimate.beta=F ensures the M-step is
  #' locked, hence not re-restimation of beta
  ############################################

  ## Construct meta observations
  P = t(sparseMatrix(i = 1:N, j = Y, x = rep(1,N)) )
  dtm_init = slam::as.simple_triplet_matrix(P %*% DTMtoSparse(dtm))

  ## An LDA to get the individual meta obervations bounds
  ## /!\ We do not reestimate beta. /!\
  lda_aggr = topicmodels::LDA(dtm_init,
                 k = K,
                  control = list(estimate.alpha = FALSE,
                                 estimate.beta = F,
                                 var = list(iter.max = control_lda_loop$var$iter.max,
                                          tol = control_lda_loop$var$tol),
                                 keep = 1),
                 model = lda_init)

  ## should not change but just in case for bound equalities throughout the rest.
  beta.init = t(exp(lda_aggr@beta) / rowSums(exp(lda_aggr@beta)))
  lbeta.init = t(lda_aggr@beta)

  ###################################################
  #' ------ Branch & Bound with topicmodels ---------
  #'
  #'
  #'   1 loop with beta fixed = B&B
  #'      1 loop over each observation = test each swaps
  #'         1 loop over q = VE-step to obtain the bound modification
  ###################################################

  ## Setup quantities used in the B&B
  lda_algo = lda_aggr
  estimate.beta.algo = control_lda_loop$estimate.beta
  lbeta = lbeta.init
  ite_cpt = 0

  ## stock the evolution of the bound every keep iterations
  bounds = c()

  ## for debugging : to keep track of the conservation of sum(Vgammas)
  theoretical.gammaSum = Q*K*alpha + sum(dtm.full)

  if (verbose > 1) cat(' ---- Begin B&B on ', N, ' documents. ----\n')
  for (ite in 1:max.iter) {

    if (verbose > 1) {
      cat('\n ---------- Epoch number ', ite, ' : \n ')
    }

    # M-step on $\pi$.
    csize = table(Y)
    PI = csize/N

    # optimisation of Y : greedy
    Ycurrent <- Y
    for (d in 1:N) {

      delta_bound = rep(0, Q) # lower bounds change induced by each possible swap
      lda_temp = list() # list of LDA models to stock the bound change for each swap

      if (csize[Y[d]] > 1) { # only apply the swap if the cluster is not going to be empty (=> numerical issues)

        for (q in setdiff(1:Q, Y[d])) {
          Y_temp = Y
          Y_temp[d] = q # temporarily swap d to cluster q

          ## temporary M-step on $\pi$
          csize_temp = table(Y_temp)
          Pi_temp = csize_temp/N

          ## Temporary meta-observations. Only construct those impacted by the swap : Y[d] & q
          P_temp = t(sparseMatrix(i = 1:N, j = Y_temp, x = rep(1,N)) )
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

        #######################
        # Select best swap
        #######################
        best = which.max(delta_bound) # find the best swap

        if (best != Y[d]) { # apply the swap : document d goes into cluster best

          ## M-step
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
          P = t(sparseMatrix(i = 1:N, j = Y, x = rep(1,N)) )
          dtm_algo = slam::as.simple_triplet_matrix(P %*% DTMtoSparse(dtm))

          if (estimate.beta.algo == TRUE) {
            # re-estimate Beta : unstable and long.
            lda_algo = topicmodels::LDA(dtm_algo, model = lda_algo, control = control_lda_init, k = K)
            lbeta  = t(lda_algo@beta)
          }



        }
      }else{
        if (verbose > 0) cat('\n Can\'t swap doc ', d, ' !\n')
      }

      ## stock bound every keep iterations
      ite_cpt = ite_cpt + 1
      if (ite_cpt %% keep) bounds = c(bounds, current_bound)
    }

    if (identical(Ycurrent, Y)) {
      #break of outer for loop (of the B&B)
      break
    }
  }

  ###################################################################
  ## Re-estimate beta on the aggregated DTM at the end of B&B ---
  ###################################################################

  lda_algo = topicmodels::LDA(dtm_algo,
                 k = K,
                 method = 'VEM',
                 model = lda_algo,
                 control = list(estimate.alpha = FALSE,
                                estimate.beta = TRUE,
                                alpha = alpha,
                                verbose = 0,
                                nstart = 1,
                                var = list(iter.max = 100, tol = fixedpoint.tol)                                                                               )
                 )
  lbeta = t(lda_algo@beta)

  ## Check for conservation of sum(Vgammas)
  if (debug) {
    Vgamma = computeVgamma(lda = lda_algo, dtm.aggr = dtm_algo)
    gammaSum = sum(Vgamma)
    if (gammaSum != theoretical.gammaSum )
      warning('Problem : total sum of pseudo counts not conserved !')
  }



  #######################
  # Compute final bound
  #######################

  final_bound = sum(meta_bounds)

  ###########################
  # compute final $\gammas$
  ###########################


  P = t(sparseMatrix(i = 1:N, j = Y, x = rep(1,N)) )
  dtm_final = slam::as.simple_triplet_matrix(P %*% DTMtoSparse(dtm))

  b <- topicmodels::dtm2ldaformat(dtm_final)
  Nq = c(); NVq = c()
  for (q in 1:Q) {
    Nq[q] = sum(b$documents[[q]][2,]) # total word count in cluster q
    NVq[q] = length(b$documents[[q]][1,]) # unique
  }

  Egama = lda_algo@gamma
  Vgamma_final = matrix(0, Q, K) # true gamma of Blei's paper
  for (q in 1:Q) {
    Vgamma_final[q,] = Egama[q,]*(alpha*K + Nq[q])
  }

  return(list(Yinit = Yinit,
              init_method = init_method,
              Y = Y,
              lda_init = lda_init,
              lda_final = lda_algo,
              final_bound = final_bound,
              bounds = bounds ,
              ARIs = ARIs,
              gamma = Vgamma_final,
              max.iter = max.iter,
              n_iterations = ite,
              alpha = alpha,
              lda_controls = list(init = control_lda_init, loop = control_lda_loop),
              Q = Q,
              K = K)
         )
}


computeVgamma = function(lda, dtm.aggr){

  alpha = lda@alpha
  Q = dim(dtm.aggr)[1]
  K = lda@k
  b <- dtm2ldaformat(dtm.aggr)
  Nq = c();
  for (q in 1:Q) {
    Nq[q] = sum(b$documents[[q]][2,]) # total word count in cluster q
  }

  Egamma = lda@gamma
  Vgamma = matrix(0, Q, K) # true gamma of Blei's paper
  for (q in 1:Q) {
    Vgamma[q,] = Egamma[q,]*(alpha*K + Nq[q])
  }

  return(Vgamma)
}

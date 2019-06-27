plot_topics <- function(res, ...) {

  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("dplyr", quietly = TRUE) ||
      !requireNamespace("tidytext", quietly = TRUE)) {
    stop("Packages : ggplot2, dplyr and tidytext are necessary to plot topics.
         Please check that they are all installed or consider installing them.")
  }

  entropy_corection_lda <- function(x, s=2){

    y = x/sum(x)
    K = length(x)
    ent = 1 + sum(y*log(y))/log(K)
    spec = x * ent^s
    return(spec)
  }

  args <- list(...)
  n_words <- if ("n_words" %in% names(args)) args$n_words else 10
  s <- if ("s" %in% names(args)) args$s else 2

  rep_topics <- tidytext::tidy(res@lda_algo, matrix = "beta")

  # entropic smoothing of topic probability
  rep_topics_cor = rep_topics %>%
    dplyr::group_by(term) %>%
    dplyr::mutate(beta = entropy_corection_lda(beta, s)) %>%
    dplyr::ungroup()


  rep_top_terms <- rep_topics_cor %>%
    dplyr::group_by(topic) %>%
    dplyr::top_n(n = n_words, beta) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(topic, -beta)

  gg = rep_top_terms  %>%
    dplyr::mutate(term = stats::reorder(term, beta)) %>%
    ggplot2::ggplot(ggplot2::aes(term, beta, fill = factor(topic))) +
    ggplot2::geom_col(show.legend = FALSE) +
    ggplot2::facet_wrap(~ topic, scales = "free") +
    ggplot2::coord_flip() +
    ggplot2::ylab("Topic matrix")

  gg
}


plot_bound <- function(res, ...) {
  if (is.null(res@logLikelihoods)) {
    stop("Bound evolution was not tracked, try restarting with the 'keep'
         argument set to a useful value.")
  }

  size = res@n_epochs*res@N
  if (res@keep > size) {
    print("The number of iteration is lower than 'keep'. Hence, the bound evolution
          can't be plotted...")
  } else if (size %% res@keep == 0) {
    df <- data.frame(iteration = seq(1, size, by = res@keep),
                   bound = res@logLikelihoods)
  } else {
    df <- data.frame(iteration = head(seq(1, size, by = res@keep), -1),
                     bound = res@logLikelihoods)
  }

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    gg <- ggplot2::ggplot(df, ggplot2::aes(x = iteration, y = bound)) +
      ggplot2::geom_point(size = 0.6) +
      ggplot2::geom_hline( yintercept = max(df$bound), show.legend = T,
                           col = 'red', alpha = 0.5) +
      ggplot2::scale_x_continuous(breaks = seq(0, size - 1, by = res@N),
                                  labels = paste0('Epoch n0 : ', 1:res@n_epochs)
                                  ) +
      ggplot2::xlab('') +
      ggplot2::ylab('Classification Evidence Lower Bound')
    return(gg)
  } else {
    plot(df)
  }
}

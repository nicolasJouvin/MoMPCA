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

plot_bound <- function(res) {
  if (is.null(res@logLikelihoods)) {
    stop("Bound evolution was not tracked, try restarting with the 'keep'
         argument set to a useful value.")
  }
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    df <- 1
  }

  df
}

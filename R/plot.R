## For CRAN check ...
## https://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when/12429344#12429344
utils::globalVariables(names = c('term', 'topic'))

#' @include mmpcaClust_class.R
NULL

#' @title Plot function for object mmpcaClust
#' @description Use ggplot2 if available.
#' @param x an S4 object of class \code{\linkS4class{mmpcaClust}}
#' @param type Either: \itemize{\item 'topics' (default): Show the top topic
#'   words of topic matrix. See \code{\link{plot_topics}} documentation for more details. \item 'bound': plot the lower bound evolution during
#'   the greedy procedure. See \code{\link{plot_bound}} documention for more details.}
#' @param ... optional argument specifying the number of words to display and
#'   the entropy correction to apply when calling \code{plot_topics}().
#' @return a plot
#' @export
setMethod(f = "plot",
          signature = signature('mmpcaClust', "missing"),
          definition = function(x, type = "topics", ...) {
            switch(type,
                   "topics" = {
                     plot_topics(x, ...)
                   },
                   "bound" = {
                     plot_bound(x)
                   },
                   stop('This plot type does not exists.')
            )
          }
)

#' @title plot_topics
#' @description Plot topic matrix
#'
#' @param res An S4 object of class \code{\linkS4class{mmpcaClust}}
#' @param s an entropy correction parameter for the topic matrix. It is applied
#'   to the beta matrix before sorting the words by highest probability. The
#'   greater, the more emphasis is put towards words contributing a lot to the
#'   entropy of a topic. Set s=1 to ignore.
#' @param n_words the number of words to display per topic.
#'
#' @return a ggplot2 object
plot_topics <- function(res, s=2, n_words = 10) {

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

  # args <- list(...)
  # n_words <- if ("n_words" %in% names(args)) args$n_words else 10
  # s <- if ("s" %in% names(args)) args$s else 2

  rep_topics <- tidytext::tidy(res@lda_algo, matrix = "beta")

  # entropic smoothing of topic probability
  rep_topics_cor = rep_topics %>%
    dplyr::group_by(`term`) %>%
    dplyr::mutate(beta = entropy_corection_lda(beta, s)) %>%
    dplyr::ungroup()


  rep_top_terms <- rep_topics_cor %>%
    dplyr::group_by(`topic`) %>%
    dplyr::top_n(n = n_words, beta) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(`topic`, -beta)

  fill_ggplot <- factor(rep_top_terms$topic)
  gg = rep_top_terms  %>%
    dplyr::mutate(term = stats::reorder(`term`, beta)) %>%
    ggplot2::ggplot(ggplot2::aes(`term`, `beta`, fill = fill_ggplot)) +
    ggplot2::geom_col(show.legend = FALSE) +
    ggplot2::facet_wrap(~ `topic`, scales = "free") +
    ggplot2::coord_flip() +
    ggplot2::ylab("Topic matrix")

  gg
}

#' @title Bound evolution plot
#' @description Plot lower bound evolution
#'
#' @param res An S4 object of class \code{\linkS4class{mmpcaClust}}
#'
#' @return a ggplot2 object if ggplot2 is available. Plot on the device
#'   otherwise.
plot_bound <- function(res) {
  if (is.null(res@logLikelihoods)) {
    stop("Bound evolution was not tracked, try restarting with the 'keep'
         argument set to a useful value.")
  }

  size = res@n_epochs*res@N
  if (res@keep > size) {
    message("The number of iteration is lower than 'keep'. Hence, the bound evolution
          can't be plotted...")
  } else if (size %% res@keep == 0) {
    df <- data.frame(iteration = seq(1, size, by = res@keep),
                   bound = res@logLikelihoods)
  } else {
    df <- data.frame(iteration = utils::head(seq(1, size, by = res@keep), -1),
                     bound = res@logLikelihoods)
  }

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    gg <- df %>%
      ggplot2::ggplot(ggplot2::aes_string(x = 'iteration', y = 'bound')) +
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

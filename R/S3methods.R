#### S3-Methods ----

#' Print an object of class hyponetResult
#'
#' @param x An object of class hyponetResult
#' @param ... Currently not used
#'
#' @export
print.hyponetResult <- function(x, ...) {
  print(x$fit)
}

#' Summarize an object of class hyponetResult
#'
#' @param object An object of class hyponetResult
#' @param ... Currently not used
#'
#' @export
summary.hyponetResult <- function(object, ...) {

  app <- object$fit[, grep('_app$', colnames(object$fit)), drop = FALSE]
  colnames(app) <- gsub('_app$', '', colnames(app))
  con <- object$fit[, grep('_con$', colnames(object$fit)), drop = FALSE]
  colnames(con) <- gsub('_con$', '', colnames(con))

  targets <- colnames(object$data)[object$target]
  if (is.null(targets)) targets <- object$target
  if (length(targets) == ncol(object$data)) targets <- NULL

  out <- list(app = app, con = con,
    resample = object$resample,
    nrep = object$nrep,
    hypothesis = object$hypothesis,
    target = targets,
    tolerance = object$tolerance)

  class(out) <- 'summary.hyponetResult'
  return(out)
}

#' Print the summary of a hyponetResult
#'
#' @param x An object of class summary.hyponetResult
#' @param ... Currently not used
#'
#' @export
print.summary.hyponetResult <- function(x, ...) {

  cat('\nSummary of the network hypothesis test\n')
  cat('======================================\n')

  if (is.null(x$target)) {
    cat('\nHypothesis type: global\n')
  } else {
    cat(paste0('\nHypothesis type: ', x$hypothesis, '\n'))
    cat(paste0('Nodes in the target system: ', paste(x$target, collapse = ', '), '\n'))
  }
  cat(paste0('Crud tolerance: ', x$tolerance, '\n'))
  cat(paste0('Results based on ', x$nrep, ifelse(x$resample == 'exchangeableSampling', ' exchangeable', ' bootstrap'), ' samples.\n'))

  cat('\nConstraint fit statistics:\n')
  print(round(x$con, 4))

  cat('\nApproximation fit statistics:\n')
  print(round(x$app, 4))

  if (x$hypothesis == 'embed') {
    cat('\nNOTE: For embedding hypotheses, the chi-square and (approximation) RMSEA cannot be computed. They are reported for the full network instead.\n')
  }

}


#' Plot the results of a network hypothesis test
#'
#' @param x An object of class hyponetResult
#' @param which Which of the fit statistics to plot. Either 'all' (the default),
#' a vector of names of the fit statistics, or a numeric vector indicating which
#' of the fit values provided in the hyponetResult object to plot.
#' @param ask A logical indicating whether to ask before drawing the next plot
#' (if more than one fit statistic is being plotted).
#' @param breaks A breaks argument as described in \code{?hist}. Defaults to
#' 'FD' for the use of \code{nclass.FD}.
#' @param mfrow A vector of two numbers setting the number of rows and columns
#' to be displayed when plotting multiple fit statistics.
#' @param ... Other arguments to be passed to \code{hist}.
#'
#' @export
plot.hyponetResult <- function(x,
  which = 'all',
  ask = FALSE,
  breaks = 'FD',
  mfrow = NULL,
  ...) {

  originalFit <- as.matrix(x$originalFit)
  resampleFit <- as.matrix(x$resampleFit)

  available <- colnames(resampleFit)

  if (is.numeric(which)) {
    which <- available[which]
  } else if (is.character(which)) {
    if (identical(which, 'all')) {
      which <- available
    }
    which <- unique(which)
    which <- match.arg(which, choices = available, several.ok = TRUE)
  }

  if (is.null(mfrow) && length(which) > 1) {
    nr <- ceiling(sqrt(length(which)))
    nc <- ceiling(length(which) / nr)
    mfrow <- c(nr, nc)
  }

  if (length(which) > 1 && !is.null(mfrow)) {
    oldPar <- graphics::par(mfrow = mfrow)
    on.exit(graphics::par(oldPar), add = TRUE)
  } else {
    oldAsk <- grDevices::devAskNewPage(ask = ask && length(which) > 1)
    on.exit(grDevices::devAskNewPage(oldAsk), add = TRUE)
  }

  for (fit in which) {
    plot_fit_histogram(
      resampled = resampleFit[, fit],
      original = originalFit[1, fit],
      fit = fit,
      breaks = breaks,
      m = x$nreps,
      resample = x$resample,
      ...
    )
  }

  invisible(x)
}

# Helper function for a single histogram
plot_fit_histogram <- function(resampled, original, fit, breaks = 'FD',
  m, resample,
  ...) {

  tit <- strsplit(fit, '_')[[1]]
  tit[2] <- ifelse(tit[2] == 'app', '(Approximation)', '(Constraint)')
  graphics::hist(
    resampled,
    breaks = breaks,
    main = 'Fit value distribution',
    xlab = tit[1],
    sub = tit[2],
    ...
  )

  graphics::abline(
    v = original,
    lty = 2,
    lwd = 2
  )
}

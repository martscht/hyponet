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

  app <- object$fit[, grep('_app$', colnames(object$fit))]
  colnames(app) <- gsub('_app$', '', colnames(app))
  con <- object$fit[, grep('_con$', colnames(object$fit))]
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
  cat(paste0('Tolerance for fuzzy testing: ', x$tolerance, '\n'))
  cat(paste0('Results based on ', x$nrep, ifelse(x$resample == 'exchangeableSampling', ' exchangeable', ' bootstrap'), ' Samples.\n'))

  cat('\nConstraint fit statistics:\n')
  print(round(x$con, 4))

  cat('\nApproximation fit statistics:\n')
  print(round(x$app, 4))

  if (x$hypothesis == 'embed') {
    cat('\nNOTE: For embedding hypotheses, the chi-square and (approximation) RMSEA cannot be computed. They are reported for the full network instead.\n')
  }

}

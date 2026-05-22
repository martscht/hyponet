#' (Local) Hypothesis Tests for GGMs
#'
#' This function tests the null hypothesis that the proposed adjacency matrix
#' holds in the sample data via resampling approaches.
#'
#' @param data The data to be tested. Can be either raw data (provided as
#' a matrix or data frame) or a \code{hyponetData} object obtained from
#' \code{exchangeableSampling} or \code{parametricBootstrap}.
#' @param adjacency The proposed adjacency matrix. Rows and columns must be in
#' the same order as the variables in the supplied data.
#' @param target The column numbers of the variables that constitute the target
#' system.
#' @param fits The fit statistics to be computed. Can be either 'all' (the
#' default) or a character vector of names. See details.
#' @param tolerance The absolute tolerance value for fuzzy hypothesis testing.
#' See details.
#' @param resampleType The type of resampling to use. Ignored if a
#' \code{hypoeNet} object is provided to \code{data}.
#' @param ... Other arguments to be passed to \code{exchangeablSampling} or
#' \code{parametricBootstrap} when providing raw data.
#'
#' @return An array of test statistics with their corresponding p-Values.
#'
#' @details
#' The currently available fit statistics can be listed with
#' \code{hyponet:::implementedFits()}. This returns a vector of all fit
#' statistics' names and the internal functions by which they are computed. Any
#' combination of these fit statistics can used simultaneously.
#'
#' The \code{toelarance} value refers to the absolute discrepancy between the
#' model-implied and the empirical network edges. Any differences that fall
#' below this threshold are not included in the computation of the fit
#' statistics.
#'
#' @export

hypothesisTest <- function(data, adjacency, target = NULL,
  fits = 'all', tolerance = 1e-6,
  resampleType = c('exchangeableSampling', 'parametricBootstrap'),
  ...) {

  # If resampling was not done ahead of time, do it now
  if (!inherits(data, 'hyponetData')) {
    if (is.null(target)) target <- 1:ncol(data)

    resampleType <- match.arg(resampleType)

    data <- switch(resampleType,
      exchangeableSampling = exchangeableSampling(
        data = data,
        adjacency = adjacency,
        target = target, ...),
      parametricBootstrap = parametricBootstrap(
        data = data,
        adjacency = adjacency,
        ...))
  }

  if (is.null(target)) target <- 1:ncol(data$original)
  if (data$type == 'exchangeableSampling' & !identical(data$target, target)) {
    warning('The target systems used during resampling and during testing do not match. Results may be severely biased.')
  }

  # determine fit statistics for original data (standardized)
  originalFit <- determineFits(data = data$original, adjacency = adjacency,
    target = target, fits = fits, tolerance = tolerance)

  # determine fit statistics for copies
  resampleFit <- determineFits(data = data$copies, adjacency = adjacency,
    target = target, fits = fits, tolerance = tolerance)

  # determine p-value
  larger <- apply(resampleFit, 1, function(x) x > originalFit) |> rowSums()
  smaller <- apply(resampleFit, 1, function(x) x < originalFit) |> rowSums()
  equal <- apply(resampleFit, 1, function(x) x == originalFit) |> rowSums()

  pvalues <- (larger + .5*equal) / (larger + smaller + equal)

  out <- rbind(originalFit, pvalues)
  rownames(out) <- c('Statistic', 'p-Value')

  # class(out) <- 'hyponetResult'

  return(out)
}

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
#' @param hypothesis The type of hypothesis to be tested. Can be either 'local'
#' (the default) or 'embed'. See details.
#' @param fits The fit statistics to be computed. Can be either 'all' (the
#' default) or a character vector of names. See details.
#' @param tolerance The absolute tolerance value for fuzzy hypothesis testing.
#' See details.
#' @param resample The type of resampling to use. Ignored if a
#' \code{hypoeNet} object is provided to \code{data}.
#' @param ... Other arguments to be passed to \code{exchangeablSampling} or
#' \code{parametricBootstrap} when providing raw data.
#'
#' @return An object of class \code{hyponetResult}. Specifically a list
#' containing \code{fit} (an array of fit measures and their p-values),
#' \code{resample} (the resampling type used), \code{hypothesis} (the hypothesis
#' type specified), \code{target} (the nodes contained in the target system),
#' and \code{tolerance} (the tolerance value used for fuzzy hypothesis testing).
#'
#' @details
#' The function can be used to test two types of hypotheses. Using
#' \code{hypothesis = 'local'} results in a test of all edges within the target
#' system #' (provided via \code{target}). Using \code{hypothesis = 'embed'}
#' results in a test #' of all edges related to nodes in the target system. Thus
#' the latter only #' excludes edges from the test which connect non-target
#' nodes to each other.
#'
#' Currently available fit statistics can be listed with
#' \code{hyponet:::implementedFits()}. This returns a vector of all fit
#' statistics' names and the internal functions by which they are computed. Any
#' combination of these fit statistics can used simultaneously.
#'
#' The \code{tolarance} value refers to the absolute discrepancy between the
#' model-implied and the empirical network edges. Any differences that fall
#' below this threshold are not included in the computation of the fit
#' statistics.
#'
#' @export

hypothesisTest <- function(data, adjacency,
  target = NULL, hypothesis = c('local', 'embed'),
  fits = 'all', tolerance = 1e-6,
  resample = c('exchangeableSampling', 'parametricBootstrap'),
  ...) {

  call <- match.call()

  # check which kind of hypothesis is being tested
  hypothesis <- match.arg(hypothesis)

  # If resampling was not done ahead of time, do it now
  if (!inherits(data, 'hyponetData')) {
    if (is.null(target)) target <- 1:ncol(data)

    resample <- match.arg(resample)

    data <- switch(resample,
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
  if (data$resample == 'exchangeableSampling' & !identical(data$target, target)) {
    warning('The target systems used during resampling and during testing do not match. Results may be severely biased.')
  }

  # determine fit statistics for original data (standardized)
  originalFit <- determineFits(data = data$original, adjacency = adjacency,
    target = target, fits = fits, tolerance = tolerance, hypothesis = hypothesis)

  # determine fit statistics for copies
  resampleFit <- determineFits(data = data$copies, adjacency = adjacency,
    target = target, fits = fits, tolerance = tolerance, hypothesis = hypothesis)

  # determine p-value
  larger <- sweep(resampleFit, 2, originalFit, `>`) |> colSums()
  smaller <- sweep(resampleFit, 2, originalFit, `<`) |> colSums()
  equal <- sweep(resampleFit, 2, originalFit, `==`) |> colSums()

  pvalues <- (larger + .5*equal) / (larger + smaller + equal)

  results <- rbind(originalFit, pvalues)
  rownames(results) <- c('Statistic', 'p-Value')

  out <- list(fit = results,
    resample = data$resample,
    hypothesis = hypothesis,
    target = target,
    tolerance = tolerance,
    nreps = length(data$copies),
    call = call,
    data = data$original)

  class(out) <- 'hyponetResult'

  return(out)
}

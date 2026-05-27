#### Collection of general helper functions ----

#' Convert correlation matrix to partial correlations
#' @noRd
#' @keywords internal
c2pc <- function(R) {

  Omega <- solve(R)
  Rho <- -stats::cov2cor(Omega)
  diag(Rho) <- 1

  return(Rho)
}

#' Determine implied correlation matrix via glasso
#' @noRd
#' @keywords internal
impliedCor <- function(empirical, adjacency, n) {

  # Rephrase constraints
  constr <- which(adjacency == 0, arr.ind = TRUE)
  constr <- constr[constr[, 1] != constr[, 2], ]

  # Estimate via glasso
  Net <- glasso::glasso(s = empirical, rho = 0, nobs = n,
    zero = constr) |> suppressWarnings() |> suppressMessages()
  SigmaHat <- stats::cov2cor(Net$w)

  return(SigmaHat)
}

#' Estimate implied and empirical networks
#' @noRd
#' @keywords internal
fullEstimation <- function(data, adjacency, target = 1:ncol(data)) {

  # determine sample size
  n <- nrow(data)
  # determine p
  p <- ncol(data)

  # get empirical and implied correlation matrices
  emp <- stats::cor(data)
  imp <- impliedCor(emp, adjacency, n)

  # same thing for partial correlation matrices
  empRho <- c2pc(emp)
  impRho <- c2pc(imp)

  # filter down to local version
  # emp <- emp[target, target]
  # imp <- imp[target, target]
  # empRho <- empRho[target, target]
  # impRho <- impRho[target, target]
  # adjacency <- adjacency[target, target]

  # combine output
  out <- list(n = n, p = p,
    emp = emp,
    empRho = empRho,
    imp = imp,
    impRho = impRho,
    adjacency = adjacency,
    target = target)

  class(out) <- 'hyponetEst'

  return(out)
}


#' Check for all available measures
#' @noRd
#' @keywords internal
implementedFits <- function(prefix = 'fit_') {

  ns <- asNamespace('hyponet')

  fitFuns <- ls(ns, all.names = TRUE, pattern = paste0('^', prefix, '[a-zA-Z]+$'))

  fits <- sub(paste0('^', prefix), '', fitFuns)

  stats::setNames(fitFuns, fits)

}

#' Determine (any collection of) fit measures for a single application
#' @noRd
#' @keywords internal
determineSingleFits <- function(est, tolerance, fits, hypothesis) {

  # check whether est is the right kind of object
  if (!inherits(est, 'hyponetEst')) {
    stop('Unable to determine fit measures.', call. = FALSE)
  }

  # Check which fit measures are implemented
  implemented <- implementedFits()
  available <- names(implemented)

  # check whether only character values were provided
  if (!is.character(fits)) {
    stop("'fits' must be a character vector or 'all'.", call. = FALSE)
  }

  # check whether they should all be computed
  if (identical(fits, 'all')) {
    fits <- available
  }

  # discard capitalization
  fits <- tolower(fits)

  # respond if something unknown was provided
  unknown <- setdiff(fits, available)
  if (length(unknown) > 0) {
    stop(
      'Unknown fit measure(s): ',
      paste(unknown, collapse = ', '),
      '\nAvailable fit measures are: ',
      paste(available, collapse = ', '),
      call. = FALSE
    )
  }

  ns <- asNamespace("hyponet")

  out <- vapply(fits, function(fit) {
    fun <- get(implemented[[fit]], envir = ns, inherits = FALSE)
    fun(est = est, tolerance = tolerance, hypothesis = hypothesis)
  },
    numeric(1)
  )

  return(out)
}

#' Determine (any collection of) fit measures for many applications
#' @noRd
#' @keywords internal
determineFits <- function(data, adjacency, target = NULL, hypothesis, fits = 'all', tolerance = 1e-6) {

  if (inherits(data, 'matrix')) {
    data <- list(data)
  }
  if (is.null(target)) target <- 1:ncol(data[[1]])

  estimated <- lapply(data, fullEstimation, adjacency = adjacency, target = target)
  fitted <- lapply(estimated, determineSingleFits, tolerance = tolerance, fits = fits, hypothesis = hypothesis)
  fitted <- do.call(rbind, fitted)

  return(fitted)
}

#' Create a filter matrix for the relevant edges
#' @noRd
#' @keywords internal
fitFilter <- function(est, tolerance = 1e-6,
  hypothesis = c("local", "embed"),
  triangle = FALSE) {

  hypothesis <- match.arg(hypothesis)

  p <- est$p
  target <- est$target
  adjacency <- est$adjacency

  # tolerance filter for fuzziness
  filt <- abs(est$empRho - est$impRho) > tolerance

  # only restricted / missing edges can contribute
  filt <- filt & adjacency == 0

  # remove diagonal
  diag(filt) <- FALSE

  if (hypothesis == "local") {

    # only pairs inside the target system
    local <- matrix(FALSE, p, p)
    local[target, target] <- TRUE
    diag(local) <- FALSE

    filt <- filt & local
  }

  if (hypothesis == "embed") {

    # rows of target nodes: target nodes as response variables
    embed <- matrix(FALSE, p, p)
    embed[target, ] <- TRUE
    diag(embed) <- FALSE

    filt <- filt & embed
  }

  # for symmetric statistics
  if (triangle) {
    filt <- filt & lower.tri(filt)
  }

  return(filt)
}


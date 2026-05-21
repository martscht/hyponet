#' Parametric Bootstrap for GGMs
#'
#' This is a function to generate parametric bootstrap samples of a dataset,
#' given a proposed adjacency matrix.
#'
#' @param data A data.frame or matrix containing all variables to include in the
#' network.
#' @param adjacency An adjacency matrix indicating whether an edge is assumed (1)
#' or not assumed (0) between two variables. Rows and columns must be in the
#' same order as the columns of \code{data}.
#' @param m The number of bootstrap samples to generate.
#' @param ncores The number of cores used in parallel computing.
#'
#' @return A list containing the (standardized) original dataset and a list of
#' all generated copies.
#'
#' @export

parametricBootstrap <- function(data, adjacency, m = 500, ncores = 1) {

  # prepare data
  if (inherits(data, 'data.frame')) {
    if (any(!sapply(data, is.numeric))) {
      stop('All variables in the data must be numeric.')
    }
  }
  data <- scale(as.matrix(data))
  attr(data, 'scaled:center') <- NULL
  attr(data, 'scaled:scale') <- NULL

  # compute implied partial correlations
  S <- stats::cor(data)
  n <- nrow(data)
  mu <- colMeans(data)

  # Fit network with adjacencies
  constr <- which(adjacency == 0, arr.ind = TRUE)
  constr <- constr[constr[, 1] != constr[, 2], ]
  Net <- glasso::glasso(s = S, rho = 0, nobs = n,
    zero = constr) |> suppressWarnings() |> suppressMessages()

  SigmaHat <- stats::cov2cor(Net$w)

  # Create copies object
  copies <- vector('list', m)

  # Parallelization setup
  cl <- parallel::makeCluster(ncores)

  on.exit(parallel::stopCluster(cl), add = TRUE)

  parallel::clusterExport(
    cl,
    varlist = c('SigmaHat', 'n'),
    envir = environment()
  )

  # fill copies object
  copies <- pbapply::pblapply(seq_len(m), cl = cl, FUN = function(mi) {
    mvtnorm::rmvnorm(n = n, sigma = SigmaHat)
  })

  # Generate output
  out <- list(original = data, copies = copies)
  class(out) <- 'hyponetData'

  return(out)

}

#' Exchangeable Resampling for GGMs
#'
#' This is a function to generate exchangeable copies of a dataset based on the
#' exchangeable sampling algorithm propose by Lin et al. (2025).
#'
#' @param data A data.frame or matrix containing all variables to include in the
#' network.
#' @param adjacency An adjacency matrix indicating whether an edge is assumed (1)
#' or not assumed (0) between two variables. Rows and columns must be in the
#' same order as the columns of \code{data}.
#' @param target The columns of \code{data} which are nodes of the target system.
#' @param m The number of exchangeable copy datasets to generate.
#' @param l The number of passes of the markov chain.
#' @param ncores The number of cores used in parallel computing.
#'
#' @return A list containing the (standardized) original dataset, the hub data, and a list of
#' all generated copies.
#'
#' @export

exchangeableSampling <- function(data, adjacency, target = 1:ncol(data),
  m = 500, l = 3, ncores = 1) {

  # prepare data
  if (inherits(data, 'data.frame')) {
    if (any(!sapply(data, is.numeric))) {
      stop('All variables in the data must be numeric.')
    }
    data <- scale(as.matrix(data))
    attr(data, 'scaled:center') <- NULL
    attr(data, 'scaled:scale') <- NULL
  }

  # Generate hub
  Xhub <- data
  for (li in seq_len(l)) {
    Xhub <- algo2step1(Xhub, target, adjacency)
  }

  # Create copies object
  Xtilde <- vector('list', m)

  # Parallelization setup
  cl <- parallel::makeCluster(ncores)

  on.exit(parallel::stopCluster(cl), add = TRUE)

  parallel::clusterExport(
    cl,
    varlist = c('Xhub', 'target', 'adjacency', 'l', 'algo2step1', 'algo1', 'orthVec'),
    envir = environment()
  )

  # fill copies object
  Xtilde <- pbapply::pblapply(seq_len(m), cl = cl, FUN = function(mi) {
    tmp <- Xhub

    for (li in seq_len(l)) {
      tmp <- algo2step1(tmp, rev(target), adjacency)
    }

    tmp
  })

  # Generate output
  out <- list(original = data, hub = Xhub, copies = Xtilde)
  class(out) <- 'hyponetData'

  return(out)
}

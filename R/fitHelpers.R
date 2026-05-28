# Collection of helper functions for determining fit statistics

#' Determine the information matrix used in computation of modindices
#' @noRd
#' @keywords internal
infoMatrix <- function(imp, params, n) {
  q <- nrow(params)
  Info <- matrix(0, q, q)

  mult <- function(a, b) {
    if (a == b) 0.5 else 1
  }

  for (u in seq_len(q)) {
    i <- params[u, 1]
    j <- params[u, 2]

    for (v in u:q) {
      k <- params[v, 1]
      l <- params[v, 2]

      if (i == j && k == l) {
        Iuv <- n / 2 * imp[i, k]^2
      } else if (i == j && k != l) {
        Iuv <- n * imp[i, k] * imp[i, l]
      } else if (i != j && k == l) {
        Iuv <- n * imp[i, k] * imp[j, k]
      } else {
        Iuv <- n * (
          imp[i, k] * imp[j, l] +
            imp[i, l] * imp[j, k]
        )
      }

      Info[u, v] <- Iuv
      Info[v, u] <- Iuv
    }
  }

  return(Info)
}

#' F-statistics helper function
#' @noRd
#' @keywords internal
phiCompute <- function(est, tolerance = 1e-6, hypothesis) {

  # extract components for easier handling
  n <- est$n
  p <- est$p
  adjacency <- est$adjacency
  emp <- est$emp

  phi <- matrix(0, p, p)

  # loop over nodes
  for (i in seq_len(p)) {
    Ni <- setdiff(which(adjacency[i, ] != 0), i)
    candidates <- setdiff(seq_len(p), c(i, Ni))

    ## RSS for baseline model: X_i ~ X_Ni
    if (length(Ni) == 0) {

      # residual if there are no neighbors (should return 1)
      res0 <- emp[i, i]

    } else {
      # compute residual sum of squares
      empNi <- emp[i, Ni, drop = FALSE]
      empNN <- emp[Ni, Ni, drop = FALSE]

      res0 <- emp[i, i] - empNi %*% solve(empNN, t(empNi))
    }

    for (j in candidates) {

      # grow model by one candidate
      A1 <- c(Ni, j)

      # compute residual sum of squares
      empAi <- emp[i, A1, drop = FALSE]
      empAA <- emp[A1, A1, drop = FALSE]

      res1 <- emp[i, i] - empAi %*% solve(empAA, t(empAi))

      df2 <- n - length(Ni) - 2

      res0 <- as.numeric(res0)
      res1 <- as.numeric(res1)

      phi[i, j] <- (res0 - res1) / (res1 / df2)
    }
  }

  # Add warning if anything is negative
  if (any(phi < 0)) {
    warning('Negative values were computed for some elements of phi while determining the F-statistic. Results may be biased.')
  }

  # return the phi matrix
  return(phi)
}


#' PRC helper function
#' @noRd
#' @keywords internal
prcCompute <- function(est, tolerance = 1e-6, hypothesis) {

  n <- est$n
  p <- est$p

  r <- pmin(pmax(est$empRho, -1 + 1e-12), 1 - 1e-12)

  zval <- abs(atanh(r)) * sqrt(n - p - 1)

  return(zval)

}

#' Create a filter matrix for the relevant edges
#' @noRd
#' @keywords internal
fitFilter <- function(est, tolerance = 1e-6,
  hypothesis = c("local", "embed", "global"),
  triangle = FALSE, constraint = TRUE) {

  hypothesis <- match.arg(hypothesis)

  p <- est$p
  target <- est$target
  adjacency <- est$adjacency

  # tolerance filter for fuzziness
  filt <- abs(est$empRho - est$impRho) > tolerance

  # only restricted / missing edges can contribute
  if (constraint) filt <- filt & adjacency == 0

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


#' Check for all available measures
#' @noRd
#' @keywords internal
implementedFits <- function(prefix = 'fit_', suffix = NULL) {

  ns <- asNamespace('hyponet')

  if (is.null(suffix)) {
    fitFuns <- ls(ns, all.names = TRUE, pattern = paste0('^', prefix, '[a-zA-Z]+'))
  } else {
    fitFuns <- ls(ns, all.names = TRUE, pattern = paste0('^', prefix, '[a-zA-Z]+', '_', suffix, '$'))
  }

  fits <- sub(paste0('^', prefix), '', fitFuns)

  stats::setNames(fitFuns, fits)

}

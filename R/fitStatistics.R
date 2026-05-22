#### Collection of test statistics ----

#' SRMR: Standardized Root Mean Square Residual
#'
#' @keywords internal
fit_srmr <- function(est, tolerance = 1e-6) {

  # get unique entries
  uni <- lower.tri(est$emp, diag = FALSE)

  # determine tolerance filter
  tolFilter <- abs(est$empRho[uni]-est$impRho[uni]) > tolerance

  # determine sum of squares (for differences outside tolerance)
  SS <- sum((est$emp[uni][tolFilter] - est$imp[uni][tolFilter])^2)

  # calculate SRMR only on off-diagonal elements
  srmr <- sqrt(SS / (est$p * (est$p - 1)/2))

  return(srmr)
}


#' (Blockwise) Chi-Square
#'
#' @keywords internal
fit_chisq <- function(est, tolerance = 1e-6) {

  # extract components for easier handling
  n <- est$n
  p <- est$p
  adjacency <- est$adjacency
  emp <- est$emp
  impRho <- est$impRho
  empRho <- est$empRho

  dRho <- empRho - impRho

  # Rerhol
  tolRho <- impRho + sign(dRho) * pmin(abs(dRho), tolerance)

  # Check PD
  if (any(eigen(tolRho, only.values = TRUE)$values < 0)) {
    tolRho <- Matrix::nearPD(tolRho, keepDiag = TRUE)$mat
  }

  diag(tolRho) <- -1
  tol <- solve(-tolRho) |> stats::cov2cor()

  chi <- n * (log(det(tol)) - log(det(emp)) + sum(diag(emp %*% solve(tol))) - p)

  return(chi)
}

#' (Blockwise) RMSEA
#'
#' @keywords internal
fit_rmsea <- function(est, tolerance = 1e-6) {

  chi <- fit_chisq(est, tolerance)

  df <- sum(est$adjacency[lower.tri(est$adjacency)] == 0)

  rmsea <- sqrt(max(((chi/df) - 1) / (est$n - 1) , 0))

  return(rmsea)
}

#' (Blockwise) AIC
#'
#' @keywords internal
# fit_aic <- function(est, tolerance = 1e-6) {
#
#   loglik <- -est$n / 2 * (log(det(est$imp)) + sum(diag(est$emp %*% solve(est$imp))))
#
#   k <- sum(est$adjacency[lower.tri(est$adjacency)] != 0)
#
#   aic <- -2 * loglik + 2 * k
#
#   return(aic)
# }

#' (Blockwise) BIC
#'
#' @keywords internal
# fit_bic <- function(est, tolerance = 1e-6) {
#
#   loglik <- -est$n / 2 * (log(det(est$imp)) + sum(diag(est$emp %*% solve(est$imp))))
#
#   k <- sum(est$adjacency[lower.tri(est$adjacency)] != 0)
#
#   bic <- -2 * loglik + log(est$n) * k
#
#   return(bic)
# }

#' F-statistics helper function
#'
#' @keywords internal
phiCompute <- function(est, tolerance = 1e-6) {

  # extract components for easier handling
  n <- est$n
  p <- est$p
  adjacency <- est$adjacency
  emp <- est$emp
  imp <- est$imp

  phi <- matrix(0, p, p)

  # determine inclusion via tolerance
  include <- abs(est$empRho - est$impRho) > tolerance


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

      # skip to next iteration if below tolerance
      if (!include[i, j]) next

      # grow model by one candidate
      A1 <- c(Ni, j)

      # compute residual sum of squares
      empAi <- emp[i, A1, drop = FALSE]
      empAA <- emp[A1, A1, drop = FALSE]

      res1 <- emp[i, i] - empAi %*% solve(empAA, t(empAi))

      df2 <- n - length(Ni) - 2

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

#' F-sum statistic
#'
#' @keywords internal
fit_fsum <- function(est, tolerance = 1e-6) {

  phi <- phiCompute(est = est, tolerance = tolerance)

  fsum <- sum(phi)

  return(fsum)
}

#' F-max statistic
#'
#' @keywords internal
fit_fmax <- function(est, tolerance = 1e-6) {

  phi <- phiCompute(est = est, tolerance = tolerance)

  fmax <- max(phi)

  return(fmax)
}

#' PRC helper function
#'
#' @keywords internal
prcCompute <- function(est, tolerance = 1e-6) {

  n <- est$n
  p <- est$p

  r <- pmin(pmax(est$empRho, -1 + 1e-12), 1 - 1e-12)

  zval <- abs(atanh(r)) * sqrt(n - p - 1)

  return(zval)

}

#' PRC sum of squares
#'
#' @keywords internal
fit_prcss <- function(est, tolerance = 1e-6) {

  zval <- prcCompute(est = est, tolerance = tolerance)

  include <- abs(est$empRho - est$impRho) > tolerance

  filt <- lower.tri(zval) &
    include &
    est$adjacency == 0 # ensures only discrepancies on restricted edges are used

  prcss <- sum(zval[filt]^2)

  return(prcss)
}

fit_prcsa <- function(est, tolerance = 1e-6) {

  zval <- prcCompute(est = est, tolerance = tolerance)

  include <- abs(est$empRho - est$impRho) > tolerance

  filt <- lower.tri(zval) &
    include &
    est$adjacency == 0 # ensures only discrepancies on restricted edges are used

  prcsa <- sum(abs(zval[filt]))

  return(prcsa)
}


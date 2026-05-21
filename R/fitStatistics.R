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

  chi <- est$n * (log(det(est$imp)) - log(det(est$emp)) + sum(diag(est$emp %*% solve(est$imp))) - est$p)

  return(chi)
}

#' (Blockwise) RMSEA
#'
#' @keywords internal
fit_rmsea <- function(est, tolerance = 1e-6) {

  chi <- fit_chisq(est, tolerance)

  df <- sum(adjacency[lower.tri(adjacency)] == 0)

  rmsea <- sqrt(max(((chi/df) - 1) / (est$n - 1) , 0))

  return(rmsea)
}

#' (Blockwise) AIC
#'
#' @keywords internal
fit_aic <- function(est, tolerance = 1e-6) {

  loglik <- -est$n / 2 * (log(det(est$imp)) + sum(diag(est$emp %*% solve(est$imp))))

  k <- sum(est$adjacency[lower.tri(est$adjacency)] != 0)

  aic <- -2 * loglik + 2 * k

  return(aic)
}

#' (Blockwise) BIC
#'
#' @keywords internal
fit_bic <- function(est, tolerance = 1e-6) {

  loglik <- -est$n / 2 * (log(det(est$imp)) + sum(diag(est$emp %*% solve(est$imp))))

  k <- sum(est$adjacency[lower.tri(est$adjacency)] != 0)

  bic <- -2 * loglik + log(est$n) * k

  return(bic)
}

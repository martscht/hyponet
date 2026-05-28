#### Collection of test statistics ----

### Constraint statistics ----
# these test only the 0-contraints implied by the model

#' SRMR: Standardized Root Mean Square Residual
#' @noRd
#' @keywords internal
fit_srmr_con <- function(est, tolerance = 1e-6, hypothesis) {

  # filter for elements to include
  filt <- fitFilter(est = est, tolerance = tolerance,
    hypothesis = hypothesis, triangle = TRUE, constraint = TRUE)

  if (!any(filt)) return(0)

  # number of elements included
  p <- sum(filt)

  # determine sum of squares (for differences outside tolerance)
  SS <- sum((est$empRho[filt] - est$impRho[filt])^2)

  # calculate SRMR only on off-diagonal elements
  srmr <- sqrt(SS / p)

  return(srmr)
}


#' Global modification indices
#' @noRd
#' @keywords internal
fit_modindex_con <- function(est, tolerance = 1e-6, hypothesis) {

  n <- est$n
  p <- est$p
  emp <- est$emp
  imp <- est$imp
  adjacency <- est$adjacency

  # currently free nuisance parameters:
  # all diagonal precision entries + all allowed undirected edges
  diag_params <- cbind(seq_len(p), seq_len(p))
  edge_params <- which(adjacency == 1 & upper.tri(adjacency), arr.ind = TRUE)

  free_params <- rbind(diag_params, edge_params)

  # candidate constrained edges from shared fit filter
  filt <- fitFilter(est = est, tolerance = tolerance,
    hypothesis = hypothesis, triangle = FALSE, constraint = TRUE)

  # MI is undirected: include edge if either direction is selected
  filt <- filt | t(filt)

  # keep each undirected edge once
  cand_params <- which(filt & upper.tri(filt), arr.ind = TRUE)

  # abort if nothing is being tested
  if (nrow(cand_params) == 0) return(0)

  all_params <- rbind(free_params, cand_params)

  n_free <- nrow(free_params)
  n_cand <- nrow(cand_params)

  Info <- infoMatrix(imp = imp, params = all_params, n = n)

  # scores for candidate precision entries
  score_cand <- n * (imp[cand_params] - emp[cand_params])

  idx_free <- seq_len(n_free)
  idx_cand <- n_free + seq_len(n_cand)

  I_ff <- Info[idx_free, idx_free, drop = FALSE]
  I_fc <- Info[idx_free, idx_cand, drop = FALSE]
  I_cf <- Info[idx_cand, idx_free, drop = FALSE]
  I_cc <- Info[idx_cand, idx_cand, drop = FALSE]

  # efficient information for candidates after adjusting for nuisance parameters
  I_eff <- I_cc - I_cf %*% solve(I_ff, I_fc)

  MI <- as.numeric(t(score_cand) %*% solve(I_eff, score_cand))

  return(MI)
}


#' Modindex-based RMSEA
#' @noRd
#' @keywords internal
fit_rmsea_con <- function(est, tolerance = 1e-6, hypothesis) {

  chi <- fit_modindex_con(est, tolerance, hypothesis)

  # candidate constrained edges from shared fit filter
  filt <- fitFilter(est = est, tolerance = -Inf,
    hypothesis = hypothesis, triangle = TRUE, constraint = TRUE)

  df <- sum(filt)

  rmsea <- sqrt(max(((chi/df) - 1) / (est$n - 1) , 0))

  return(rmsea)
}

#' F-sum statistic
#' @noRd
#' @keywords internal
fit_fsum_con <- function(est, tolerance = 1e-6, hypothesis) {

  phi <- phiCompute(est = est, tolerance = tolerance)

  # filter for elements to include
  filt <- fitFilter(est = est, tolerance = tolerance,
    hypothesis = hypothesis, triangle = FALSE, constraint = TRUE)

  if (!any(filt)) return(0)

  fsum <- sum(phi[filt])

  return(fsum)
}

#' F-max statistic
#' @noRd
#' @keywords internal
fit_fmax_con <- function(est, tolerance = 1e-6, hypothesis) {

  phi <- phiCompute(est = est, tolerance = tolerance)

  # filter for elements to include
  filt <- fitFilter(est = est, tolerance = tolerance,
    hypothesis = hypothesis, triangle = FALSE, constraint = TRUE)

  if (!any(filt)) return(0)

  fmax <- max(phi[filt])

  return(fmax)
}

#' PRC sum of squares
#' @noRd
#' @keywords internal
fit_prcss_con <- function(est, tolerance = 1e-6, hypothesis) {

  zval <- prcCompute(est = est, tolerance = tolerance)

  # filter for elements to include
  filt <- fitFilter(est = est, tolerance = tolerance,
    hypothesis = hypothesis, triangle = TRUE, constraint = TRUE)

  if (!any(filt)) return(0)

  prcss <- sum(zval[filt]^2)

  return(prcss)
}

#' PRC sum of absolutes
#' @noRd
#' @keywords internal
fit_prcsa_con <- function(est, tolerance = 1e-6, hypothesis) {

  zval <- prcCompute(est = est, tolerance = tolerance)

  # filter for elements to include
  filt <- fitFilter(est = est, tolerance = tolerance,
    hypothesis = hypothesis, triangle = TRUE, constraint = TRUE)

  if (!any(filt)) return(0)

  prcsa <- sum(abs(zval[filt]))

  return(prcsa)
}

### Approximation statistics ----
# these test all elements to check how close the implied network
# approximates the empirical network

#' (Blockwise) Chi-Square
#' @noRd
#' @keywords internal
fit_chisq_app <- function(est, tolerance = 1e-6, hypothesis) {

  # filter for elements to include
  if (hypothesis == 'local') {
    filt <- fitFilter(est = est, tolerance = tolerance,
      hypothesis = hypothesis, triangle = FALSE, constraint = FALSE)
  } else {
    # warning('Chi-Square cannot be computed for embedding hypotheses. Values reported relate to the entire network.')
    filt <- fitFilter(est = est, tolerance = tolerance,
      hypothesis = 'global', triangle = FALSE, constraint = FALSE)
  }

  if (!any(filt)) return(0)

  # extract components for easier handling
  n <- est$n
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
  tol <- solve(-tolRho)

  tol <- tol[apply(filt, 1, any), apply(filt, 2, any)]
  emp <- emp[apply(filt, 1, any), apply(filt, 2, any)]

  p <- ncol(tol)

  chi <- n * (log(det(tol)) - log(det(emp)) + sum(diag(emp %*% solve(tol))) - p)

  return(chi)
}


#' (Blockwise) RMSEA
#' @noRd
#' @keywords internal
fit_rmsea_app <- function(est, tolerance = 1e-6, hypothesis) {

  # if (hypothesis == 'embed') {
  #   warning('RMSEA cannot be computed for embedding hypotheses. Values reported relate to the entire network.')
  # }

  hypothesis <- ifelse(hypothesis == 'embed', 'global', hypothesis)

  chi <- suppressWarnings(fit_chisq_app(est, tolerance, hypothesis))

  # candidate constrained edges from shared fit filter
  filt <- fitFilter(est = est, tolerance = tolerance,
    hypothesis = hypothesis, triangle = TRUE)

  if (!any(filt)) return(0)

  df <- sum(filt)

  rmsea <- sqrt(max(((chi/df) - 1) / (est$n - 1) , 0))

  return(rmsea)
}

#' SRMR: Standardized Root Mean Square Residual
#' @noRd
#' @keywords internal
fit_srmr_app <- function(est, tolerance = 1e-6, hypothesis) {

  # filter for elements to include
  filt <- fitFilter(est = est, tolerance = tolerance,
    hypothesis = hypothesis, triangle = TRUE, constraint = FALSE)

  if (!any(filt)) return(0)

  # number of elements included
  p <- sum(filt)

  # determine sum of squares (for differences outside tolerance)
  SS <- sum((est$emp[filt] - est$imp[filt])^2)

  # calculate SRMR only on off-diagonal elements
  srmr <- sqrt(SS / p)

  return(srmr)
}


#' PRC sum of squares
#' @noRd
#' @keywords internal
fit_prcss_app <- function(est, tolerance = 1e-6, hypothesis) {

  zval <- prcCompute(est = est, tolerance = tolerance)

  # filter for elements to include
  filt <- fitFilter(est = est, tolerance = tolerance,
    hypothesis = hypothesis, triangle = TRUE, constraint = FALSE)

  if (!any(filt)) return(0)

  prcss <- sum(zval[filt]^2)

  return(prcss)
}

#' PRC sum of absolutes
#' @noRd
#' @keywords internal
fit_prcsa_app <- function(est, tolerance = 1e-6, hypothesis) {

  zval <- prcCompute(est = est, tolerance = tolerance)

  # filter for elements to include
  filt <- fitFilter(est = est, tolerance = tolerance,
    hypothesis = hypothesis, triangle = TRUE, constraint = FALSE)

  if (!any(filt)) return(0)

  prcsa <- sum(abs(zval[filt]))

  return(prcsa)
}

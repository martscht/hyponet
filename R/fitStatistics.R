#### Collection of test statistics ----

#' SRMR: Standardized Root Mean Square Residual
#' @noRd
#' @keywords internal
fit_srmr <- function(est, tolerance = 1e-6, hypothesis) {

  # filter for elements to include
  filt <- fitFilter(est = est, tolerance = tolerance,
    hypothesis = hypothesis, triangle = TRUE)

  if (!any(filt)) return(NA)

  # number of elements included
  p <- sum(filt)

  # determine sum of squares (for differences outside tolerance)
  SS <- sum((est$emp[filt] - est$imp[filt])^2)

  # calculate SRMR only on off-diagonal elements
  srmr <- sqrt(SS / p)

  return(srmr)
}

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


#' Global modification indices
#' @noRd
#' @keywords internal
fit_modindex <- function(est, tolerance = 1e-6, hypothesis) {

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
  filt <- fitFilter(est = est, tolerance = tolerance, hypothesis = hypothesis, triangle = FALSE)

  # MI is undirected: include edge if either direction is selected
  filt <- filt | t(filt)

  # keep each undirected edge once
  cand_params <- which(filt & upper.tri(filt), arr.ind = TRUE)

  # abort if nothing is being tested
  if (nrow(cand_params) == 0) return(NA)

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


#' (Blockwise) Chi-Square
#' @noRd
#' @keywords internal
fit_chisq <- function(est, tolerance = 1e-6, hypothesis) {

  # filter for elements to include
  # filt <- fitFilter(est = est, tolerance = tolerance,
  #   hypothesis = hypothesis, triangle = TRUE)

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

  # if (hypothesis == 'local') {
  #   tol <- tol[filt]
  #   emp <- emp[filt]
  # } else {
  #   warning('Chi-Square cannot be computed for embedding hypotheses. Values reported relate to the entire network.')
  # }

  chi <- n * (log(det(tol)) - log(det(emp)) + sum(diag(emp %*% solve(tol))) - p)

  return(chi)
}

#' (Blockwise) RMSEA
#' @noRd
#' @keywords internal
fit_rmsea <- function(est, tolerance = 1e-6, hypothesis) {

  chi <- fit_modindex(est, tolerance, hypothesis)

  # candidate constrained edges from shared fit filter
  filt <- fitFilter(est = est, tolerance = tolerance, hypothesis = hypothesis, triangle = FALSE)
  df <- sum(filt)

  rmsea <- sqrt(max(((chi/df) - 1) / (est$n - 1) , 0))

  return(rmsea)
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

#' F-sum statistic
#' @noRd
#' @keywords internal
fit_fsum <- function(est, tolerance = 1e-6, hypothesis) {

  phi <- phiCompute(est = est, tolerance = tolerance)

  # filter for elements to include
  filt <- fitFilter(est = est, tolerance = tolerance,
    hypothesis = hypothesis, triangle = FALSE)

  fsum <- sum(phi[filt])

  return(fsum)
}

#' F-max statistic
#' @noRd
#' @keywords internal
fit_fmax <- function(est, tolerance = 1e-6, hypothesis) {

  phi <- phiCompute(est = est, tolerance = tolerance)

  # filter for elements to include
  filt <- fitFilter(est = est, tolerance = tolerance,
    hypothesis = hypothesis, triangle = FALSE)

  fmax <- max(phi[filt])

  return(fmax)
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

#' PRC sum of squares
#' @noRd
#' @keywords internal
fit_prcss <- function(est, tolerance = 1e-6, hypothesis) {

  zval <- prcCompute(est = est, tolerance = tolerance)

  # filter for elements to include
  filt <- fitFilter(est = est, tolerance = tolerance,
    hypothesis = hypothesis, triangle = TRUE)

  prcss <- sum(zval[filt]^2)

  return(prcss)
}

#' PRC sum of absolutes
#' @noRd
#' @keywords internal
fit_prcsa <- function(est, tolerance = 1e-6, hypothesis) {

  zval <- prcCompute(est = est, tolerance = tolerance)

  # filter for elements to include
  filt <- fitFilter(est = est, tolerance = tolerance,
    hypothesis = hypothesis, triangle = TRUE)

  prcsa <- sum(abs(zval[filt]))

  return(prcsa)
}


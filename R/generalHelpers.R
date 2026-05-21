#### Collection of general helper functions ----

# Convert correlation matrix to partial correlations
c2pc <- function(R) {

  Omega <- solve(R)
  Rho <- -stats::cov2cor(Omega)
  diag(Rho) <- 1

  return(Rho)
}

impliedCor <- function(empirical, adjacency, n, ...) {

  # Rephrase constraints
  constr <- which(adjacency == 0, arr.ind = TRUE)
  constr <- constr[constr[, 1] != constr[, 2], ]

  # Estimate via glasso
  Net <- glasso::glasso(s = empirical, rho = 0, nobs = n,
    zero = constr) |> suppressWarnings() |> suppressMessages()
  SigmaHat <- stats::cov2cor(Net$w)

  return(SigmaHat)
}

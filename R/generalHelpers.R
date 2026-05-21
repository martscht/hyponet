#### Collection of general helper functions ----

# Convert correlation matrix to partial correlations
c2pc <- function(R) {

  Omega <- solve(R)
  Rho <- -stats::cov2cor(Omega)
  diag(Rho) <- 1

  return(Rho)
}

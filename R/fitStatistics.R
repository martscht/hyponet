#### Collection of test statistics ----

SRMR <- function(data, adjacency, target = 1:ncol(data),
  tolerance = 1e-6) {

  # determine sample size
  n <- nrow(data)

  # determine p
  p <- length(target)

  # get empirical and implied correlation matrices
  emp <- stats::cor(data)

  constr <- which(adjacency == 0, arr.ind = TRUE)
  constr <- constr[constr[, 1] != constr[, 2], ]
  Net <- glasso::glasso(s = emp, rho = 0, nobs = n,
    zero = constr) |> suppressWarnings() |> suppressMessages()
  imp <- stats::cov2cor(Net$w)

  # filter down to local version
  emp <- emp[target, target]
  imp <- imp[target, target]

  # same thing for partial correlation matrices
  empRho <- c2pc(emp)
  impRho <- -stats::cov2cor(Net$wi)
  diag(impRho) <- 1

  # get unique entries
  uni <- lower.tri(emp, diag = FALSE)

  # determine tolerance filter
  tolFilter <- abs(empRho[uni]-impRho[uni]) > tolerance

  # determine sum of squares (for differences outside tolerance)
  SS <- sum((emp[uni][tolFilter] - imp[uni][tolFilter])^2)

  # calculate SRMR only on off-diagonal elements
  srmr <- sqrt(SS / (p * (p - 1)/2))

  return(srmr)
}

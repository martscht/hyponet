#### Collection of test statistics ----

# SRMR: Standardized Root Mean Square Residual
SRMR <- function(data, adjacency, target = 1:ncol(data),
  tolerance = 1e-6) {

  # determine sample size
  n <- nrow(data)

  # determine p
  p <- length(target)

  # get empirical and implied correlation matrices
  emp <- stats::cor(data)
  imp <- impliedCor(emp, adjacency, n)

  # same thing for partial correlation matrices
  empRho <- c2pc(emp)
  impRho <- c2pc(imp)

  # filter down to local version
  emp <- emp[target, target]
  imp <- imp[target, target]
  empRho <- empRho[target, target]
  impRho <- impRho[target, target]

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


# (Blockwise) Chi²
ChiSq <- function(data, adjacency, target = 1:ncol(data),
  tolerance = 1e-6) {

  # determine sample size
  n <- nrow(data)
  # determine p
  p <- length(target)

  # get empirical and implied correlation matrices
  emp <- stats::cor(data)
  imp <- impliedCor(emp, adjacency, n)

  # same thing for partial correlation matrices
  empRho <- c2pc(emp)
  impRho <- c2pc(imp)

  # filter down to local version
  emp <- emp[target, target]
  imp <- imp[target, target]
  empRho <- empRho[target, target]
  impRho <- impRho[target, target]

  chi <- n * (log(det(imp)) - log(det(emp)) + sum(diag(emp %*% solve(imp))) - p)

  return(chi)
}

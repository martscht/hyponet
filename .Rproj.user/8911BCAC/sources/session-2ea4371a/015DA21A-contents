#### Helper functions for exchangeable sampling ----

# Orthogonal Vector for a matrix (algorithm 1)
orthVec <- function(X) {

  # QR decomposition
  QR <- qr(cbind(X, rnorm(nrow(X))))

  # Generate vector
  vec <- qr.Q(QR)[, ncol(X) + 1]
  vec <- vec * sample(c(-1, 1), 1)

  # Output
  return(vec)
}

# Generate x tilde for a single variable
algo1 <- function(data, i, adjacency) {

  # Determine neighbors of i
  N <- which(adjacency[i, ] != 0) |> setdiff(i)

  # Isolate relevant variable
  Xi <- data[, i]

  # Base rotation on sample mean (no neighbors)
  if (length(N) == 0) {
    Xi <- mean(Xi) + orthVec(matrix(1, nrow = length(Xi))) * sqrt(sum((Xi - mean(Xi))^2))
  } else if (length(Xi) >= 2 + length(N)) {
    tmp <- lm(Xi ~ data[, N])
    Xi <- predict(tmp) + sqrt(sum(resid(tmp)^2)) * orthVec(cbind(1, data[, N]))
  }

  # Output new variable
  return(Xi)

}

# Generate X tilde for entire data frame
algo2step1 <- function(data, target, adjacency) {

  Xtilde <- data

  # Replace each x with x tilde in sequence
  for (i in target) {
    Xtilde[, i] <- algo1(Xtilde, i, adjacency)
  }

  return(Xtilde)
}


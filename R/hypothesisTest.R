hypothesisTest <- function(data, adjacency, target = NULL,
  fits = 'all', tolerance = 1e-6,
  resampleType = c('exchangeableSampling', 'parametricBootstrap'),
  ...) {

  # If resampling was not done ahead of time, do it now
  if (!inherits(data, 'hyponetData')) {
    if (is.null(target)) target <- 1:ncol(data)

    resampleType <- match.arg(resampleType)

    data <- switch(resampleType,
      exchangeableSampling = exchangeableSampling(
        data = data,
        adjacency = adjacency,
        target = target, ...),
      parametricBootstrap = parametricBootstrap(
        data = data,
        adjacency = adjacency,
        ...))
  }

  if (is.null(target)) target <- 1:ncol(data$original)

  # determine fit statistics for original data (standardized)
  originalFit <- determineFits(data = data$original, adjacency = adjacency,
    target = target, fits = fits, tolerance = tolerance)

  # determine fit statistics for copies
  resampleFit <- determineFits(data = data$copies, adjacency = adjacency,
    target = target, fits = fits, tolerance = tolerance)

  # determine p-value
  larger <- apply(resampleFit, 1, function(x) x > originalFit) |> rowSums()
  smaller <- apply(resampleFit, 1, function(x) x < originalFit) |> rowSums()
  equal <- apply(resampleFit, 1, function(x) x == originalFit) |> rowSums()

  pvalues <- (larger + .5*equal) / (larger + smaller + equal)

  out <- rbind(originalFit, pvalues)
  rownames(out) <- c('Statistic', 'p-Value')

  # class(out) <- 'hyponetResult'

  return(out)
}

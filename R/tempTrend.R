#' Long-term local climatic trends
#'
#' Function to calculate temporal trend from a raster series
#' of a climatic variable. This trend is to be used for the calculation of the
#' gradient-based climate velocity using gVoCC.
#'
#' @param r \code{RasterStack} containing a time series of (annual, seasonal, monthly...) values of
#' the climatic variable for the period of interest.
#' @param th \code{Integer} minimum number of observations in the series needed to
#' calculate the trend at each cell.
#'
#' @return A \code{RasterStack} containing the cell-specific temporal trends
#' extracted from simple linear regressions of the climatic variable against time
#' ("slpTrends" in degree Celsius per year), together with their standard
#' errors ("seTrends") and statistical significance ("sigTrends").
#'
#' @seealso{\code{\link{spatGrad}}, \code{\link{gVoCC}}}
#'
#' @export
#' @examples
#'
#' HSST <- terra::rast(system.file("extdata", "HadiSST.tif", package = "VoCCdata"))
#'
#' yrSST <- sumSeries(HSST,
#'   p = "1969-01/2009-12", yr0 = "1955-01-01", l = terra::nlyr(HSST),
#'   fun = function(x) colMeans(x, na.rm = TRUE), freqin = "months", freqout = "years"
#' )
#'
#' # Mean annual SST trend (minimum threshold of 10 years of data), with SE and p-values.
#'
#' tr <- tempTrend(yrSST, th = 10)
#'
#' terra::plot(tr)
#'
tempTrend <- function(r, th) {
  y <- terra::values(r)
  ocean <- which(rowSums(is.na(y)) != ncol(y)) # remove land cells
  y <- t(y[ocean, ])

  # OPTIMIZED: Vectorized NA counting using colSums instead of apply
  N <- colSums(!is.na(y))
  ind <- which(N >= th)
  y <- y[, ind, drop = FALSE] # drop cells with less than th observations
  N <- N[ind] # Use pre-calculated N values

  # OPTIMIZED: More efficient x matrix creation using outer product
  n_layers <- terra::nlyr(r)
  n_cells <- ncol(y)

  # Create time vector once and replicate efficiently
  time_vec <- seq_len(n_layers)
  x <- matrix(rep(time_vec, n_cells), nrow = n_layers, ncol = n_cells)

  # put NA values into the x values so they correspond with y
  x[is.na(y)] <- NA

  # OPTIMIZED: Vectorized sum calculations using colSums instead of apply
  sx <- colSums(x, na.rm = TRUE)
  sy <- colSums(y, na.rm = TRUE)
  sxx <- colSums(x^2, na.rm = TRUE)
  syy <- colSums(y^2, na.rm = TRUE)
  sxy <- colSums(x * y, na.rm = TRUE)

  # OPTIMIZED: Vectorized slope calculations with better numerical stability
  denominator <- N * sxx - sx^2
  # Add small epsilon to prevent division by zero
  denominator[denominator == 0] <- .Machine$double.eps

  slope <- (N * sxy - (sx * sy)) / denominator
  sres <- (N * syy - sy^2 - slope^2 * denominator) / (N * (N - 2))
  SE <- suppressWarnings(sqrt((N * sres) / denominator))
  Test <- slope / SE

  # OPTIMIZED: Vectorized p-value calculation with bounds checking
  df_vals <- pmax(N - 2, 1)  # Ensure df >= 1
  p <- 2 * stats::pt(abs(Test), df = df_vals, lower.tail = FALSE)

  # OPTIMIZED: Direct raster creation and assignment
  template_raster <- terra::rast(r[[1]])
  slpTrends <- seTrends <- sigTrends <- template_raster

  # Use vectorized assignment
  ocean_ind <- ocean[ind]
  slpTrends[ocean_ind] <- slope
  seTrends[ocean_ind] <- SE
  sigTrends[ocean_ind] <- p

  output <- c(slpTrends, seTrends, sigTrends)
  names(output) <- c("slpTrends", "seTrends", "sigTrends")

  return(output)
}

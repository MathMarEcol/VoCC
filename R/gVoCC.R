#' Gradient-based climate velocity
#'
#' Function to calculate the velocity of climate change after Burrows et al. (2011)
#' based on local climatic temporal trends and spatial gradients.
#'
#' @param tempTrend The output from the tempTrend function containing the long-term linear climatic trends.
#' @param spatGrad The output from the spatGrad function containing the magnitudes and angles for the spatial climatic gradient.
#'
#' @return A \code{RasterStack} containing the climate velocity magnitude ("voccMag",
#' km/yr for unprojected rasters and spatial unit/year for projected rasters) and angle("voccAng" in
#' degrees north: 0N, 90E, 180S and 270W).
#'
#' @references \href{http://science.sciencemag.org/content/334/6056/652}{Burrows et al. 2011}. The pace of shifting climate
#' in marine and terrestrial ecosystems. Science, 334, 652-655.
#'
#' @seealso{\code{\link{tempTrend}}, \code{\link{spatGrad}}}
#'
#' @export
#'
#' @examples
#'
#' HSST <- terra::rast(system.file("extdata", "HadiSST.tif", package = "VoCCdata"))
#'
#' yrSST <- sumSeries(HSST,
#'   p = "1960-01/2009-12", yr0 = "1955-01-01", l = terra::nlyr(HSST),
#'   fun = function(x) colMeans(x, na.rm = TRUE), freqin = "months", freqout = "years"
#' )
#' tr <- tempTrend(yrSST, th = 10)
#' sg <- spatGrad(yrSST, th = 0.0001, projected = FALSE)
#'
#' # Magnitude and angle of the climate velocity (km/yr) 1960-2009
#'
#' v <- gVoCC(tr, sg)
#' terra::plot(v)
#' }
#'
gVoCC <- function(tempTrend, spatGrad) {
  VoCC <- tempTrend[[1]] / spatGrad[[1]]

  # OPTIMIZATION: Extract values once and use vectorized operations
  VoCC_values <- terra::values(VoCC)
  spatGrad_ang_values <- terra::values(spatGrad[[2]])

  # Vectorized angle calculation with proper NA handling
  warming_cells <- VoCC_values > 0 & !is.na(VoCC_values)
  VoCCang_values <- spatGrad_ang_values
  VoCCang_values[warming_cells] <- spatGrad_ang_values[warming_cells] + 180

  # Handle angle wrapping with NA protection
  needs_wrapping <- !is.na(VoCCang_values) & VoCCang_values >= 360
  VoCCang_values[needs_wrapping] <- VoCCang_values[needs_wrapping] - 360

  # Create output raster efficiently
  VoCCang <- spatGrad[[2]]
  terra::values(VoCCang) <- VoCCang_values

  output <- c(VoCC, VoCCang)
  names(output) <- c("voccMag", "voccAng")
  return(output)
}

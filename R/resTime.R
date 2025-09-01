#' Climatic residence time of a polygon
#'
#' Function to calculate VoCC-based residence time of isotherms within a polygon after Loaire et al. (2009)
#'
#' @param pg \code{sf} object or \code{terra::vect} object containing the polygons for which
#' the residence time is to be calculated. The polygons must be on the same coordinate system as vel.
#' @param vel \code{raster} with climate velocity (km/year) for the period of interest.
#' @param areapg \code{vector} with the area (in km2) of the polygons. Use NA (default) to calculate internally if field not avilable.
#'
#' @return a \code{data.frame} containing for each polygon its ID, mean velocity (km/yr),
#' diameter of the equivalent circle (km), and residence time (years) as the ratio D/vel.
#'
#' @references \href{https://www.nature.com/articles/nature08649}{Loarie et al. 2009}. The velocity of climate change. Nature, 462, 1052-1055.
#'
#' @seealso{\code{\link{gVoCC}}}
#'
#' @export
#' @author Jorge Garcia Molinos
#' @examples
#'
#' # Load example Exclusive Economic Zone polygon
#' \dontrun{
#' EEZ <- VoCC_get_data("EEZ.gpkg")
#' HSST <- VoCC_get_data("HSST.tif")
#'
#' yrSST <- sumSeries(HSST,
#'   p = "1969-01/2009-12", yr0 = "1955-01-01", l = terra::nlyr(HSST),
#'   fun = function(x) colMeans(x, na.rm = TRUE),
#'   freqin = "months", freqout = "years"
#' )
#' tr <- tempTrend(yrSST, th = 10)
#' sg <- spatGrad(yrSST, th = 0.0001, projected = FALSE)
#' v <- gVoCC(tr, sg)
#' vel <- v[[1]]
#'
#' # Calculating area internally
#' a1 <- resTime(EEZ, vel, areapg = NA)
#' a1
#'
#' # Using the area field from the polygon data table
#' a2 <- resTime(EEZ, vel, areapg = as.numeric(as.numeric(levels(EEZ$Area_km2))[EEZ$Area_km2]))
#' a2
#'
#' # Using a user defined polygon
#' x_coord <- c(-28, -20, -20.3, -25.5)
#' y_coord <- c(60, 61, 63, 62)
#' coords <- matrix(c(x_coord, y_coord), ncol = 2)
#' poly_sf <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(coords))))
#' a3 <- resTime(poly_sf, vel, areapg = NA)
#'
#' terra::plot(vel)
#' plot(sf::st_geometry(EEZ), add = TRUE)
#' plot(sf::st_geometry(poly_sf), add = TRUE)
#' }
resTime <- function(pg, vel, areapg = NA) {

  resTim <- v <- d <- NULL # Fix devtools check warnings

  # Handle both sf and terra::vect objects
  if (inherits(pg, "sf")) {
    pg_vect <- terra::vect(pg)
    n_features <- nrow(pg)
  } else if (inherits(pg, "SpatVector")) {
    pg_vect <- pg
    n_features <- nrow(pg)
  } else {
    stop("pg must be an sf object or a terra::vect (SpatVector) object")
  }

  # Create ID sequence based on number of features
  RT <- data.table::data.table(ID = seq_len(n_features))

  # Extract velocity values using terra::vect for efficiency
  extracted_values <- terra::extract(vel, pg_vect, fun = mean, na.rm = TRUE)
  
  # Handle case where extraction returns NULL or empty results
  if (is.null(extracted_values) || nrow(extracted_values) == 0) {
    stop("No values could be extracted from the raster. Check that polygons overlap with the raster and have the same coordinate system.")
  }
  
  # Extract the mean values, handling potential column name variations
  if ("mean" %in% names(extracted_values)) {
    RT[, v := extracted_values$mean]
  } else if (ncol(extracted_values) >= 2) {
    # If no 'mean' column, use the second column (first is usually ID)
    RT[, v := extracted_values[[2]]]
  } else {
    stop("Unexpected format from terra::extract(). Please check your input data.")
  }

  # If not provided, calculate the area of the polygon
  if (all(is.na(areapg))) {
    if (inherits(pg, "sf")) {
      # Calculate area in km2 using sf::st_area
      areapg <- as.numeric(sf::st_area(pg)) / 1000000 # Convert from m2 to km2
    } else {
      # Calculate area in km2 using terra::expanse
      areapg <- terra::expanse(pg_vect, unit = "km")
    }
  }

  RT[, d := 2 * sqrt((areapg / pi))]
  RT[, resTim := abs(d / v)]
  return(RT[])
}

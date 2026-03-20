#' Internal. Angle associated to the spatial gradient
#' @param dx \code{numeric} giving the longitudinal gradient component
#' @param dy \code{numeric} giving the latitudinal gradient component
#' angulo()

angulo <- function(dx, dy) {
  # OPTIMIZATION: Fully vectorized angle calculation - eliminates apply() loop
  # Convert atan to degrees once for all values
  atan_deg <- CircStats::deg(atan(dx / dy))

  # Vectorized conditional logic - much faster than apply()
  angle <- ifelse(dy < 0, 180 + atan_deg,
                 ifelse(dx < 0, 360 + atan_deg, atan_deg))

  return(angle)
}

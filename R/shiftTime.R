#' Shift in timing of seasonal climatology
#'
#' Function to calculate the seasonal shift in the arriving of typical seasonal climates
#' for a given period of interest as per Burrows et al. (2011).
#'
#' @param r \code{stack} with monthly values of the climatic
#' variable for the period of interest.
#' @param yr0 \code{integer} specifying the first year in the series.
#' @param yr1 \code{integer} specifying the initial year for the period of interest.
#' @param yr2 \code{integer} specifying the end year for the period of interest.
#' @param th \code{integer} minimum number of non NAs in the series needed to
#' calculate the trend (default 3).
#' @param m \code{integer} number (1-12) of the month for which the shift is to be calculated
#'
#' @return a \code{stack} with the long-term monthly trend (C/year for temperature in degrees; "mTrend"),
#' seasonal rate of change (C/month; "seaRate"), and seasonal shift (day/decade; "seaShift").
#'
#' @references \href{http://science.sciencemag.org/content/334/6056/652}{Burrows et al. 2011}. The pace of
#' shifting climate in marine and terrestrial ecosystems. Science, 334, 652-655.
#'
#' @export
#' @examples
#'
#' HSST <- terra::rast(system.file("extdata", "HadiSST.tif", package = "VoCCdata"))
#' Apr <- shiftTime(HSST, yr1 = 1960, yr2 = 2009, yr0 = 1955, th = 10, m = 4)
#'
#' terra::plot(Apr)
#'
shiftTime <- function(r, yr1, yr2, yr0, th, m) {
  # OPTIMIZATION: Pre-calculate constants and month indices
  yr_offset <- (yr1 - yr0) * 12

  # 1. Long term trends in monthly values (e.g. deg/year if temperature)
  m1 <- yr_offset + m
  m2 <- ((yr2 - yr0) * 12) + m
  r1 <- r[[seq(m1, m2, by = 12)]]
  trend <- tempTrend(r1, th)[[1]]

  # OPTIMIZATION: Pre-calculate preceding and following months with bounds checking
  prev_month <- if (m == 1) 12 else (m - 1)
  next_month <- if (m == 12) 1 else (m + 1)

  # 2. seasonal rate of shift in temperature centered on each month (deg/month)
  # preceding month
  m1_prev <- yr_offset + prev_month
  m2_prev <- ((yr2 - yr0) * 12) + prev_month
  x2 <- r[[seq(m1_prev, m2_prev, by = 12)]]

  # following month
  m1_next <- yr_offset + next_month
  m2_next <- ((yr2 - yr0) * 12) + next_month
  x3 <- r[[seq(m1_next, m2_next, by = 12)]]

  # OPTIMIZATION: Pre-calculate conversion factor
  days_per_decade_factor <- 3652.5 / 12  # 304.375

  # slope
  x4 <- terra::mean((x3 - x2) / 2, na.rm = TRUE)

  # 3. seasonal shifts (month/year) converted to days per decade
  sShift <- (trend / x4) * days_per_decade_factor
  sShift[sShift == Inf | sShift == -Inf] <- NA

  r2 <- c(trend, x4, sShift)
  names(r2) <- c("mTrend", "seaRate", "seaShift")
  return(r2)
}

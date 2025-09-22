#' Climate velocity trajectory spatial lines
#'
#' Create a spatial line data frame object from trajectory points.
#'
#' @param x \code{data.frame} containing the coordinates (x, y) of the constituent
#' points and identification number (trajIDs) for each trajectory as returned by VoCCTraj.
#' @param projx \code{CRS} detailing the coordinate reference system of the input data
#' (default geographic CRS).
#'
#' @return A \code{sf} object with one line per trajectory as specified in x.
#' To avoid artifacts, trajectories crossing the date line need to be split into two segments.
#' Where the trajectory on one side of the date line is only composed of a single point,
#' the trajectory won't be displayed (no line object created). The function assumes
#' a -180 to 180 longitudinal arrangement.
#'
#' @seealso{\code{\link{voccTraj}}}
#'
#' @export
#' @author Jorge Garcia Molinos
#' @examples
#' \dontrun{
#' HSST <- VoCC_get_data("HSST.tif")
#'
#' yrSST <- sumSeries(HSST,
#'   p = "1969-01/2009-12", yr0 = "1955-01-01", l = terra::nlyr(HSST),
#'   fun = function(x) colMeans(x, na.rm = TRUE), freqin = "months", freqout = "years"
#' )
#' tr <- tempTrend(yrSST, th = 10)
#' sg <- spatGrad(yrSST, th = 0.0001, projected = FALSE)
#' v <- gVoCC(tr, sg)
#' vel <- v[[1]]
#' ang <- v[[2]]
#'
#' # calculate the annual SST mean over the period
#' mn <- terra::mean(yrSST, na.rm = TRUE)
#'
#' # get the set of starting cells for the trajectories
#' lonlat <- stats::na.omit(data.frame(
#'   terra::xyFromCell(vel, 1:terra::ncell(vel)),
#'   vel[], ang[], mn[]
#' ))[, 1:2]
#'
#' # Calculate trajectories.
#' traj <- voccTraj(lonlat, vel, ang, mn, tyr = 50, correct = TRUE)
#'
#' # create a spatial line data frame from traj
#' lns <- trajLine(x = traj)
#' terra::plot(mn)
#' terra::plot(lns, add = TRUE)
#'
#' # Export as ESRI shape file
#' terra::writeVector(lns, filename = "velTraj", filetype = "ESRI Shapefile")
#' }
trajLine <- function(x, projx = "EPSG:4326") {

  x %>%
    dplyr::group_split(.data$ID) %>%
    furrr::future_map(get_trajLine, proj_x = projx,
                      .options = furrr::furrr_options(seed = TRUE),
                      .progress = TRUE) %>%
    purrr::list_rbind() %>%
    sf::st_sf()

}


#' @noRd
get_trajLine <- function(x, proj_x){

  # Get distance between first and last points
  d <- dplyr::bind_rows(dplyr::slice(x, 1), dplyr::slice(x, dplyr::n())) |>
    sf::st_as_sf(coords = c("lon", "lat")) |>
    sf::st_set_crs(proj_x) |>
    sf::st_distance() |>
    max()

  # Get number of steps; anything <240 means that the tracer died on land
  steps <- max(x$Steps)

  # Get distance along the tracer
  line_string <- x |>
    sf::st_as_sf(coords = c("lon", "lat")) |>
    sf::st_set_crs(proj_x) |>
    dplyr::summarise(do_union = FALSE) |>
    sf::st_cast(to = "LINESTRING")

# Seperate this out to allow use of sf::st_length(line_string)
    line_string |>
    dplyr::mutate(steps = steps,
                  line_distance = d,
                  line_length = sf::st_length(line_string),
                  ID = dplyr::first(x$ID),
                  lon = dplyr::first(x$lon),
                  lat = dplyr::first(x$lat))

}

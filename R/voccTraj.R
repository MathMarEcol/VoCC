#' Climate velocity trajectories
#'
#' Function to calculate vocc trajectories after Burrows et al (2014). Trajectories
#' are calculated by propagating climatic isopleths using the magnitude and direction of
#' local (cell) velocities. This is a slightly modified version of the original
#' Burrows et al. (2014) approach in that iterations of a trajectory are based on
#' cumulative time traveled instead of using fixed time steps.
#'
#' @param lonlat \code{data.frame} with the longitude and latitude (in decimal degrees)
#' of the points to project.
#' @param vel \code{raster} with the magnitude of gradient-based climate velocity.
#' @param ang \code{raster} with velocity angles in degrees.
#' @param mn \code{raster} with the overall mean climatic value over the period of interest.
#' @param x_res Numeric. Resolution of the grid in longitude direction (degrees or km).
#' @param y_res Numeric. Resolution of the grid in latitude direction (degrees or km).
#' @param tstep Numeric. Timestep for each trajectory iteration (usually decimal year).
#' @param tyr Integer. Temporal length of the period of interest (years).
#' @param grid_resolution Character. "coarse" (default) or "fine". Controls how land crossings are handled and allows for higher resolution grids.
#'
#' @return a \code{tibble} containing the coordinates ("lon", "lat") of the constituent
#' points, time step ("Steps"), identification number ("ID") for each trajectory, and cell IDs for start and end cells.
#'
#' @references \href{https://www.nature.com/articles/nature12976}{Burrows et al. 2014}. Geographical limits to species-range shifts are suggested by climate velocity. Nature, 507, 492-495.
#'
#' @seealso{\code{\link{gVoCC}}, \code{\link{trajClas}}}
#' @export
#' @author Jorge Garcia Molinos, David S. Schoeman and Michael T. Burrows
#' @examples
#' \dontrun{
#' yrSST <- sumSeries(HSST,
#'   p = "1960-01/2009-12", yr0 = "1955-01-01", l = terra::nlyr(HSST),
#'   fun = function(x) colMeans(x, na.rm = TRUE),
#'   freqin = "months", freqout = "years"
#' )
#'
#' # Long-term local climatic trends
#' tr <- tempTrend(yrSST, th = 10)
#'
#' # Local spatial climatic gradients
#' sg <- spatGrad(yrSST, th = 0.0001, projected = FALSE)
#'
#' # Gradient-based climate velocity
#' v <- gVoCC(tr, sg)
#' vel <- v[[1]]
#' ang <- v[[2]]
#'
#' # Calculate the annual SST mean over the period
#' mn <- terra::mean(yrSST, na.rm = TRUE)
#'
#' # Get the set of starting cells for the trajectories
#' lonlat <- stats::na.omit(data.frame(
#'   terra::xyFromCell(vel, 1:terra::ncell(vel)),
#'   vel[], ang[], mn[]
#' ))[, 1:2]
#'
#' # Calculate trajectories
#' # The following throws an error due to the trajectories moving beyond the raster extent
#' traj <- voccTraj(lonlat, vel, ang, mn, tyr = 50)
#'
#' # This accounts for the extent issue
#' traj <- voccTraj(lonlat, vel, ang, mn, tyr = 50, correct = TRUE)
#' }
#'
voccTraj <- function(lonlat, # Starting points
                     vel, ang, mn, # Components of velocity: speed, angle and climatology
                     x_res, y_res, # Resolution of the grid at which velocity was computed
                     tstep, # Timestep (usually decimal year)
                     tyr = 20, # Number of years to run for
                     grid_resolution = "coarse") { # Set to "fine" if you have disaggregated to original velocity field to a finer resolution





  # Setup -------------------------------------------------------------------

  # A base raster with original resolution and extent
  r_base <- terra::rast(res = c(x_res, y_res)) %>%
    terra::crop(vel)
  # Constrain max velocity to avoid stepping over grid squares
  max_vel <- 111.325 * x_res / tstep
  vel[vel > max_vel] <- max_vel
  vel[vel < -max_vel] <- -max_vel
  # Sort out start points
  lonlat <- lonlat %>%
    dplyr::select("x", "y") %>% # Collect just lon and lat (in case there's anything else there)
    as.data.frame()
  # Get initial descriptors
  tcells <- terra::cellFromXY(vel, lonlat) # # Cell IDs of starting cells
  n <- nrow(lonlat) # Get number of cells in your sequence
  # Set up variables to catch results
  llon <- llat <- cellIDend <- Steps <- cellIDs <- trajIDs <- list() # Initiate lists to catch results
  # Populate the first slots with starting points
  cellIDs[[1]] <- tcells
  trajIDs[[1]] <- 1:n
  llon[[1]] <- lonlat[, 1]
  llat[[1]] <- lonlat[, 2]
  cellIDend[[1]] <- tcells
  Steps[[1]] <- rep(0, n)
  # Set up objects that keep track of things
  sflags <- rep(NA, n)
  trj_id <- trajIDs[[1]]


  # Loop --------------------------------------------------------------------

  # Loop through the trajectories
  for (i in 1:(tyr / tstep)) { # 1:(tyr/tstep)
    llold <- lonlat # Take a copy of lonlat

    vell <- terra::extract(vel, llold) %>%
      dplyr::pull(.data$voccMag)
    angg <- terra::extract(ang, llold) %>%
      dplyr::pull(.data$voccAng)

    fcells <- terra::cellFromXY(vel, llold) # Get the cells that the trajectories start in
    # Get new locations
    lonlat <- destcoords(vell, angg, tstep, llold, y_res, x_res) # Extract lon and lat of landing point
    tcells <- terra::cellFromXY(vel, lonlat) # Get the cells that the trajectories end in
    sflags[which(is.na(tcells))] <- 1 # Sets the trajectory to "stuck" if it exits the velocity field (i.e., tcells == NA)
    # Remove out-of-bounds cells
    in_bounds <- which(is.na(sflags))
    llold <- llold[in_bounds, ]
    fcells <- fcells[in_bounds]
    lonlat <- lonlat[in_bounds, ]
    tcells <- tcells[in_bounds]
    sflags <- sflags[in_bounds]
    trj_id <- trj_id[in_bounds]
    # Deal with cells that end on land
    onland <- which(is.na(vel[tcells]))
    if (length(onland) > 0) { # Only bother if there is at least one cell that returns a NA velocity = land
      if (grid_resolution == "fine") { #*** Fine switch
        # Collect the stuff we need here for cells that are onland
        fpos <- llold[onland, ]
        tpos <- lonlat[onland, ]
        SFlags <- sflags[onland]
        fcell <- fcells[onland]
        tcell <- tcells[onland]
        ft <- data.frame(fcell = fcell, tcell = tcell, fx = purrr::pluck(fpos, 1), fy = purrr::pluck(fpos, 2), code = paste(fcell, tcell, sep = " "), ref = 1:length(fcell))
        ttcell <- apply(ft[, 1:4], 1, get_dest_cell_fine) #*** fine switch
        # Filter "stuck" flags here
        stuck <- which(is.na(ttcell)) # This is set in get_dest_cell(), where no cell in the "catchment" has a suitable sst to facilitate movement
        unstuck <- which(!is.na(ttcell))
        SFlags[stuck] <- 1 # Adjust flags
        # Make data frame to catch data needed to find new positions
        ttpos <- data.frame(x = rep(NA, length(onland)), y = rep(NA, length(onland)))
        ttpos[stuck, ] <- fpos[stuck, ] # If they're stuck, pass on starting points
        ttcell[stuck] <- fcell[stuck] # If they're stuck, pass on starting cells

        ttpos[unstuck, ] <- terra::xyFromCell(mn, ttcell[unstuck])
        # Collect results
        lonlat[onland, ] <- ttpos
        tcells[onland] <- ttcell
        sflags[onland] <- SFlags
      } else {
        # Old onland loop here
        # Collect the stuff we need here
        fpos <- llold[onland, ]
        tpos <- lonlat[onland, ]
        SFlags <- sflags[onland]
        fcell <- fcells[onland]
        tcell <- tcells[onland]
        # ft <- data.frame(fcell = fcell, tcell = tcell, code = paste(fcell, tcell, sep = " "), ref = 1:length(fcell))
        # ttcell <- apply(ft[,1:2], 1, get_dest_cell)
        ft <- data.frame(
          fcell = fcell, tcell = tcell,
          fx = fpos %>% purrr::pluck(1),
          fy = fpos %>% purrr::pluck(2),
          code = paste(fcell, tcell, sep = " "), ref = 1:length(fcell)
        )
        ttcell <- apply(ft[, 1:4], 1, get_dest_cell_coarse)
        # Filter "stuck" flags here
        stuck <- which(is.na(ttcell)) # This is set in get_dest_cell(), where no cell in the "catchment" has a suitable sst to facilitate movement
        unstuck <- which(!is.na(ttcell))
        SFlags[stuck] <- 1 # Adjust flags
        # Make data frame to catch data needed to find new positions #***done up to here
        ttpos <- data.frame(x = rep(NA, length(onland)), y = rep(NA, length(onland)))
        ttpos[stuck, ] <- fpos[stuck, ] # If they're stuck, pass on starting points
        ttcell[stuck] <- fcell[stuck] # If they're stuck, pass on starting cells
        if (length(unstuck) > 0) {
          tt_original_cells <- terra::cellFromXY(mn, fpos[unstuck, ]) # Departure cells in the resolution of the original velocity field
          ttdat <- tibble::tibble(ttcell = ttcell[unstuck]) %>% # Destination cells (nearest cell with appropriate sst)
            dplyr::mutate(
              loncent = terra::xFromCell(mn, ttcell), # Coordinates of destination cell
              latcent = terra::yFromCell(mn, ttcell)
            ) %>% # Coordinates of destination cell
            dplyr::mutate(e = NA, w = NA, n = NA, s = NA, dlon = NA, dlat = NA) # To facilitate finding corners of the cells
          corner_block_size <- 0.25 * x_res # The "corner" is set to 0.25 of the grid square at the original resolution
          # Send trajectory to the nearest corner of the appropriate cell, where corner is a quarter of the grid size. Position is "fuzzed" within this corner at random.
          ttdat$e <- ttdat$loncent + (0.5 * x_res) - (stats::runif(1, -1, 1) * corner_block_size) # The centre of the "corner" +- random fuzz...
          ttdat$w <- ttdat$loncent - (0.5 * x_res) + (stats::runif(1, -1, 1) * corner_block_size) # The centre of the "corner" +- random fuzz...
          ttdat$n <- ttdat$latcent + (0.5 * y_res) - (stats::runif(1, -1, 1) * corner_block_size) # The centre of the "corner" +- random fuzz...
          ttdat$s <- ttdat$latcent - (0.5 * y_res) + (stats::runif(1, -1, 1) * corner_block_size) # The centre of the "corner" +- random fuzz...
          coords <- with(ttdat, cbind(n, e, n, w, s, w, s, e)) # NE, NW, SW, SE corners' coordinates

          corners <- data.frame(ne = get_dist(fpos[unstuck, 2], fpos[unstuck, 1], coords[, 1], coords[, 2])) %>%
            dplyr::mutate(
              nw = get_dist(fpos[unstuck, 2], fpos[unstuck, 1], coords[, 3], coords[, 4]),
              sw = get_dist(fpos[unstuck, 2], fpos[unstuck, 1], coords[, 5], coords[, 6]),
              se = get_dist(fpos[unstuck, 2], fpos[unstuck, 1], coords[, 7], coords[, 8])
            )
          # Select nearest corner
          cornset <- apply(corners, 1, mfind) * 2 # Identify which corners for each onland cell. Have to mul by 2 to shift along correctly.
          cornset <- cbind(cornset, coords) # Add in coordinates
          ttdat[, 8:9] <- data.frame(t(apply(cornset, 1, mplace))) # Add in coordinates of correct corner point
          ttpos[unstuck, ] <- terra::xyFromCell(mn, ttcell[unstuck])
          ttcell[unstuck] <- terra::cellFromXY(mn, ttpos[unstuck, ])
        }
        # Collect results
        lonlat[onland, ] <- ttpos
        tcells[onland] <- ttcell
        sflags[onland] <- SFlags
      }
    }
    # Pass on only those cells that are not stuck
    # if(sum(is.na(lonlat[,1])) > 0) {break}
    cells_to_keep <- which(is.na(sflags))
    lonlat <- lonlat[cells_to_keep, ]
    # Collect results
    llon[[i + 1]] <- lonlat[, 1] # Add lon to the list
    llat[[i + 1]] <- lonlat[, 2] # Add lat to the list
    cellIDs[[i + 1]] <- fcells[cells_to_keep] # Add departure cellIDs to the list
    cellIDend[[i + 1]] <- tcells[cells_to_keep] # Add arrival cellIDs to the list
    Steps[[i + 1]] <- rep(i, length(tcells[cells_to_keep])) # Capture the time_step
    trajIDs[[i + 1]] <- trj_id[cells_to_keep]
    sflags <- sflags[cells_to_keep] # Adjust the flags list
    trj_id <- trj_id[cells_to_keep] # Adjust the trajectory IDs
    print(c(i, length(onland), round(100 * i / (tyr / tstep), 3), base::date())) # Keep a reference of progress
  }
  return(cbind(
    Steps = purrr::reduce(Steps, c),
    lon = purrr::reduce(llon, c),
    lat = purrr::reduce(llat, c),
    ID = purrr::reduce(trajIDs, c),
    start_cell = purrr::reduce(cellIDs, c),
    end_cell = purrr::reduce(cellIDend, c)
  ) %>%
    tibble::as_tibble())
}

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
#' @param bfr Numeric. Distance around a cell to look for a cooler/warmer cell.
#' @param grid_resolution Character. "coarse" (default) or "fine". Controls how land crossings are handled and allows for higher resolution grids.
#' @param seed Integer. Random seed for reproducibility. If NULL (default), no seed is set.
#'
#'
#' @return a \code{tibble} containing the coordinates ("lon", "lat") of the constituent
#' points, time step ("Steps"), identification number ("ID") for each trajectory, and cell IDs for start and end cells.
#'
#' @references \href{https://www.nature.com/articles/nature12976}{Burrows et al. 2014}. Geographical limits to species-range shifts are suggested by climate velocity. Nature, 507, 492-495.
#'
#' @seealso{\code{\link{gVoCC}}, \code{\link{trajClas}}}
#' @export
#' @examples
#' \dontrun{
#' HSST <- terra::rast(system.file("extdata", "HadiSST.tif", package = "VoCCdata"))
#'
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
#' traj <- voccTraj(lonlat, vel, ang, mn, tstep = 1/4, tyr = 50)
#' }
#'
voccTraj <- function(lonlat, # Starting points
                     vel, ang, mn, # Components of velocity: speed, angle and climatology
                     tstep, # Timestep (usually decimal year)
                     x_res = NULL, y_res = NULL, # Resolution of the grid at which velocity was computed
                     tyr = 20, # Number of years to run for
                     bfr = 75,
                     grid_resolution = "coarse", # Set to "fine" if you have disaggregated to original velocity field to a finer resolution
                     seed = NULL) { # Random seed for reproducibility


  # Setup -------------------------------------------------------------------

  # Validation checks using assertthat
  assertthat::assert_that(
    # lonlat must be a data.frame with at least 2 columns
    is.data.frame(lonlat),
    msg = "lonlat must be a data.frame"
  )

  assertthat::assert_that(
    ncol(lonlat) >= 2,
    msg = "lonlat must have at least 2 columns (longitude and latitude)"
  )

  assertthat::assert_that(
    nrow(lonlat) > 0,
    msg = "lonlat must contain at least one row"
  )

  assertthat::assert_that(
    # vel, ang, and mn must be terra rasters
    inherits(vel, "SpatRaster"),
    msg = "vel must be a terra SpatRaster"
  )

  assertthat::assert_that(
    inherits(ang, "SpatRaster"),
    msg = "ang must be a terra SpatRaster"
  )

  assertthat::assert_that(
    inherits(mn, "SpatRaster"),
    msg = "mn must be a terra SpatRaster"
  )

  assertthat::assert_that(
    # tstep must be positive numeric
    is.numeric(tstep) && length(tstep) == 1 && tstep > 0,
    msg = "tstep must be a positive numeric value"
  )

  assertthat::assert_that(
    # x_res and y_res must be positive if provided
    is.null(x_res) || (is.numeric(x_res) && length(x_res) == 1 && x_res > 0),
    msg = "x_res must be a positive numeric value if provided"
  )

  assertthat::assert_that(
    is.null(y_res) || (is.numeric(y_res) && length(y_res) == 1 && y_res > 0),
    msg = "y_res must be a positive numeric value if provided"
  )

  assertthat::assert_that(
    # When either x_res or y_res is given, the other must also be given
    (is.null(x_res) && is.null(y_res)) || (!is.null(x_res) && !is.null(y_res)),
    msg = "When either x_res or y_res is provided, both must be provided"
  )

  assertthat::assert_that(
    # tyr must be positive numeric
    is.numeric(tyr) && length(tyr) == 1 && tyr > 0,
    msg = "tyr must be a positive numeric value"
  )

  assertthat::assert_that(
    # bfr must be positive numeric
    is.numeric(bfr) && length(bfr) == 1 && bfr > 0,
    msg = "bfr must be a positive numeric value"
  )

  assertthat::assert_that(
    # grid_resolution can only be "coarse" or "fine"
    is.character(grid_resolution) && length(grid_resolution) == 1 &&
    grid_resolution %in% c("coarse", "fine"),
    msg = "grid_resolution must be either 'coarse' or 'fine'"
  )

  assertthat::assert_that(
    # When both x_res and y_res are NULL, grid_resolution must be "coarse"
    !(is.null(x_res) && is.null(y_res)) || grid_resolution == "coarse",
    msg = "When both x_res and y_res are NULL, grid_resolution must be 'coarse'"
  )

  assertthat::assert_that(
    # seed must be numeric if provided
    is.null(seed) || (is.numeric(seed) && length(seed) == 1),
    msg = "seed must be a numeric value if provided"
  )

  if (is.null(x_res) | is.null(y_res)){
    vel_res <- terra::res(vel)

    x_res <- vel_res[1]
    y_res <- vel_res[2]

  }


  # Set seed for reproducibility if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # A base raster with original resolution and extent
  r_base <- terra::rast(res = c(x_res, y_res)) %>%
    terra::crop(vel)

  # MEMORY LEAK FIX: Don't modify original raster, work with values directly
  max_vel <- 111.325 * x_res / tstep
  vel_values <- terra::values(vel)
  vel_values[vel_values > max_vel] <- max_vel
  vel_values[vel_values < -max_vel] <- -max_vel

  # Sort out start points
  lonlat <- lonlat %>%
    dplyr::select("x", "y") %>% # Collect just lon and lat (in case there's anything else there)
    as.data.frame()

  # Get initial descriptors
  tcells <- terra::cellFromXY(vel, lonlat) # # Cell IDs of starting cells
  n <- nrow(lonlat) # Get number of cells in your sequence

  # OPTIMIZATION: Use the already constrained vel_values and extract ang_values
  ang_values <- terra::values(ang)
  mn_values <- terra::values(mn)  # Pre-extract mn values too

  # MEMORY OPTIMIZATION: Pre-allocate list structure (not content) to avoid dynamic growth
  # This creates a list of NULL pointers - no contiguous memory required
  max_steps <- ceiling(tyr / tstep) + 1  # Maximum possible steps + initial

  # Pre-allocate list structure only - each element will be allocated independently
  llon <- vector("list", max_steps)
  llat <- vector("list", max_steps)
  cellIDend <- vector("list", max_steps)
  Steps <- vector("list", max_steps)
  cellIDs <- vector("list", max_steps)
  trajIDs <- vector("list", max_steps)

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
  n_steps <- ceiling(tyr / tstep)

  # Add safety check for reasonable number of steps
  if (n_steps > 1000) {
    warning("Large number of trajectory steps (", n_steps, "). Consider reducing tyr or increasing tstep.")
  }

  # Track actual steps used for efficient final processing
  actual_steps_used <- 1  # Start with 1 (initial step)

  for (i in seq_len(n_steps)) {

    # Safety check: if no trajectories remain, break early
    if (nrow(lonlat) == 0) {
      message("All trajectories terminated at step ", i-1)
      break
    }

    llold <- lonlat # Take a copy of lonlat

    # OPTIMIZATION: Get cell IDs first, then extract values by indexing (much faster)
    fcells <- terra::cellFromXY(vel, llold) # Get the cells that the trajectories start in

    # MEMORY LEAK FIX: Add bounds checking to prevent invalid indexing
    valid_fcells <- !is.na(fcells) & fcells > 0 & fcells <= length(vel_values)
    if (!any(valid_fcells)) {
      message("All trajectories moved out of bounds at step ", i)
      break
    }

    # Extract values by direct indexing (much faster than terra::extract)
    vell <- vel_values[fcells[valid_fcells]]
    angg <- ang_values[fcells[valid_fcells]]

    # Update tracking variables for valid cells only
    if (!all(valid_fcells)) {
      llold <- llold[valid_fcells, , drop = FALSE]
      fcells <- fcells[valid_fcells]
      trj_id <- trj_id[valid_fcells]
      sflags <- sflags[valid_fcells]
    }

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

    # MEMORY LEAK FIX: Use pre-extracted values instead of raster indexing
    onland <- which(is.na(vel_values[tcells]))

    if (length(onland) > 0) { # Only bother if there is at least one cell that returns a NA velocity = land
      if (grid_resolution == "fine") { #*** Fine switch

        # Collect the stuff we need here for cells that are onland
        fpos <- llold[onland, ]
        tpos <- lonlat[onland, ]
        SFlags <- sflags[onland]
        fcell <- fcells[onland]
        tcell <- tcells[onland]

        ft <- data.frame(fcell = fcell,
                         tcell = tcell,
                         fx = purrr::pluck(fpos, 1),
                         fy = purrr::pluck(fpos, 2),
                         code = paste(fcell, tcell, sep = " "),
                         ref = 1:length(fcell))

        # MEMORY LEAK FIX: Pass values instead of full raster objects
        ttcell <- apply(ft[, 1:4], 1, get_dest_cell_fine, x_res = x_res, y_res = y_res, bfr = bfr,
                       vel_values = vel_values, mn_values = mn_values, vel_raster = vel, mn_raster = mn)

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
        # MEMORY LEAK FIX: Pass values instead of full raster objects
        ttcell <- apply(ft[, 1:4], 1, get_dest_cell_coarse, x_res = x_res, y_res = y_res, bfr = bfr,
                       vel_values = vel_values, mn_values = mn_values, vel_raster = vel, mn_raster = mn)

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
          # Use reproducible random numbers for package checking
          n_points <- nrow(ttdat) #TODO Check that this is still 1. Not sure why we need it.
          ttdat$e <- ttdat$loncent + (0.5 * x_res) - (stats::runif(n_points, -1, 1) * corner_block_size) # The centre of the "corner" +- random fuzz...
          ttdat$w <- ttdat$loncent - (0.5 * x_res) + (stats::runif(n_points, -1, 1) * corner_block_size) # The centre of the "corner" +- random fuzz...
          ttdat$n <- ttdat$latcent + (0.5 * y_res) - (stats::runif(n_points, -1, 1) * corner_block_size) # The centre of the "corner" +- random fuzz...
          ttdat$s <- ttdat$latcent - (0.5 * y_res) + (stats::runif(n_points, -1, 1) * corner_block_size) # The centre of the "corner" +- random fuzz...
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

    # MEMORY OPTIMIZATION: Direct assignment to pre-allocated list slots
    # Each element is stored independently - no contiguous memory required
    step_index <- i + 1
    llon[[step_index]] <- lonlat[, 1]
    llat[[step_index]] <- lonlat[, 2]
    cellIDs[[step_index]] <- fcells[cells_to_keep]
    cellIDend[[step_index]] <- tcells[cells_to_keep]
    Steps[[step_index]] <- rep(i, length(tcells[cells_to_keep]))
    trajIDs[[step_index]] <- trj_id[cells_to_keep]

    # Update tracking variables
    sflags <- sflags[cells_to_keep]
    trj_id <- trj_id[cells_to_keep]
    actual_steps_used <- step_index

    # Progress reporting - only in interactive sessions or when explicitly requested
    if (interactive() && getOption("VoCC.verbose", FALSE)) {
      message("Step ", i, "/", tyr/tstep, " (", round(100 * i / (tyr / tstep), 1), "%) - ",
              length(onland), " cells on land - ", Sys.time())
    }

  }

  # MEMORY LEAK FIX: Clean up large objects before final processing
  rm(vel_values, ang_values, mn_values)

  # MEMORY OPTIMIZATION: Only process used slots and clean up progressively
  # Trim to actual used length - unused slots remain as NULL (minimal memory)
  if (actual_steps_used < max_steps) {
    llon <- llon[1:actual_steps_used]
    llat <- llat[1:actual_steps_used]
    Steps <- Steps[1:actual_steps_used]
    trajIDs <- trajIDs[1:actual_steps_used]
    cellIDs <- cellIDs[1:actual_steps_used]
    cellIDend <- cellIDend[1:actual_steps_used]
  }

  # MEMORY LEAK FIX: Progressive cleanup to minimize peak memory usage
  # Each unlist operation works on independent memory chunks
  steps_vec <- unlist(Steps, use.names = FALSE)
  rm(Steps)

  lon_vec <- unlist(llon, use.names = FALSE)
  rm(llon)

  lat_vec <- unlist(llat, use.names = FALSE)
  rm(llat)

  id_vec <- unlist(trajIDs, use.names = FALSE)
  rm(trajIDs)

  start_cell_vec <- unlist(cellIDs, use.names = FALSE)
  rm(cellIDs)

  end_cell_vec <- unlist(cellIDend, use.names = FALSE)
  rm(cellIDend)

  # Create final result with cleaned vectors
  result <- tibble::tibble(
    Steps = steps_vec,
    lon = lon_vec,
    lat = lat_vec,
    ID = id_vec,
    start_cell = start_cell_vec,
    end_cell = end_cell_vec
  )

  # Clean up final vectors
  rm(steps_vec, lon_vec, lat_vec, id_vec, start_cell_vec, end_cell_vec)

  # Clean up any remaining large objects from the loop
  if (exists("lonlat")) rm(lonlat)
  if (exists("llold")) rm(llold)
  if (exists("fcells")) rm(fcells)
  if (exists("tcells")) rm(tcells)
  if (exists("sflags")) rm(sflags)
  if (exists("trj_id")) rm(trj_id)

  # Force garbage collection to free memory immediately
  gc()

  return(result)
}

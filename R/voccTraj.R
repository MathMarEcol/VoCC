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
#' @param vel_c
#' @param ang_c
#' @param mn_c
#' @param fine_coast
#' @param lk_up
#' @param tstep
#' @param tyr \code{integer} temporal length of the period of interest.
#'
#' @return a \code{data.frame} containing the coordinates ("x", "y") of the constituent
#' points and identification number ("trajIDs") for each trajectory.
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
voccTraj <- function(lonlat,
                     vel, ang, mn,
                     # vel_c, ang_c, mn_c,
                     #fine_coast, lk_up,
                     tstep, tyr = 50) {

  # Setup -------------------------------------------------------------------

  y_res <- x_res <- terra::res(vel)[1] # Set resolution of operations

  # Constrain max velocity to avoid stepping over grid squares
  max_vel = 111.325*x_res/tstep
  vel[vel > max_vel] <- max_vel
  vel[vel < -max_vel] <- -max_vel
  # vel_c[vel_c > max_vel] <- max_vel
  # vel_c[vel_c < -max_vel] <- -max_vel

  lonlat <- lonlat %>%
    dplyr::select(tidyselect::all_of(c("x", "y"))) %>%  # Collect just lon and lat (in case there's anything else there)
    as.data.frame()

  tcells <- terra::cellFromXY(vel, lonlat) # # Cell IDs of starting cells
  n <- nrow(lonlat) # Get number of cells in your sequence
  sflags <- rep(NA, n) # Set a string of cells that change from NA to 1 where the trajectory sticks

  # Set up variables to catch results, allocating the right amount of memory
  llon <- numeric((n * tyr / tstep) + n) # Needs to have a set of starting lons, plus one more set for each time step
  llat <- numeric((n * tyr / tstep) + n) # Needs to have a set of starting lats, plus one more set for each time step
  # cellIDs <- rep(tcells, ((tyr / tstep) + 1)) # Needs to have a set of starting cells, plus one more set for each time step—just as a reference
  cellIDs <- rep(1:n, ((tyr / tstep) + 1)) # Needs to have a set of starting cells, plus one more set for each time step—just as a reference
  cellIDend <- numeric((n * tyr / tstep) + n) # Needs to have a set of ending cells, plus one more set for each time step
  # coast <- llon != 0 # Needs to have a set of starting coastal flags, plus one more set for each time step...set up as boolean
  flags <- numeric((n * tyr / tstep) + n) # Needs to have a set of starting flags, plus one more set for each time step
  Steps <- numeric((n * tyr / tstep) + n) # Needs to have a set of starting steps, plus one more set for each time step

  # Populate the first n slots with starting points
  llon[1:n] <- lonlat[,1]
  llat[1:n] <- lonlat[,2]
  cellIDend[1:n] <- tcells
  # coast[1:n] <- scoast
  flags[1:n] <- NA
  Steps[1:n] <- 0

  # Get coarse cellIDs that exist in fine raster
  # fn_lk_up <- lk_up[] %>%
  #   na.omit() %>%
  #   unique()

  # Helper functions --------------------------------------------------------

  # Trajectory helper functions

  # Function for to grep min distance
  mfind <- function(rw){
    X <- which.min(rw)
    return(X)
  }

  # Function for to extract corner coordinates per row
  mplace <- function(rw){
    X <- rw[(rw[1]+1)]
    Y <- rw[rw[1]]
    return(c(X, Y))
  }

  # Simple circular functions
  deg2rad <- function(deg) return(deg*pi/180)
  rad2deg <- function(rad) return(180*rad/pi)

  # A function to find the nearest coastal cell
  get_nearest_coastal_cell <- function(x, y, tccell, lk_up, ...) { # JDE added lk_up to arguments
    t_block <- which(lk_up[] == tccell) # The cellIDs of lk_up
    cst_pts <- terra::xyFromCell(lk_up, t_block) %>%
      as.data.frame() %>%
      sf::st_as_sf(coords = c("x", "y"), crs = terra::crs(terra::rast()))
    # Here, we can't search for corners and add random "fuzz" to we collect 10 random points (if there are 10), and select the nearest of those, instead. If there are 10 or fewer, just pick the nearest
    if(nrow(cst_pts) > 10) {
      cst_pts <- cst_pts[sample(1:nrow(cst_pts), 10, replace = FALSE),]
    } # The unspoken "else" is that we just retain the cst_pts we have
    pt <- sf::st_as_sf(data.frame(x = x, y = y), coords = c("x", "y"), crs = terra::crs(terra::rast())) # The point departed from
    nearest <- sf::st_distance(cst_pts, pt) %>%
      which.min()
    out <- cst_pts[nearest,] %>%
      sf::st_coordinates() %>%
      as.data.frame() %>%
      dplyr::rename(x = 1, y = 2)
    return(out)
  }

  # Find new destination, given velocity (km/yr), angle (º), time step (yr) and initial coordinates (ºlon, ºlat); max allowed jump is 1 deg
  # vell = vel[fcells] %>% pull(1); angg = ang[fcells] %>% pull(1); timestep = tstep; ll = llold
  destcoords <- function(vell, angg, timestep, ll, y_res, x_res){

    latshift <- (abs(vell) * timestep * cos(deg2rad(angg))) / 111.325 # Calculate shift in lat
    latnew <- ll[,2] + latshift # Find new lat...first approximation
    lonshift <- (abs(vell) * timestep * sin(deg2rad(angg))) / (111.325 * cos(deg2rad(latnew))) # Shift in lon
    # Limit large longitudinal jumps at high latitudes

    # Because we constrain velocity to be no more than 12 * y_res, all problems will be with lonshift
    # Limit lonshift to at most 1 cell
    lonshift[lonshift > x_res] <- x_res
    lonshift[lonshift < -x_res] <- -x_res

    # Now on that basis, adjust latshift
    x_gt <- which(lonshift == x_res) # Indices for adjusted lon shifts
    latshift[x_gt] <- ((x_res*111.325 * cos(deg2rad(ll[x_gt,2])))/tan(deg2rad(angg[x_gt])))/111.325 # Using trig on distances
    x_lt <- which(lonshift == -x_res) # Indices for adjusted lon shifts
    latshift[x_lt] <- ((x_res*111.325 * cos(deg2rad(ll[x_lt,2])))/tan(deg2rad(angg[x_lt])))/111.325 # Using trig on distances
    latnew <- ll[,2] + latshift # Find new lat by adding the adjusted lats
    # Stop new lat from jumping the poles
    latnew[latnew > 90] <- 90
    latnew[latnew < -90] <- -90
    # Adjust lon
    lonnew <- ll[,1] + lonshift # Find new lon...first approximation
    # Adjust for dateline jumps
    lonnew <- lonnew - (360 * floor((lonnew + 180) / 360))
    return(data.frame(lonnew, latnew) %>%
             stats::setNames(c("dlon", "dlat")))
  }

  # Function for to find closest cooler/warmer cell within the union of two levels of 8-cell adjacency from the "departed" and "destination" cells
  # TODO JDE - The arguments are not correct. Only 1, but it should include. `mn`
  get_dest_cell <- function(rw) { # A to-from cell pair and the buffer around a cell on land that can be searched for a suitable cell


    x_res <- terra::res(mn)[1]
    y_res <- terra::res(mn)[2]

    # Find clumps of ocean; the rule is that you can't jump from one clump to another, because this would mean passing over unbroken land
    xy <- terra::xyFromCell(mn, as.vector(as.matrix(rw))) %>%
      as.data.frame()

    bfr = ((y_res*111.325) + 1) * 1000 # Set buffer to be 1 grid-square width at the equator + 1 km

    xy <- xy %>%
      sf::st_as_sf(coords = c("x", "y"), crs = "EPSG:4326")

    sp_buffer <- sf::st_buffer(sf::st_as_sf(xy), bfr) # Remembering that buffer is in metres

    buffer_zone <- terra::extract(mn, sp_buffer, cells = TRUE, xy = TRUE, touches = TRUE) %>%
      dplyr::select(-"ID") %>%
      dplyr::distinct() %>%
      dplyr::rename(sst = mean) %>% # JDE Changed from climatology to mean
      dplyr::select(tidyselect::all_of(c("x", "y", "sst", "cell"))) %>%
      tidyr::drop_na(.data$sst)

    # JDE Changes. When all cells are in a straight line (only 1 unique x or y),
    # rast() will not work. So here we add a more robust approach to creating
    # a rast() regardless of the data

    # clumped <- buffer_zone %>%
    #   dplyr::select(-cell) %>%
    #   rast(crs = "EPSG:4326") %>%
    #   terra::patches(directions = 8, allowGaps = FALSE)


    clumped_rast <- terra::rast(
      xmin = min(buffer_zone$x) - x_res/2,
      xmax = max(buffer_zone$x) + x_res/2,
      ymin = min(buffer_zone$y) - y_res/2,
      ymax = max(buffer_zone$y) + y_res/2,
      resolution = c(x_res, y_res),
      crs = "EPSG:4326"
    )

    clumped <- terra::rasterize(
      x = buffer_zone %>% # Needs to be SpatVector to add the sst
        dplyr::select(-"cell") %>%
        terra::vect(geom = c("x", "y"), crs = "EPSG:4326"),
      y = clumped_rast, # The template to rasterize onto
      field = "sst") %>% # The data
      terra::patches(directions = 8, allowGaps = FALSE)

    # JDE END CHANGES

    # Which clump did I start in?
    r1 <- rw[1] %>%
      as.vector() # We use this cell a lot, so let's just make it an object

    from_clump <- terra::extract(clumped, terra::xyFromCell(mn, r1)) %>%
      unlist() %>%
      as.vector()

    # What are the coordinates of cells within the search area that fall in the clump I came from?
    search_xy <- terra::xyFromCell(clumped, which(clumped[] == from_clump)) %>%
      as.data.frame() #****
    search_cells <- terra::cellFromXY(mn, search_xy) # Which cells are these

    # Check if any of these are NOT in in EITHER fine coast or mn, and eliminate, if needed
    to_keep <- c(which(!is.na(terra::extract(mn, search_xy, ID = FALSE) %>%
                                dplyr::pull(1)))) %>% #, # Coords in the coarse grid
      # which(search_cells %in% fn_lk_up)) %>%
      unique()

    search_xy1 <- search_xy[to_keep,] %>%
      as.data.frame() %>%
      dplyr::distinct()

    # Get the ssts in the cells to search
    or <- terra::extract(mn, search_xy1, cells = TRUE, xy = TRUE) %>%
      dplyr::rename(sst = 2, x = 4, y = 5)


    # If velocity is positive, find nearest cell that is cooler, if there is one
    if(unlist(vel[r1]) > 0) {
      o <- or %>%
        dplyr::filter(.data$sst < unlist(mn[r1])) %>%
        na.omit()
      if(nrow(o) == 0) {
        dest_cell <- NA # Set condition with which to ID stuck cells
      } else {
        pt <- terra::xyFromCell(mn, r1) %>%
          data.frame()

        dest_cell <- sf::st_distance(sf::st_as_sf(pt, coords = c("x", "y"), crs = terra::crs(terra::rast())),
                                     sf::st_as_sf(data.frame(o), coords = c("x", "y"), crs = terra::crs(terra::rast()))) %>%
          which.min() %>%
          o$cell[.]
      }
    } else { # Otherwise, find nearest cell that is warmer, if there is one
      o <- or %>%
        dplyr::filter(.data$sst > unlist(mn[r1])) %>%
        na.omit()
      if(nrow(o) == 0) {
        dest_cell <- NA # Set condition with which to ID stuck cells
      } else {
        pt <- terra::xyFromCell(mn, r1) %>%
          data.frame()

        dest_cell <- sf::st_distance(sf::st_as_sf(pt, coords = c("x", "y"), crs = terra::crs(terra::rast())),
                                     sf::st_as_sf(o, coords = c("x", "y"), crs = terra::crs(terra::rast()))) %>%
          which.min() %>%
          o$cell[.]
      }
    }
    return(dest_cell)
  }

  # Loop --------------------------------------------------------------------

  # Loop through the trajectories
  for(i in 1:(tyr/tstep)) { # 1:(tyr/tstep)
    llold <- lonlat # Take a copy of lonlat
    fcells <- terra::cellFromXY(vel,llold) # Get the cells that the trajectories start in

    # Pull velocity and angle from the "inland" versions, because we always ensure that we can find destination cells for departures within the coastal "blind spot"
    # TIN: changed this from vel_c to vel and ang_c to ang because we don't have coarse/fine versions
    vc <- vel[fcells] %>% dplyr::pull(1)
    ac <- ang[fcells] %>% dplyr::pull(1)
    llc <- llold

    # Get new locations
    lonlat <- destcoords(vc, ac, tstep, llc, y_res, x_res) # Extract lon and lat of landing point
    tcells <- terra::cellFromXY(vel,lonlat) # Get the cells that the trajectories end in
    sflags[which(is.na(tcells))] <-  1 # Sets the trajectory to "stuck" if it exits the velocity field (i.e., tcells == NA)

    # Bounce stuck cells
    stuck_cells <- which(!is.na(sflags))
    lonlat[stuck_cells,] <- llold[stuck_cells,] # Bounce the corresponding coordinates across for "stuck" cells
    tcells[stuck_cells] <- fcells[stuck_cells] # Bounce the corresponding cells across for "stuck" cells

    # Deal with cells that end on land
    onland <- which(is.na(vel[tcells])) %>%  # Identify which rows of velend are on land, provided that they are not stuck
      dplyr::setdiff(stuck_cells) # Ignoring stuck cells

    # if(length(onland) > 0){ # Only bother if there is at least one cell that returns a NA velocity = land
    # Here, you need to check whether it really *IS* on land, or whether it is just in the coastal "blind spot"
    # fn <- terra::extract(fine_coast, lonlat[onland,], ID = FALSE) %>%
    #   pull(1)
    # onland <- onland[which(is.na(fn))]
    if(length(onland) > 0) {
      # Collect the stuff we need here for cells that are really onland
      fpos <- llold[onland,]
      tpos <- lonlat[onland,]
      SFlags <- sflags[onland]
      fcell <- fcells[onland]
      tcell <- tcells[onland]
      ft <- data.frame(fcell = fcell, tcell = tcell, code = paste(fcell, tcell, sep = " "), ref = 1:length(fcell))

      # Get the destination cells of those that are hitting the coast
      ttcell <- apply(ft[,1:2], 1, get_dest_cell) # TODO using apply in places and map in others. Can we consistently move to purrr?

      # Filter "stuck" flags here
      stuck <- which(is.na(ttcell))
      unstuck <- which(!is.na(ttcell))
      SFlags[stuck] <- 1 # Adjust flags

      # Make data frame to catch data needed to find new positions
      ttpos <- data.frame(x = rep(NA, length(onland)), y = rep(NA, length(onland)))
      ttpos[stuck,] <- fpos[stuck,] # If they're stuck, pass on starting points
      ttcell[stuck] <- fcell[stuck] # If they're stuck, pass on starting cells

      if(length(unstuck) > 0) {
        ttdat <- data.frame(ttcell = ttcell[unstuck],
                            loncent = terra::xFromCell(vel, ttcell[unstuck]),
                            latcent = terra::yFromCell(vel, ttcell[unstuck])) %>%  # Start building coordinates
          dplyr::mutate(e = NA, w = NA, n = NA, s = NA, dlon = NA, dlat = NA) # To facilitate finding corners, if we need them

        # Separate cells in the coastal blind spot from those that are not
        # which_coast <- which(ttdat$ttcell %in% fn_lk_up)
        # which_not_coast <- which(!(ttdat$ttcell %in% fn_lk_up))
        # For non-coastal cells

        idx <- 1:nrow(ttdat)

        if(nrow(ttdat) > 0) {
          corner_block_size <- 0.25*x_res
          nc_ttdat <- ttdat # Copying data
          nc_ttdat$e <- nc_ttdat$loncent + (0.5 * x_res) - (stats::runif(1, 0, 1) * corner_block_size)
          nc_ttdat$w <- nc_ttdat$loncent - (0.5 * x_res) + (stats::runif(1, 0, 1) * corner_block_size)
          nc_ttdat$n <- nc_ttdat$latcent + (0.5 * y_res) - (stats::runif(1, 0, 1) * corner_block_size)
          nc_ttdat$s <- nc_ttdat$latcent - (0.5 * y_res) + (stats::runif(1, 0, 1) * corner_block_size)
          coords <- with(nc_ttdat, cbind(n, e, n, w, s, w, s, e)) # NE, NW, SW, SE corners' coordinates
          # Find distances from departure point to corners of ttcell
          get_dist <- function(y1, x1, y2, x2) {
            pt1 <- sf::st_as_sf(data.frame(x = x1, y = y1), coords = c("x", "y"), crs = terra::crs(terra::rast()))
            pt2 <- sf::st_as_sf(data.frame(x = x2, y = y2), coords = c("x", "y"), crs = terra::crs(terra::rast()))
            out <- sf::st_distance(pt1, pt2, by_element = TRUE) %>%
              as.vector()
            return(out)
          }

          corners <- data.frame(ne = get_dist(fpos[idx,2], fpos[idx,1], coords[,1], coords[,2])) %>%
            dplyr::mutate(nw = get_dist(fpos[idx,2], fpos[idx,1], coords[,3], coords[,4]),
                   sw = get_dist(fpos[idx,2], fpos[idx,1], coords[,5], coords[,6]),
                   se = get_dist(fpos[idx,2], fpos[idx,1], coords[,7], coords[,8]))
          cornset <- apply(corners, 1, mfind)*2 # Identify which corners for each onland cell. Have to mul by 2 to shift along correctly.
          cornset <- cbind(cornset, coords) # Add in coordinates
          ttdat[,8:9] <- data.frame(t(apply(cornset, 1, mplace))) # Extract coordinates of correct corner
        }
        # 	          # For coastal cells
        #               if(length(which_coast) > 0) {
        #                   coastal_dat <- ttdat[which_coast,2:3] %>%
        #                     mutate(x = fpos[which_coast,1], y = fpos[which_coast,2],
        #                            tccell = ttdat[which_coast,]$ttcell)
        #                   # JDE TODO - Argument to get_nearest_coast_cell are not correct - `get_nearest_coastal_cell <- function(x, y, tccell, lk_up, ...)`
        #                   ttdat[which_coast,8:9] <- pmap_df(coastal_dat, get_nearest_coastal_cell)
        #                 }
        ttpos[unstuck,] <- ttdat[,8:9]
        ttcell[unstuck] <- terra::cellFromXY(mn, ttpos[unstuck,])
      }
      # Collect results
      lonlat[onland,] <- ttpos
      tcells[onland] <- ttcell
      sflags[onland] <- SFlags
    }
    # }
    llon[((i * n) + 1): ((i * n) + n)] <- lonlat[,1] # Add lon to the list
    llat[((i * n) + 1): ((i * n) + n)] <- lonlat[,2]# Add lat to the list
    cellIDend[((i * n) + 1): ((i * n) + n)] <- tcells # Add cellIDs to the list
    flags[((i * n) + 1): ((i * n) + n)] <- sflags # Add flags to the list
    Steps[((i * n) + 1): ((i * n) + n)] <- i # Capture the time_step
    # print(c(i, length(onland), round(100*i/(tyr/tstep), 3), date())) # Keep a reference of progress

  }
  return(cbind(Steps, llon, llat, cellIDs, cellIDend, flags) %>%
           tibble::as_tibble())
}


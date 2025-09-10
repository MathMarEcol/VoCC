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
voccTraj <- function(lonlat, # Starting points
                                     vel, ang, mn, # Components of velocity: speed, angle and climatology
                                     x_res, y_res, # Resolution of the grid at which velocity was computed
                                     tstep, # Timestep (usually decimal year)
                                     tyr = 20, # Number of years to run for
                                     grid_resolution = "coarse") { # Set to "fine" if you have disaggregated to original velocity field to a finer resolution

  # Internal helper functions ----------------------------------------------------

  # Function for finding minimum distance
  mfind <- function(rw) {
    X <- which.min(rw)
    return(X)
  }

  # Function to extract corner coordinates per row
  mplace <- function(rw) {
    X <- rw[(rw[1]+1)]
    Y <- rw[rw[1]]
    return(c(X, Y))
  }

  # Simple circular functions
  deg2rad <- function(deg) return(deg*pi/180)
  rad2deg <- function(rad) return(180*rad/pi)


  # Find new destination, given velocity (km/yr), angle (º), time step (yr) and initial coordinates (ºlon, ºlat); max allowed jump is 1 deg
  # vell = vel[fcells] %>% pull(1); angg = ang[fcells] %>% pull(1); timestep = tstep; ll = llold
  destcoords <- function(vell, angg, timestep, ll, y_res, x_res) {
    latshift <- (abs(vell) * timestep * cos(deg2rad(angg))) / 111.325 # Calculate shift in lat
    latnew <- ll[,2] + latshift # Find new lat...first approximation
    lonshift <- (abs(vell) * timestep * sin(deg2rad(angg))) / (111.325 * cos(deg2rad(latnew))) # Shift in lon
    # Limit lonshift to at most 1 x resolution cell
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
             setNames(c("dlon", "dlat")))
  }

  # Function to find closest cooler/warmer cell within the union of two levels of 8-cell adjacency from the "departed" and "destination" cells
  get_dest_cell_fine <- function(rw) { # A to-from cell pair and the buffer around a cell on land that can be searched for a suitable cell
    # Find clumps of ocean; the rule is that you can't jump from one clump to another, because this would mean passing over unbroken land
    pos_depart <- data.frame(x = rw %>% pluck(3), y = rw %>% pluck(4))
    xy <- terra::xyFromCell(mn, as.vector(as.matrix(rw[1:2]))) %>% # Coordinates of cell centres for start cell and end cell %>%
      as.data.frame()
    bfr = ((y_res*111.325) + 1) * 1000 # Set buffer to be 1 grid-square width at the equator + 1 km
    xy <- xy %>%
      st_as_sf(coords = c("x", "y"), crs = "EPSG:4326")
    sp_buffer <- st_buffer(st_as_sf(xy), bfr) # Remembering that buffer is in metres
    buffer_zone <- terra::extract(mn, sp_buffer, cells = TRUE, xy = TRUE, touches = TRUE) %>%
      dplyr::select(-ID) %>%
      distinct() %>%
      dplyr::rename(sst = 1) %>% #*** rename "climatology", if needed
      dplyr::select(x, y, sst, cell) %>%
      drop_na(sst)
    clumped <- buffer_zone %>%
      dplyr::select(-cell) %>%
      rast(crs = "EPSG:4326") %>%
      terra::patches(directions = 8, allowGaps = FALSE)
    # Which clump did I start in?
    r1 <- rw %>%
      pluck(1) # We use this cell a lot, so let's just make it an object
    from_clump <- terra::extract(clumped, terra::xyFromCell(mn, r1)) %>%
      pluck(1)
    # What are the coordinates of cells within the search area that fall in the clump I came from?
    search_xy <- terra::xyFromCell(clumped, which(clumped[] == from_clump)) %>%
      as.data.frame()
    search_cells <- terra::cellFromXY(mn, search_xy) # Which cells are these
    # Get the ssts in the cells to search
    or <- terra::extract(mn, search_xy, cells = TRUE, xy = TRUE) %>%
      dplyr::rename(sst = 2, x = 4, y = 5)
    # If velocity is positive, find nearest cells that are cooler, if there are any
    if(unlist(vel[r1]) > 0) {
      # Find all cells in the search area that meet the sst criterion
      o <- or %>%
        dplyr::filter(sst < unlist(mn[r1])) %>%
        na.omit()
      if(nrow(o) == 0) {
        dest_cell <- NA # Set condition with which to ID stuck cells (i.e., if there are no suitable cells in the search area)
      } else {
        potential_dest_cells <- o %>%
          st_as_sf(coords = c("x", "y"), crs = crs(rast())) %>%
          mutate(distances = st_distance(st_as_sf(pos_depart, coords = c("x", "y"), crs = crs(rast())), .) %>%
                   as.vector()) %>%
          arrange(distances)
        n_cells <- 25 # Number of "close" cells to pick among *** Can change this, or set as an argument (just remove the hard-code here)
        if(nrow(potential_dest_cells) >= n_cells) {
          dest_cell <- potential_dest_cells %>%
            slice_head(n = n_cells) %>% # Get the closest cells
            slice_sample(n = 1) %>% # Pick one at random
            pull(cell)
        } else { # If fewer than n_cells, sample 10% of what you have
          dest_cell <- potential_dest_cells %>%
            slice_head(n = ceiling(nrow(potential_dest_cells)/10)) %>%
            slice_sample(n = 1) %>%
            pull(cell)
        }
      }
    } else { # Otherwise, find nearest cells that is warmer, if there are any
      # Find all cells in the search area that meet the sst criterion
      o <- or %>%
        dplyr::filter(sst < unlist(mn[r1])) %>%
        na.omit()
      if(nrow(o) == 0) {
        dest_cell <- NA # Set condition with which to ID stuck cells (i.e., if there are no suitable cells in the search area)
      } else {
        potential_dest_cells <- o %>%
          st_as_sf(coords = c("x", "y"), crs = crs(rast())) %>%
          mutate(distances = st_distance(st_as_sf(pos_depart, coords = c("x", "y"), crs = crs(rast())), .) %>%
                   as.vector()) %>%
          arrange(distances)
        n_cells <- 25 # Number of "close" cells to pick among *** Can change this, or set as an argument (just remove the hard-code here)
        if(nrow(potential_dest_cells) >= n_cells) {
          dest_cell <- potential_dest_cells %>%
            slice_head(n = n_cells) %>% # Get the closest cells
            slice_sample(n = 1) %>% # Pick one at random
            pull(cell)
        } else { # If fewer than n_cells, sample 10% of what you have
          dest_cell <- potential_dest_cells %>%
            slice_head(n = ceiling(nrow(potential_dest_cells)/10)) %>%
            slice_sample(n = 1) %>%
            pull(cell)
        }
      }
    }
    return(dest_cell)
  }

  # Function for to find closest cooler/warmer cell within the union of two levels of 8-cell adjacency from the "departed" and "destination" cells
  # Use when working with the coarse grid!
  get_dest_cell_coarse <- function(rw, bfr = 75) { # A from-to cell pair and coordinates of departure point; bfr is a search radius in km for nearest cooler cell, if you end on land
    # Find clumps of ocean; the rule is that you can't jump from one clump to another, because this would mean passing over unbroken land
    pos_depart <- data.frame(x = rw %>% pluck(3), y = rw %>% pluck(4)) %>% # Departure coordinates
      st_as_sf(coords = c("x", "y"), crs = "EPSG:4326")
    # xy <- terra::xyFromCell(mn, as.vector(as.matrix(rw))) %>% # Coordinates of cell centres for start cell and end cell
    #     st_as_sf(coords = c("x", "y"), crs = "EPSG:4326")
    sp_buffer <- st_buffer(pos_depart, bfr*1000) # Buffer around departure point, remembering that buffer is in metres
    buffer_zone <- terra::extract(mn, sp_buffer, cells = TRUE, xy = TRUE) %>%
      dplyr::select(-ID) %>%
      distinct() %>%
      dplyr::rename(sst = 1) %>% #*** rename "climatology", if needed
      dplyr::select(x, y, sst, cell)
    clumped <- buffer_zone %>%
      dplyr::select(-cell) %>%
      rast(crs = "EPSG:4326") %>%
      terra::patches(directions = 8, allowGaps = FALSE)
    # Which clump did I start in?
    r1 <- rw %>%
      pluck(1) # We use this cell a lot, so let's just make it an object
    from_clump <- terra::extract(clumped, terra::xyFromCell(mn, r1)) %>%
      unlist()
    # What are the coordinates of cells within the search area that fall in the clump I came from?
    search_xy <- terra::xyFromCell(clumped, which(clumped[] == from_clump)) %>%
      as.data.frame()
    # Get the ssts in the cells to search
    or <- terra::extract(mn, search_xy, cells = TRUE, xy = TRUE) %>%
      dplyr::rename(sst = 2, x = 4, y = 5)
    # If velocity is positive, find nearest cooler cell, if there is one
    if(vel[r1] > 0) {
      o <- or %>%
        dplyr::filter(sst < unlist(mn[r1])) %>%
        na.omit()
      # %>%
      #   bind_cols(., terra::xyFromCell(mn, .$cell)) #*** Do we need this line?
      if(nrow(o) == 0) {
        dest_cell <- NA
      } else {
        # pt <- terra::xyFromCell(mn, r1) %>%
        #   data.frame()
        dest_cell <- st_distance(pos_depart,
                                 st_as_sf(data.frame(o), coords = c("x", "y"), crs = "EPSG:4326")) %>%
          which.min() %>%
          o$cell[.] # The closest cell with appropriate sst
      }
    } else { # Otherwise, find nearest warmer cell, if there is one
      o <- or %>%
        dplyr::filter(sst < unlist(mn[r1])) %>%
        na.omit()
      # %>%
      #   bind_cols(., terra::xyFromCell(mn, .$sst)) #*** Do we need this line?
      if(nrow(o) == 0) {
        dest_cell <- NA
      } else {
        # pt <- terra::xyFromCell(mn, r1) %>%
        #   data.frame()
        dest_cell <- st_distance(pos_depart,
                                 st_as_sf(data.frame(o), coords = c("x", "y"), crs = "EPSG:4326")) %>%
          which.min() %>%  # The closest cell with appropriate sst
          o$cell[.]
      }
    }
    return(dest_cell)
  }


  # Setup -------------------------------------------------------------------

  # A base raster with original resolution and extent
  r_base <- rast(res = c(x_res, y_res)) %>%
    crop(vel)
  # Constrain max velocity to avoid stepping over grid squares
  max_vel = 111.325*x_res/tstep
  vel[vel > max_vel] <- max_vel
  vel[vel < -max_vel] <- -max_vel
  # Sort out start points
  lonlat <- lonlat %>%
    dplyr::select(x, y) %>%  # Collect just lon and lat (in case there's anything else there)
    as.data.frame()
  # Get initial descriptors
  tcells <- terra::cellFromXY(vel, lonlat) # # Cell IDs of starting cells
  n <- nrow(lonlat) # Get number of cells in your sequence
  # Set up variables to catch results
  llon <- llat <- cellIDend <- Steps <- cellIDs <- trajIDs <-  list() # Initiate lists to catch results
  # Populate the first slots with starting points
  cellIDs[[1]] <- tcells
  trajIDs[[1]] <- 1:n
  llon[[1]] <- lonlat[,1]
  llat[[1]] <- lonlat[,2]
  cellIDend[[1]] <- tcells
  Steps[[1]] <- rep(0, n)
  # Set up objects that keep track of things
  sflags <- rep(NA, n)
  trj_id <- trajIDs[[1]]


  # Loop --------------------------------------------------------------------

  # Loop through the trajectories
  for(i in 1:(tyr/tstep)) { # 1:(tyr/tstep)
    #   for(i in x) { # 1:(tyr/tstep)
    llold <- lonlat # Take a copy of lonlat
    vell <- terra::extract(vel, llold) %>%
      pull(voccMag)
    angg <- terra::extract(ang, llold) %>%
      pull(voccAng)
    fcells <- terra::cellFromXY(vel,llold) # Get the cells that the trajectories start in
    # Get new locations
    lonlat <- destcoords(vell, angg, tstep, llold, y_res, x_res) # Extract lon and lat of landing point
    tcells <- terra::cellFromXY(vel,lonlat) # Get the cells that the trajectories end in
    sflags[which(is.na(tcells))] <- 1 # Sets the trajectory to "stuck" if it exits the velocity field (i.e., tcells == NA)
    # Remove out-of-bounds cells
    in_bounds <- which(is.na(sflags))
    llold <- llold[in_bounds,]
    fcells <- fcells[in_bounds]
    lonlat <- lonlat[in_bounds,]
    tcells <- tcells[in_bounds]
    sflags <- sflags[in_bounds]
    trj_id <- trj_id[in_bounds]
    # Deal with cells that end on land
    onland <- which(is.na(vel[tcells]))
    if(length(onland) > 0){ # Only bother if there is at least one cell that returns a NA velocity = land
      if(grid_resolution == "fine") { #*** Fine switch
        # Collect the stuff we need here for cells that are onland
        fpos <- llold[onland,]
        tpos <- lonlat[onland,]
        SFlags <- sflags[onland]
        fcell <- fcells[onland]
        tcell <- tcells[onland]
        ft <- data.frame(fcell = fcell, tcell = tcell, fx = fpos %>% pluck(1), fy = fpos %>% pluck(2), code = paste(fcell, tcell, sep = " "), ref = 1:length(fcell))
        # ft <- data.frame(fcell = fcell, tcell = tcell, code = paste(fcell, tcell, sep = " "), ref = 1:length(fcell))
        ttcell <- apply(ft[,1:4], 1, get_dest_cell_fine) #*** fine switch
        # Filter "stuck" flags here
        stuck <- which(is.na(ttcell)) # This is set in get_dest_cell(), where no cell in the "catchment" has a suitable sst to facilitate movement
        unstuck <- which(!is.na(ttcell))
        SFlags[stuck] <- 1 # Adjust flags
        # Make data frame to catch data needed to find new positions
        ttpos <- data.frame(x = rep(NA, length(onland)), y = rep(NA, length(onland)))
        ttpos[stuck,] <- fpos[stuck,] # If they're stuck, pass on starting points
        ttcell[stuck] <- fcell[stuck] # If they're stuck, pass on starting cells

        ttpos[unstuck,] <- xyFromCell(mn, ttcell[unstuck])
        # Collect results
        lonlat[onland,] <- ttpos
        tcells[onland] <- ttcell
        sflags[onland] <- SFlags
      } else {
        # Old onland loop here
        # Collect the stuff we need here
        fpos <- llold[onland,]
        tpos <- lonlat[onland,]
        SFlags <- sflags[onland]
        fcell <- fcells[onland]
        tcell <- tcells[onland]
        # ft <- data.frame(fcell = fcell, tcell = tcell, code = paste(fcell, tcell, sep = " "), ref = 1:length(fcell))
        # ttcell <- apply(ft[,1:2], 1, get_dest_cell)
        ft <- data.frame(fcell = fcell, tcell = tcell, fx = fpos %>% pluck(1), fy = fpos %>% pluck(2), code = paste(fcell, tcell, sep = " "), ref = 1:length(fcell))
        ttcell <- apply(ft[,1:4], 1, get_dest_cell_coarse)
        # Filter "stuck" flags here
        stuck <- which(is.na(ttcell)) # This is set in get_dest_cell(), where no cell in the "catchment" has a suitable sst to facilitate movement
        unstuck <- which(!is.na(ttcell))
        SFlags[stuck] <- 1 # Adjust flags
        # Make data frame to catch data needed to find new positions #***done up to here
        ttpos <- data.frame(x = rep(NA, length(onland)), y = rep(NA, length(onland)))
        ttpos[stuck,] <- fpos[stuck,] # If they're stuck, pass on starting points
        ttcell[stuck] <- fcell[stuck] # If they're stuck, pass on starting cells
        if(length(unstuck) > 0) {
          tt_original_cells <- terra::cellFromXY(mn, fpos[unstuck,]) # Departure cells in the resolution of the original velocity field
          ttdat <- tibble(ttcell = ttcell[unstuck]) %>% # Destination cells (nearest cell with appropriate sst)
            mutate(loncent = terra::xFromCell(mn, ttcell), # Coordinates of destination cell
                   latcent = terra::yFromCell(mn, ttcell)) %>%  # Coordinates of destination cell
            mutate(e = NA, w = NA, n = NA, s = NA, dlon = NA, dlat = NA) # To facilitate finding corners of the cells
          corner_block_size <- 0.25*x_res # The "corner" is set to 0.25 of the grid square at the original resolution
          # Send trajectory to the nearest corner of the appropriate cell, where corner is a quarter of the grid size. Position is "fuzzed" within this corner at random.
          ttdat$e <- ttdat$loncent + (0.5 * x_res) - (runif(1, -1, 1) * corner_block_size) # The centre of the "corner" +- random fuzz...
          ttdat$w <- ttdat$loncent - (0.5 * x_res) + (runif(1, -1, 1) * corner_block_size) # The centre of the "corner" +- random fuzz...
          ttdat$n <- ttdat$latcent + (0.5 * y_res) - (runif(1, -1, 1) * corner_block_size) # The centre of the "corner" +- random fuzz...
          ttdat$s <- ttdat$latcent - (0.5 * y_res) + (runif(1, -1, 1) * corner_block_size) # The centre of the "corner" +- random fuzz...
          coords <- with(ttdat, cbind(n, e, n, w, s, w, s, e)) # NE, NW, SW, SE corners' coordinates
          # Find distances from departure point to corners of ttcell
          get_dist <- function(y1, x1, y2, x2) {
            pt1 <- st_as_sf(data.frame(x = x1, y = y1), coords = c("x", "y"), crs = crs(rast()))
            pt2 <- st_as_sf(data.frame(x = x2, y = y2), coords = c("x", "y"), crs = crs(rast()))
            out <- st_distance(pt1, pt2, by_element = TRUE) %>%
              as.vector()
            return(out)
          }
          corners <- data.frame(ne = get_dist(fpos[unstuck,2], fpos[unstuck,1], coords[,1], coords[,2])) %>%
            mutate(nw = get_dist(fpos[unstuck,2], fpos[unstuck,1], coords[,3], coords[,4]),
                   sw = get_dist(fpos[unstuck,2], fpos[unstuck,1], coords[,5], coords[,6]),
                   se = get_dist(fpos[unstuck,2], fpos[unstuck,1], coords[,7], coords[,8]))
          # Select nearest corner
          cornset <- apply(corners, 1, mfind)*2 # Identify which corners for each onland cell. Have to mul by 2 to shift along correctly.
          cornset <- cbind(cornset, coords) # Add in coordinates
          ttdat[,8:9] <- data.frame(t(apply(cornset, 1, mplace))) # Add in coordinates of correct corner point
          ttpos[unstuck,] <- xyFromCell(mn, ttcell[unstuck])
          ttcell[unstuck] <- terra::cellFromXY(mn, ttpos[unstuck,])
        }
        # Collect results
        lonlat[onland,] <- ttpos
        tcells[onland] <- ttcell
        sflags[onland] <- SFlags
      }
    }
    # Pass on only those cells that are not stuck
    # if(sum(is.na(lonlat[,1])) > 0) {break}
    cells_to_keep <- which(is.na(sflags))
    lonlat <- lonlat[cells_to_keep,]
    # Collect results
    llon[[i + 1]] <- lonlat[,1] # Add lon to the list
    llat[[i + 1]] <- lonlat[,2]# Add lat to the list
    cellIDs[[i + 1]] <- fcells[cells_to_keep] # Add departure cellIDs to the list
    cellIDend[[i + 1]] <- tcells[cells_to_keep] # Add arrival cellIDs to the list
    Steps[[i + 1]] <- rep(i, length(tcells[cells_to_keep])) # Capture the time_step
    trajIDs[[i + 1]] <- trj_id[cells_to_keep]
    sflags <- sflags[cells_to_keep] # Adjust the flags list
    trj_id <- trj_id[cells_to_keep] # Adjust the trajectory IDs
    print(c(i, length(onland), round(100*i/(tyr/tstep), 3), date())) # Keep a reference of progress
  }
  return(cbind(Steps = Steps %>%
                 reduce(c),
               lon = llon %>%
                 reduce(c),
               lat = llat %>%
                 reduce(c),
               ID = trajIDs %>%
                 reduce(c),
               start_cell = cellIDs %>%
                 reduce(c),
               end_cell = cellIDend %>%
                 reduce(c)) %>%
           as_tibble())
}

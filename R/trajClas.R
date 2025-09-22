#' Climate velocity trajectory classification
#'
#' Function for the spatial classification of cells based on VoCC trajectories after Burrows et al. (2014). The function performs
#' a hierarchical sequential classification based on length of trajectories, geographical features, and the relative abundance of
#' trajectories ending in, starting from and flowing through each cell. Essentially, cells are first classified as non-moving,
#' slow-moving and fast-moving relative to the distance a trajectory will cover over the projection period based on local climate velocities.
#' Two types of climate sinks are then identified among the fast-moving cells: (i) boundary (e.g., coastal) cells disconnected from cooler (warmer)
#' neighbouring cells under a locally warming (cooling) climate, and (ii) locations with endorheic spatial gradients where the velocity angles of
#' neighbouring cells converge towards their central point of intersection. Finally, the remaining cells are classified by reference to the total
#' number of trajectories per cell based on the proportions of the number of trajectories starting from (Nst), ending in (Nend), and flowing
#' through (NFT) a cell over the period. Based on these proportions, cells are classified into five classes: (1) climate sources, when no
#' trajectories end in a cell (Nend = 0); (2) relative climate sinks, when the relative number of trajectories ending in a cell is high and the
#' proportion of starting trajectories is low; (3) corridors as cells with a high proportion of trajectories passing through; and (4) divergence
#' and (5) convergence cells identified from the remaining cells as those where fewer/more trajectories ended than started in that
#' cell, respectively.
#'
#' @param traj \code{data.frame} as retuned by voccTraj containing the coordinates
#' and identification number for each trajectory.
#' @param vel \code{SpatRaster} with the magnitude of gradient-based climate velocity.
#' @param ang \code{SpatRaster} with velocity angles.
#' @param mn \code{SpatRaster} with mean climatic values for the study period.
#' @param trajSt \code{integer} number of trajectories starting from each cell or spatial unit.
#' @param tyr \code{integer} number of years comprising the projected period.
#' @param nmL \code{numeric} upper threshold (distance units as per vel object) up to which
#' a trajectory is considered to have traveled a negligible distance over the study period (non-moving).
#' @param smL \code{numeric} upper threshold up to which a trajectory is considered to have traveled a small
#' distance over the study period (slow-moving).
#' @param Nend \code{numeric} the percentage of trajectories ending to be used as threshold in the classification.
#' @param Nst \code{numeric} the percentage of trajectories starting to be used as threshold in the classification.
#' @param NFT \code{numeric} the percentage of trajectories flowing through to be used as threshold in the classification.
#' @param DateLine \code{logical} does the raster extent cross the international date line? (default "FALSE").
#'
#' @return A \code{SpatRaster} containing the trajectory classification ("TrajClas"),
#' as well as those based on trajectory length ("ClassL"; 1 non-moving, 2 slow-moving, 3 fast-moving cells),
#' boundrary ("BounS") and internal sinks ("IntS"), and the proportion of trajectories ending("PropEnd"),
#' flowing through ("PropFT") and starting ("PropSt"). The trajectory classes ("TrajClas") are (1) non-moving,
#' (2) slow-moving, (3) internal sinks, (4) boundary sinks, (5) sources, (6) relative sinks, (7) corridors,
#' (8) divergence and (9) convergence.
#'
#' @references \href{https://www.nature.com/articles/nature12976}{Burrows et al. 2014}. Geographical limits to species-range shifts are suggested by climate velocity. Nature, 507, 492-495.
#'
#' @seealso{\code{\link{voccTraj}}}
#'
#' @export
#' @examples
#'
#' \dontrun{
#' EEZ <- terra::vect(system.file("extdata", "EEZ.gpkg", package = "VoCC"))
#'
#' HSST <- terra::rast(system.file("extdata", "HadiSST.tif", package = "VoCCdata")) %>%
#'   terra::crop(terra::ext(EEZ) + 10)
#'
#' # input raster layers
#' yrSST <- sumSeries(HSST,
#'   p = "1960-01/2009-12", yr0 = "1955-01-01", l = terra::nlyr(HSST),
#'   fun = function(x) colMeans(x, na.rm = TRUE), freqin = "months", freqout = "years"
#' )
#'
#' mn <- terra::mean(yrSST, na.rm = TRUE)
#' tr <- tempTrend(yrSST, th = 10)
#' sg <- spatGrad(yrSST, th = 0.0001, projected = FALSE)
#' v <- gVoCC(tr, sg)
#' vel <- v[[1]]
#' ang <- v[[2]]
#'
#' # Get the set of starting cells for the trajectories and calculate trajectories
#' lonlat <- stats::na.omit(data.frame(
#'   terra::xyFromCell(vel, 1:terra::ncell(vel)),
#'   terra::values(vel), terra::values(ang), terra::values(mn)
#' ))[, 1:2]
#'
#' traj <- voccTraj(lonlat, vel, ang, mn, tstep = 1/4, tyr = 20, seed = 23)
#'
#'
#' # Generate the trajectory-based classification
#' clas <- trajClas(traj, vel, ang, mn,
#'   trajSt = 16, tyr = 50, nmL = 20, smL = 100,
#'   Nend = 45, Nst = 15, NFT = 70, DateLine = FALSE
#' )
#'
#' # Classify raster / build attribute table
#' clasr <- terra::as.factor(clas[[7]])
#' rat_r <- data.frame(ID = sort(unique(terra::values(clas[[7]]))),
#'                     trajcat = c("N-M", "S-M", "IS", "BS", "Srce",
#'                                "RS", "Cor", "Div", "Con")[sort(unique(terra::values(clas[[7]])))])
#' levels(clasr) <- rat_r
#' terra::plot(clasr)
#' }
#'
trajClas <- function(traj, vel, ang, mn, trajSt, tyr, nmL, smL, Nend, Nst, NFT, DateLine = FALSE) {

  ang1 <- ang2 <- ang3 <- ang4 <- d1 <- d2 <- d3 <- d4 <- NULL # Fix devtools check warnings
  isink <- .SD <- .N <- cid <- coastal <- val <- ID <- NULL # Fix devtools check warnings

  # MEMORY LEAK FIX: Create template raster once and reuse
  template_raster <- terra::rast(ang)

  # Force loading of ang values to avoid lazy loading warnings
  invisible(terra::values(ang))

  TrajEnd <- terra::rast(template_raster)
  terra::values(TrajEnd) <- 0
  TrajFT <- terra::rast(template_raster)
  terra::values(TrajFT) <- 0
  TrajSt <- terra::rast(template_raster)
  terra::values(TrajSt) <- 0
  IntS <- terra::rast(template_raster)
  terra::values(IntS) <- 0
  BounS <- terra::rast(template_raster)
  terra::values(BounS) <- 0
  TrajClas <- terra::rast(template_raster)
  terra::values(TrajClas) <- 0

  # add cell ID to the data frame
  traj <- data.table::data.table(traj)
  traj$cid <- terra::cellFromXY(ang, traj[, 1:2])

  # A. Number of traj starting from each cell
  # Set trajSt for all non-NA cells in ang
  terra::values(TrajSt)[!is.na(terra::values(ang))] <- trajSt

  # B. Number of traj ending in each cell
  tr <- traj[, data.table::.SD[.N], by = ID] # subset last point of each trajectory
  enTraj <- tr[, .N, by = cid]
  TrajEnd <- terra::mask(TrajEnd, vel)
  terra::values(TrajEnd)[enTraj$cid] <- enTraj$N

  # C. Number of traj flowing through each cell
  cxtrj <- unique(traj, by = c("ID", "cid"))
  TotTraj <- cxtrj[, .N, by = cid] # total number of trajectories per cell
  TrajFT <- terra::mask(TrajFT, vel)
  terra::values(TrajFT)[TotTraj$cid] <- TotTraj$N
  TrajFT <- TrajFT - TrajEnd - TrajSt # subtract traj starting and ending to get those actually transversing the cell
  terra::values(TrajFT)[terra::values(TrajFT) < 0] <- 0 # to avoid negative values in ice covered cells

  # C. Identify cell location for internal sinks (groups of 4 endorheic cells with angles pointing inwards)
  ll <- data.table::data.table(terra::xyFromCell(ang, 1:terra::ncell(ang)))
  ll[, 1:2] <- ll[, 1:2] + 0.1 # add small offset to the centroid

  # Terra replacement for raster::fourCellsFromXY()
  # The original function returned 4 cells in a 2x2 arrangement around each point
  # We recreate this by getting cells at the four corners of a 2x2 grid

  coords <- as.matrix(ll[, 1:2])
  res_x <- terra::xres(ang)
  res_y <- terra::yres(ang)

  # Create the four corner positions for a 2x2 grid centered on each coordinate
  # Offsets for: bottom-left, bottom-right, top-left, top-right
  offsets <- matrix(c(-res_x/2, -res_y/2,  # bottom-left
                      res_x/2, -res_y/2,   # bottom-right
                      -res_x/2, res_y/2,   # top-left
                      res_x/2, res_y/2),   # top-right
                    ncol = 2, byrow = TRUE)

  # OPTIMIZED: Vectorized corner coordinate calculation
  # Get the four cells for each coordinate point
  n_coords <- nrow(coords)
  a <- matrix(as.numeric(NA), nrow = n_coords, ncol = 4)

  # Only process valid coordinates
  valid_idx <- !is.na(coords[,1]) & !is.na(coords[,2])
  valid_coords <- coords[valid_idx, , drop = FALSE]

  if (nrow(valid_coords) > 0) {
    # Vectorized corner calculation for all valid coordinates at once
    n_valid <- nrow(valid_coords)

    # Create all corner coordinates in one operation
    all_corners <- array(NA, dim = c(n_valid, 4, 2))
    for(k in 1:4) {
      # Fix dimension mismatch by ensuring proper matrix dimensions
      offset_matrix <- matrix(rep(offsets[k,], each = n_valid), nrow = n_valid, ncol = 2)
      all_corners[, k, ] <- offset_matrix + valid_coords
    }

    # Get cell numbers for all corners at once
    for(k in 1:4) {
      corner_matrix <- all_corners[, k, , drop = FALSE]
      dim(corner_matrix) <- c(n_valid, 2)
      cell_nums <- terra::cellFromXY(ang, corner_matrix)
      a[valid_idx, k] <- as.numeric(cell_nums)
    }
  }

  # OPTIMIZED: Vectorized row sorting using apply
  a <- t(apply(a, 1, sort, na.last = TRUE))

  # If date line crossing, correct sequences on date line
  if (DateLine == TRUE) {
    a[seq(terra::ncol(ang), by = terra::ncol(ang), length = terra::nrow(ang)), ] <- t(apply(a[seq(terra::ncol(ang), by = terra::ncol(ang), length = terra::nrow(ang)), ], 1, function(x) {
      x[c(2, 1, 4, 3)]
    }))
  }

  # OPTIMIZED: Vectorized angle extraction
  # Extract the angles for each group of 4 cells
  b <- matrix(NA, nrow = nrow(a), ncol = 4)

  # Get all angle values once to avoid repeated calls
  ang_values <- terra::values(ang)
  max_cell <- length(ang_values)

  # Vectorized angle extraction
  valid_cells <- !is.na(a) & a > 0 & a <= max_cell
  b[valid_cells] <- ang_values[a[valid_cells]]
  # Ensure a and b have the same number of rows before combining
  if(nrow(a) != nrow(b)) {
    # Pad b to match a if needed
    if(nrow(b) < nrow(a)) {
      b_padded <- matrix(NA, nrow = nrow(a), ncol = 4)
      b_padded[1:nrow(b), ] <- b
      b <- b_padded
    }
  }
  ll[, c("c1", "c2", "c3", "c4", "ang1", "ang2", "ang3", "ang4") := data.frame(a, b)]

  # now look if the 4 angles point inwards (internal sink)
  ll[, c("d1", "d2", "d3", "d4") := list(((ang1 - 180) * (90 - ang1)), ((ang2 - 270) * (180 - ang2)), ((ang3 - 90) * (0 - ang3)), ((ang4 - 360) * (270 - ang4)))]
  ll[, isink := 0L]
  ll[d1 > 0 & d2 > 0 & d3 > 0 & d4 > 0, isink := 1L]

  # get the cids for the cells contained in the sinks
  celdas <- ll[isink == 1, 3:6]
  IntS <- terra::mask(IntS, vel)

  # Convert data.table columns to vectors and combine
  if(nrow(celdas) > 0) {
    # Convert to numeric vectors to avoid list issues
    sink_cells <- c(as.numeric(celdas$c1), as.numeric(celdas$c2),
                    as.numeric(celdas$c3), as.numeric(celdas$c4))
    sink_cells <- sink_cells[!is.na(sink_cells)]  # Remove NA values
    if(length(sink_cells) > 0) {
      terra::values(IntS)[sink_cells] <- 1
    }
  }

  # D. Identify cell location for boundary sinks (coastal cells which are disconected from cooler climates under warming or warmer climates under cooling)
  # detect coastal cells
  coast <- suppressWarnings(terra::boundaries(vel, inner = TRUE)) # terra uses 'inner' parameter instead of 'type'

  # OPTIMIZED: Pre-extract raster values once to avoid repeated terra::values() calls
  vel_values <- terra::values(vel)
  mn_values <- terra::values(mn)
  coast_values <- terra::values(coast)

  # make a list of vel values and SST values for each coastal cells and their marine neighbours
  cc <- stats::na.omit(data.table::data.table(cid = 1:terra::ncell(vel), coast = coast_values))
  ad <- terra::adjacent(vel, cc$cid, directions = 8, include = TRUE) # matrix with adjacent cells
  # Sort the adjacency matrix to ensure consistent ordering (replaces sorted=TRUE from raster package)
  ad <- ad[order(ad[, 1], ad[, 2]), ]
  ad <- data.table::data.table(
    coastal = ad[, 1],
    adjacent = ad[, 2],
    cvel = vel_values[ad[, 1]],
    ctemp = mn_values[ad[, 1]],
    atemp = mn_values[ad[, 2]]
  )

  # locate the sinks
  ad <- stats::na.omit(ad[ad$cvel != 0, ]) # remove cells with 0 velocity (ice) and with NA (land neighbours)
  j <- ad[, ifelse(.SD$cvel > 0, all(.SD$ctemp <= .SD$atemp), all(.SD$ctemp >= .SD$atemp)), by = coastal]
  data.table::setkey(j)
  j <- unique(j)
  BounS <- terra::mask(BounS, vel)
  boundary_cells <- unique(subset(j$coastal, j$V == TRUE))
  if(length(boundary_cells) > 0) {
    terra::values(BounS)[as.numeric(boundary_cells)] <- 1
  }

  # OPTIMIZED: Pre-extract ang values and minimize raster operations
  ang_values <- terra::values(ang)

  # Total number of trajectories per cell and proportions per cell
  TrajTotal <- TrajSt + TrajFT + TrajEnd
  terra::values(TrajTotal)[is.na(ang_values)] <- NA
  PropTrajEnd <- (TrajEnd / TrajTotal) * 100
  PropTrajFT <- (TrajFT / TrajTotal) * 100
  PropTrajSt <- (TrajSt / TrajTotal) * 100

  # OPTIMIZED: Direct classification using pre-extracted values - avoid creating intermediate rasters
  abs_vel_values <- abs(vel_values)

  # Vectorized classification instead of terra::classify
  ClassMov_values <- ifelse(abs_vel_values <= (nmL / tyr), 1,
                           ifelse(abs_vel_values <= (smL / tyr), 2, 3))
  ClassMov_values[is.na(vel_values)] <- NA

  # OPTIMIZED: Pre-extract all raster values for classification
  IntS_values <- terra::values(IntS)
  BounS_values <- terra::values(BounS)
  PropTrajEnd_values <- terra::values(PropTrajEnd)
  PropTrajSt_values <- terra::values(PropTrajSt)
  PropTrajFT_values <- terra::values(PropTrajFT)

  # Classify the cells using vectorized operations
  TrajClas_values <- rep(0, terra::ncell(TrajClas))
  TrajClas_values[is.na(vel_values)] <- NA

  # capture non-moving (1)
  TrajClas_values[ClassMov_values == 1] <- 1

  # capture slow-moving (2)
  TrajClas_values[ClassMov_values == 2] <- 2

  # capture internal (3) and (4) boundary sinks
  TrajClas_values[IntS_values == 1] <- 3
  TrajClas_values[BounS_values == 1] <- 4

  # OPTIMIZED: Vectorized classification for remaining cells
  remaining_cells <- which(TrajClas_values == 0)
  if (length(remaining_cells) > 0) {
    Nend_vals <- PropTrajEnd_values[remaining_cells]
    Nst_vals <- PropTrajSt_values[remaining_cells]
    NFT_vals <- PropTrajFT_values[remaining_cells]

    # Vectorized classification logic
    clas_vals <- ifelse(Nend_vals == 0, 5,
                       ifelse(Nend_vals > Nend & Nst_vals < Nst, 6,
                             ifelse(NFT_vals > NFT, 7,
                                   ifelse(Nend_vals < Nst_vals, 8, 9))))
    TrajClas_values[remaining_cells] <- clas_vals
  }

  # Set final values
  terra::values(TrajClas) <- TrajClas_values

  # Create ClassMov raster from the values before cleanup
  ClassMov <- terra::rast(template_raster)
  terra::values(ClassMov) <- ClassMov_values

  # MEMORY LEAK FIX: Clean up large intermediate objects before final assembly
  rm(TrajClas_values, ClassMov_values, IntS_values, BounS_values,
     PropTrajEnd_values, PropTrajSt_values, PropTrajFT_values,
     vel_values, mn_values, coast_values, ang_values, template_raster,
     abs_vel_values)

  # return raster
  s <- c(PropTrajEnd, PropTrajFT, PropTrajSt, ClassMov, IntS, BounS, TrajClas)
  names(s) <- c("PropEnd", "PropFT", "PropSt", "ClassL", "IntS", "BounS", "TrajClas")

  # Force garbage collection
  gc()

  return(s)
}

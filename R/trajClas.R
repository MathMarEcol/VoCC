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
#' @author Jorge Garcia Molinos
#' @examples
#' \dontrun{
#' HSST <- VoCC_get_data("HSST.tif")
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
#' # at 1/4-deg resolution (16 trajectories per 1-deg cell)
#' mnd <- terra::disagg(mn, 4)
#' veld <- terra::disagg(vel, 4)
#' angd <- terra::disagg(ang, 4)
#' lonlat <- stats::na.omit(data.frame(
#'   terra::xyFromCell(veld, 1:terra::ncell(veld)),
#'   terra::values(veld), terra::values(angd), terra::values(mnd)
#' ))[, 1:2]
#'
#' traj <- voccTraj(lonlat, vel, ang, mn, tyr = 50, correct = TRUE)
#'
#' # Generate the trajectory-based classification
#' clas <- trajClas(traj, vel, ang, mn,
#'   trajSt = 16, tyr = 50, nmL = 20, smL = 100,
#'   Nend = 45, Nst = 15, NFT = 70, DateLine = FALSE
#' )
#'
#' # Define first the colour palette for the full set of categories
#' my_col <- c(
#'   "gainsboro", "darkseagreen1", "coral4", "firebrick2", "mediumblue", "darkorange1",
#'   "magenta1", "cadetblue1", "yellow1"
#' )
#' # Keep only the categories present in our raster
#' my_col <- my_col[sort(unique(terra::values(clas[[7]])))]
#'
#' # Classify raster / build attribute table
#' clasr <- terra::as.factor(clas[[7]])
#' rat_r <- data.frame(ID = sort(unique(terra::values(clas[[7]]))),
#'                     trajcat = c("N-M", "S-M", "IS", "BS", "Srce",
#'                                "RS", "Cor", "Div", "Con")[sort(unique(terra::values(clas[[7]])))])
#' terra::cats(clasr) <- rat_r
#' # Produce the plot using the rasterVis levelplot function
#' rasterVis::levelplot(clasr,
#'   col.regions = my_col,
#'   xlab = NULL, ylab = NULL, scales = list(draw = FALSE)
#' )
#' }
trajClas <- function(traj, vel, ang, mn, trajSt, tyr, nmL, smL, Nend, Nst, NFT, DateLine = FALSE) {

  ang1 <- ang2 <- ang3 <- ang4 <- d1 <- d2 <- d3 <- d4 <- NULL # Fix devtools check warnings
  isink <- .SD <- .N <- cid <- coastal <- val <- trajIDs <- NULL # Fix devtools check warnings

  # Force loading of ang values to avoid lazy loading warnings
  invisible(terra::values(ang))

  TrajEnd <- TrajFT <- TrajSt <- IntS <- BounS <- TrajClas <- terra::rast(ang)

  # add cell ID to the data frame
  traj <- data.table::data.table(traj)
  traj$cid <- terra::cellFromXY(ang, traj[, 1:2])

  # A. Number of traj starting from each cell
  # Set trajSt for all non-NA cells in ang
  terra::values(TrajSt)[!is.na(terra::values(ang))] <- trajSt

  # B. Number of traj ending in each cell
  tr <- traj[, data.table::.SD[.N], by = trajIDs] # subset last point of each trajectory
  enTraj <- tr[, .N, by = cid]
  terra::values(TrajEnd) <- 0
  TrajEnd <- terra::mask(TrajEnd, vel)
  terra::values(TrajEnd)[enTraj$cid] <- enTraj$N

  # C. Number of traj flowing through each cell
  cxtrj <- unique(traj, by = c("trajIDs", "cid"))
  TotTraj <- cxtrj[, .N, by = cid] # total number of trajectories per cell
  terra::values(TrajFT) <- 0
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

  # Get the four cells for each coordinate point
  # Initialize as numeric matrix explicitly
  a <- matrix(as.numeric(NA), nrow = nrow(coords), ncol = 4)

  for(i in 1:nrow(coords)) {
    if(!is.na(coords[i,1]) && !is.na(coords[i,2])) {
      # Calculate corner coordinates
      corner_coords <- sweep(offsets, 2, coords[i,], "+")
      # Get cell numbers for each corner - handle one at a time
      for(k in 1:4) {
        cell_num <- terra::cellFromXY(ang, corner_coords[k, , drop = FALSE])
        a[i, k] <- as.numeric(cell_num[1])
      }
    }
  }

  # Sort each row as in original code
  for(i in 1:nrow(a)) {
    a[i,] <- sort(a[i,])
  }

  # If date line crossing, correct sequences on date line
  if (DateLine == TRUE) {
    a[seq(terra::ncol(ang), by = terra::ncol(ang), length = terra::nrow(ang)), ] <- t(apply(a[seq(terra::ncol(ang), by = terra::ncol(ang), length = terra::nrow(ang)), ], 1, function(x) {
      x[c(2, 1, 4, 3)]
    }))
  }

  # Extract the angles for each group of 4 cells
  # Use direct indexing approach that works better with terra
  b <- matrix(NA, nrow = nrow(a), ncol = 4)

  # Get all angle values once to avoid repeated calls
  ang_values <- terra::values(ang)

  for(i in 1:nrow(a)) {
    for(j in 1:4) {
      cell_val <- a[i,j]  # Should now be numeric
      if(!is.na(cell_val) && cell_val > 0 && cell_val <= length(ang_values)) {
        b[i,j] <- ang_values[cell_val]
      }
    }
  }
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
  terra::values(IntS) <- 0
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

  # make a list of vel values and SST values for each coastal cells and their marine neighbours
  cc <- stats::na.omit(data.table::data.table(cid = 1:terra::ncell(vel), coast = terra::values(coast)))
  ad <- terra::adjacent(vel, cc$cid, directions = 8, include = TRUE) # matrix with adjacent cells
  # Sort the adjacency matrix to ensure consistent ordering (replaces sorted=TRUE from raster package)
  ad <- ad[order(ad[, 1], ad[, 2]), ]
  ad <- data.table::data.table(
    coastal = ad[, 1],
    adjacent = ad[, 2],
    cvel = terra::values(vel)[ad[, 1]],
    ctemp = terra::values(mn)[ad[, 1]],
    atemp = terra::values(mn)[ad[, 2]]
  )

  # locate the sinks
  ad <- stats::na.omit(ad[ad$cvel != 0, ]) # remove cells with 0 velocity (ice) and with NA (land neighbours)
  j <- ad[, ifelse(.SD$cvel > 0, all(.SD$ctemp <= .SD$atemp), all(.SD$ctemp >= .SD$atemp)), by = coastal]
  data.table::setkey(j)
  j <- unique(j)
  terra::values(BounS) <- 0
  BounS <- terra::mask(BounS, vel)
  boundary_cells <- unique(subset(j$coastal, j$V == TRUE))
  if(length(boundary_cells) > 0) {
    terra::values(BounS)[as.numeric(boundary_cells)] <- 1
  }

  # Total number of trajectories per cell and proportions per cell
  TrajTotal <- TrajSt + TrajFT + TrajEnd
  terra::values(TrajTotal)[is.na(terra::values(ang))] <- NA
  PropTrajEnd <- (TrajEnd / TrajTotal) * 100
  PropTrajFT <- (TrajFT / TrajTotal) * 100
  PropTrajSt <- (TrajSt / TrajTotal) * 100

  # reclassify by traj length
  rclM <- matrix(c(0, (nmL / tyr), 1, (nmL / tyr), (smL / tyr), 2, (smL / tyr), Inf, 3), ncol = 3, byrow = TRUE)
  v <- terra::rast(vel)
  terra::values(v) <- abs(terra::values(vel))
  ClassMov <- terra::classify(v, rclM)

  # Classify the cells
  terra::values(TrajClas) <- 0
  TrajClas <- terra::mask(TrajClas, vel)

  # capture non-moving (1)
  terra::values(TrajClas)[terra::values(ClassMov) == 1] <- 1

  # capture slow-moving (2)
  terra::values(TrajClas)[terra::values(ClassMov) == 2] <- 2

  # capture internal (3) and (4) boundary sinks
  terra::values(TrajClas)[terra::values(IntS) == 1] <- 3
  terra::values(TrajClas)[terra::values(BounS) == 1] <- 4

  # Classify remaining cells into sources(5), rel sinks(6), corridors(7), divergence(8) and convergence(9)
  d <- stats::na.omit(data.table::data.table(cid = 1:terra::ncell(TrajClas), val = terra::values(TrajClas)))
  d <- d[val == 0, 1]
  d[, Nend := terra::values(PropTrajEnd)[d$cid]]
  d[, Nst := terra::values(PropTrajSt)[d$cid]]
  d[, NFT := terra::values(PropTrajFT)[d$cid]]
  d$clas <- ifelse(d$Nend == 0, 5, ifelse(d$Nend > Nend & d$Nst < Nst, 6, ifelse(d$NFT > NFT, 7, ifelse(d$Nend < d$Nst, 8, 9))))
  terra::values(TrajClas)[d$cid] <- d$clas

  # return raster
  s <- c(PropTrajEnd, PropTrajFT, PropTrajSt, ClassMov, IntS, BounS, TrajClas)
  names(s) <- c("PropEnd", "PropFT", "PropSt", "ClassL", "IntS", "BounS", "TrajClas")
  return(s)
}

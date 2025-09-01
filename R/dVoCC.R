#' Distance-based velocity based on geographically closest climate analogue
#'
#' Function to calculate the geographically closest climate analogues and related distance-based velocity. Cell analogues
#' are identified by comparing the baseline climatic conditions at each focal cell with those existing for all
#' other (target) cells in the future by reference to a specified climatic threshold. The function allows for the
#' specification of search distances and incorporates both least-cost path and Great Circle (as-the-crow-flies) distances.
#'
#' @importFrom foreach %dopar%
#'
#' @param clim \code{data.frame} with the value for the climatic parameters (columns) by cell (rows), arranged as follows (see examples below):
#' The first 2n columns must contain the present and future values for each of the n climatic variables (V1p, V1f, V2p, V2f,...).
#' Where cell-specific analogue thresholds (see "variable" in "method" below) are to be calculated, the next (2n+1:3n) columns
#' should contain the standard deviation (or any other measure of climatic variability) of each variable for the baseline period.
#' These columns are not required if using the "Single" method. The last three columns of the table should contain an identifyier and centroid coordinates of each cell.
#' @param n \code{integer} defining the number of climatic variables.
#' @param tdiff \code{integer} defining the number of years (or other temporal unit) between periods.
#' @param method \code{character string} specifying the analogue method to be used. 'Single': a constant, single analogue threshold
#' for each climate variable is applied to all cells (Ohlemuller et al. 2006, Hamann et al. 2015); climate analogy corresponds to target cells
#' with values below the specified threshold for each climatic variable. 'Variable': a cell-specific climate threshold is used for each climatic variable
#' to determine the climate analogues associated with each cell by reference to its baseline climatic variability (Garcia Molinos et al. 2017).
#' @param climTol \code{numeric} a vector of length n giving the tolerance threshold defining climate analogue conditions for each climatic variable.
#' If a cell-specific threshold is being used, this function parameter should be passed as NA.
#' @param geoTol \code{integer} impose a geographical distance threshold (in km for lat/lon or map units if projected).
#' If used, the pool of potential climate analogues will be limited to cells within that distance from the focal cell.
#' @param distfun \code{character string} specifying the function to be used for estimating distances
#' between focal and target cells. Either 'Euclidean', 'GreatCircle' (Great Circle Distances).
#' 'LeastCost' (Least Cost Path Distances) is NOT CURRENTLY IMPLEMENTED. LeastCost requires a transition matrix supplied to the function via de 'trans' argument.
#' @param trans \code{TransitionLayer} NOT CURRENTLY IMPLEMENTED gdistance object to be used for the analogue search if distfun = 'LeastCost'.
#' @param lonlat \code{logical} is the analysis to be done in unprojected (lon/lat) coordinates?
#'
#' @return A \code{data.frame} containing the cell id of the future analogue for each focal cell (NA = no analogue available),
#' together with the climatic ("climDis") and geographical ("geoDis") distances in input units,
#' the bearing ("ang", degrees North), and resulting climate velocity ("vel", km/yr). Mean climatic distances are returned for multivariate analogues.
#'
#' @references \href{https://onlinelibrary.wiley.com/doi/full/10.1111/j.1466-822X.2006.00245.x}{Ohlemuller et al. 2006}. Towards European climate risk surfaces: the extent and distribution of analogous and non-analogous climates 1931-2100. Global Ecology and Biogeography, 15, 395-405. \cr
#'  \href{https://onlinelibrary.wiley.com/doi/full/10.1111/gcb.12736}{Hamann et al. 2015}. Velocity of climate change algorithms for guiding conservation and management. Global Change Biology, 21, 997-1004. \cr
#'  \href{https://onlinelibrary.wiley.com/doi/full/10.1111/gcb.13665}{Garcia Molinos et al. 2017}. Improving the interpretability of climate landscape metrics: An ecological risk analysis of Japan's Marine Protected Areas. Global Change Biology, 23, 4440-4452.
#'
#' @seealso{\code{\link{climPCA}}, \code{\link{climPlot}}}
#'
#' @export
#' @author Jorge Garcia Molinos
#' @examples
#' \dontrun{
#' JapTC <- VoCC_get_data("JapTC.tif")
#'
#' # Create a data frame with the necessary variables in the required order
#' clim <- stats::na.omit(data.frame(terra::values(JapTC), cid = 1:terra::ncell(JapTC)))
#' clim[, c("x", "y")] <- terra::xyFromCell(JapTC, clim$cid)
#'
#' # Constant threshold, distance-restricted velocity based on geographical distances
#' avocc1 <- dVoCC(clim,
#'   n = 3, tdiff = 40, method = "Single", climTol = c(10, 0.1, 0.1),
#'   geoTol = 160, distfun = "GreatCircle", trans = NA, lonlat = TRUE
#' )
#'
#' r1 <- JapTC[[1]]
#' r1[avocc1$focal] <- avocc1$vel
#' terra::plot(r1)
#'
#' # Cell-specific, distance-unrestricted climate analogue velocity based on least-cost path distances
#' # First, create the conductance matrix (all land cells considered to have conductance of 1)
#' r <- JapTC[[1]]
#' r[!is.na(JapTC[[1]])] <- 1
#' h8 <- gdistance::transition(r, transitionFunction = mean, directions = 8)
#' h8 <- gdistance::geoCorrection(h8, type = "c")
#'
#' # Now calculate the analogue velocity using the baseline SD for each variable as analogue threshold
#' avocc2 <- dVoCC(clim,
#'   n = 3, tdiff = 40, method = "Variable", climTol = NA, geoTol = Inf,
#'   distfun = "LeastCost", trans = h8, lonlat = TRUE
#' )
#'
#' # Plot results
#' r1 <- r2 <- JapTC[[1]]
#' r1[avocc1$focal] <- avocc1$vel
#' r2[avocc2$focal] <- avocc2$vel
#' terra::plot(c(r1, r2))
#' }
#'
dVoCC <- function(clim, n, tdiff, method = "Single", climTol, geoTol,
                  distfun = "GreatCircle", trans = NA, lonlat = TRUE) {

  geoDis <- climDis <- ang <- vel <- target <- cid <- a <- NULL # Fix devtools check warnings

  if (distfun == "Euclidean" && lonlat == TRUE) {
    print("Error: Euclidean distances specified for unprojected coordinates")
    stop()
  }

  if (distfun == "LeastCost") {
    stop("LeastCost distances are not currently supported. Use 'Euclidean' or 'GreatCircle' instead.")
  }

  assertthat::assert_that(
    all(distfun %in% c("Euclidean", "GreatCircle")),
    is.na(trans)
  )

  dat <- stats::na.omit(data.table::data.table(clim))

  # matrix with the future climatic values for all cells
  fut <- dat[, seq(2, (2 * n), by = 2), with = FALSE]

  # Determine optimal number of cores, ensuring we don't exceed data rows
  ncores <- parallelly::availableCores(constraints = "connections", omit = 2)
  ncores <- min(ncores, nrow(dat))  # Don't use more cores than data rows
  ncores <- max(ncores, 1)          # Ensure at least 1 core

  # Only use parallel processing if we have multiple cores and sufficient data
  if (ncores > 1 && nrow(dat) > ncores) {
    cuts <- cut(seq_len(nrow(dat)), ncores, labels = FALSE)
    cl <- parallelly::makeClusterPSOCK(ncores, autoStop = TRUE)

    doParallel::registerDoParallel(cl)

    result <- foreach::foreach(a = seq_len(ncores),
                               .combine = rbind,
                               .packages = c("terra", "gdistance", "geosphere", "data.table"),
                               .multicombine = TRUE) %dopar% {

                                 Dat <- dat[cuts == a, ]

                                 resu <- data.table::data.table(
                                   focal = Dat$cid,
                                   target = as.integer(NA),
                                   climDis = as.double(NA),
                                   geoDis = as.double(NA),
                                   ang = as.double(NA),
                                   vel = as.double(NA)
                                 )

                                 i <- 0
                                 while (i <= nrow(Dat)) {
                                   i <- i + 1

                                   # for each focal cell subset target cell analogues (within ClimTol)
                                   pres <- as.numeric(Dat[i, seq(1, (2 * n), by = 2), with = FALSE])
                                   dif <- data.table::data.table(sweep(fut, 2, pres, "-"))

                                   # Identify future analogue cells
                                   if (method == "Single") { # Ohlemuller et al 2006 / Hamann et al 2015
                                     upper <- colnames(dif)
                                     l <- lapply(upper, function(x) call("<", call("abs", as.name(x)), climTol[grep(x, colnames(dif))]))
                                     ii <- Reduce(function(c1, c2) substitute(.c1 & .c2, list(.c1 = c1, .c2 = c2)), l)
                                     anacid <- dat$cid[dif[eval(ii), which = TRUE]] # cids analogue cells
                                   }

                                   if (method == "Variable") { # Garcia Molinos et al. 2017
                                     climTol <- as.numeric(Dat[i, ((2 * n) + 1):(3 * n), with = FALSE]) # focal cell tolerance
                                     upper <- colnames(dif)
                                     l <- lapply(upper, function(x) call("<", call("abs", as.name(x)), climTol[grep(x, colnames(dif))]))
                                     ii <- Reduce(function(c1, c2) substitute(.c1 & .c2, list(.c1 = c1, .c2 = c2)), l)
                                     anacid <- dat$cid[dif[eval(ii), which = TRUE]] # cids analogue cells
                                   }

                                   # LOCATE CLOSEST ANALOGUE
                                   if (length(anacid) > 0) {
                                     # check which of those are within distance and get the analogue at minimum distance
                                     if (distfun == "Euclidean") {
                                       d <- stats::dist(cbind(Dat$x[i], Dat$y[i]), cbind(dat$x[dat$cid %in% anacid], dat$y[dat$cid %in% anacid]))
                                     } # in x/y units
                                     if (distfun == "GreatCircle") {
                                       d <- (geosphere::distHaversine(cbind(Dat$x[i], Dat$y[i]), cbind(dat$x[dat$cid %in% anacid], dat$y[dat$cid %in% anacid]))) / 1000
                                     } # in km

                                     # LeastCost distances not supported - error is thrown at function start

                                     an <- anacid[d < geoTol] # cids analogue cells within search radius
                                     dis <- d[d < geoTol] # distance to candidate analogues
                                     if (length(an) > 0) {
                                       resu[i, target := an[which.min(dis)]] # cid of geographically closest climate analogue
                                       if (method == "Single") {
                                         resu[i, climDis := mean(as.numeric(dif[which(anacid == resu[i, target]), ]))]
                                       } # mean clim difference for the closest analogue
                                       resu[i, geoDis := min(dis)]
                                       resu[i, ang := geosphere::bearing(Dat[i, c("x", "y")], dat[cid == resu[i, target], c("x", "y")])]
                                       resu[i, vel := resu$geoDis[i] / tdiff]
                                     }
                                   }
                                 }
                                 return(resu)
                               }
  } else {
    # Sequential processing for small datasets or limited cores
    result <- data.table::data.table(
      focal = dat$cid,
      target = as.integer(NA),
      climDis = as.double(NA),
      geoDis = as.double(NA),
      ang = as.double(NA),
      vel = as.double(NA)
    )

    for (i in seq_len(nrow(dat))) {
      # for each focal cell subset target cell analogues (within ClimTol)
      pres <- as.numeric(dat[i, seq(1, (2 * n), by = 2), with = FALSE])
      dif <- data.table::data.table(sweep(fut, 2, pres, "-"))

      # Identify future analogue cells
      if (method == "Single") { # Ohlemuller et al 2006 / Hamann et al 2015
        upper <- colnames(dif)
        l <- lapply(upper, function(x) call("<", call("abs", as.name(x)), climTol[grep(x, colnames(dif))]))
        ii <- Reduce(function(c1, c2) substitute(.c1 & .c2, list(.c1 = c1, .c2 = c2)), l)
        anacid <- dat$cid[dif[eval(ii), which = TRUE]] # cids analogue cells
      }

      if (method == "Variable") { # Garcia Molinos et al. 2017
        climTol <- as.numeric(dat[i, ((2 * n) + 1):(3 * n), with = FALSE]) # focal cell tolerance
        upper <- colnames(dif)
        l <- lapply(upper, function(x) call("<", call("abs", as.name(x)), climTol[grep(x, colnames(dif))]))
        ii <- Reduce(function(c1, c2) substitute(.c1 & .c2, list(.c1 = c1, .c2 = c2)), l)
        anacid <- dat$cid[dif[eval(ii), which = TRUE]] # cids analogue cells
      }

      # LOCATE CLOSEST ANALOGUE
      if (length(anacid) > 0) {
        # check which of those are within distance and get the analogue at minimum distance
        if (distfun == "Euclidean") {
          d <- stats::dist(cbind(dat$x[i], dat$y[i]), cbind(dat$x[dat$cid %in% anacid], dat$y[dat$cid %in% anacid]))
        } # in x/y units
        if (distfun == "GreatCircle") {
          d <- (geosphere::distHaversine(cbind(dat$x[i], dat$y[i]), cbind(dat$x[dat$cid %in% anacid], dat$y[dat$cid %in% anacid]))) / 1000
        } # in km

        an <- anacid[d < geoTol] # cids analogue cells within search radius
        dis <- d[d < geoTol] # distance to candidate analogues
        if (length(an) > 0) {
          result[i, target := an[which.min(dis)]] # cid of geographically closest climate analogue
          if (method == "Single") {
            result[i, climDis := mean(as.numeric(dif[which(anacid == result[i, target]), ]))]
          } # mean clim difference for the closest analogue
          result[i, geoDis := min(dis)]
          result[i, ang := geosphere::bearing(dat[i, c("x", "y")], dat[cid == result[i, target], c("x", "y")])]
          result[i, vel := result$geoDis[i] / tdiff]
        }
      }
    }
  }

  return(result)
}

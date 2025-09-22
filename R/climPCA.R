#' Reduce dimensionality of climate predictors via Principal Component Analysis
#'
#' Function to extract the first n principal components explaining a predefined total amount of variance among climatic variables.
#' These components can subsequently be used as synthetic climatic variables to reduce dimensionality in climate-analogue methods.
#'
#' @param climp \code{raster.stack} with one layer for each climatic variable with the values for present or baseline conditions.
#' @param climf \code{raster.stack} with one layer for each climatic variable with the values for future conditions.
#' @param trans \code{function} specifying the type of transformation to be applied prior to the PCA.
#' Specify NA where no transformation is required (default log(x)).
#' @param cen \code{logical} should the variables be centered prior to the PCA? (default TRUE).
#' @param sc \code{logical} should the variables be scaled prior to the PCA? (default TRUE).
#' @param th \code{numeric} threshold giving the minimum amount of total variance that should be explained by the principal components extracted.
#'
#' @return a \code{list} containing (i) the output from the PCA (call to 'prcomp'), and
#' (ii) a table with the present/future cell values for the principal components accounting
#' for the specified percentage of total variance (th).
#'
#' @seealso{\code{\link{dVoCC}}, \code{\link{climPlot}}}
#'
#' @importFrom data.table ":="
#'
#' @export
#'
#' @examples
#'
#' JapTC <- terra::rast(system.file("extdata", "JapTC.tif", package = "VoCC"))
#'
#' comp <- climPCA(JapTC[[c(1, 3, 5)]], JapTC[[c(2, 4, 6)]],
#'                 trans = NA, cen = TRUE, sc = TRUE, th = 0.85)
#' summary(comp[[1]]) # first two components explain >90% of variance
#' # Create a data frame with the necessary variables in the required order (see climAna? for details)
#' clim <- comp[[2]][, c(2, 4, 3, 5, 1)]
#' clim[, c("x", "y")] <- terra::xyFromCell(JapTC[[1]], clim$cid)
#'
climPCA <- function(climp, climf, trans = function(x) log(x), cen = TRUE, sc = TRUE, th = 0.8) {

  .SD <- NULL

  # OPTIMIZATION: Pre-calculate cell count to avoid repeated calls
  n_cells <- terra::ncell(climp[[1]])

  # get a data table with the pooled values (current/future) of the clim variables
  clim <- data.table::data.table(rbind(terra::values(climp), terra::values(climf)))
  clim <- stats::na.omit(clim[, c("cid", "p") := list(rep(1:n_cells, times = 2), rep(c("present", "future"), each = n_cells))])

  # OPTIMIZATION: Store column indices to avoid repeated column selection
  data_cols <- !names(clim) %in% c("cid", "p")

  # apply transformation if required
  if (!is.na(trans)) {
    clim_data <- clim[, .SD, .SDcols = data_cols]
    clim_data <- trans(clim_data)
    # Rebuild clim with transformed data
    clim <- cbind(clim_data, clim[, c("cid", "p")])
  }

  # apply PCA
  clim.pca <- stats::prcomp(clim[, .SD, .SDcols = data_cols], center = cen, scale. = sc)

  # OPTIMIZATION: Vectorized variance calculation
  sdev_squared <- clim.pca$sdev^2
  cumvar_prop <- cumsum(sdev_squared) / sum(sdev_squared)
  a <- which(cumvar_prop >= th)[1]

  val.pca <- clim.pca$x[, 1:a, drop = FALSE]
  val <- data.frame(val.pca, cid = clim$cid, p = clim$p)

  # put it back in wide form
  v <- stats::reshape(val, idvar = "cid", timevar = "p", direction = "wide")
  return(list(clim.pca, v))
}

#' Binned scatter plot for 2-dimensional climate space
#'
#' Function to create a binned scatter plot of two climate variables.
#'
#' @param xy \code{data.frame} with cells as rows and 4 columns representing the present and future local values for the two variables (V1p, V1f, V2p, V2f).
#' @param x.binSize \code{numeric} the bin size for the first variable.
#' @param y.binSize \code{numeric} the bin size for the second variable.
#' @param x.name \code{character} the variable name for the first variable. Used to label the plot.
#' @param y.name \code{character} the variable name for the second variable. Used to label the plot.
#'
#' @return A series of \code{plot} objects displaying the (i) present and (ii) future
#' cell frequency for each combination of local climates,
#' and (iii) the location of remnant, novel and disappearing climates between both periods.
#'
#' @seealso{\code{\link{dVoCC}}, \code{\link{climPCA}}}
#'
#' @export
#' @examples
#'
#' JapTC <- terra::rast(system.file("extdata", "JapTC.tif", package = "VoCC"))
#'
#' # Plot climate space for the two first variables(annual precipitation and maximum temperature)
#' xy <- stats::na.omit(data.frame(
#'   terra::values(JapTC[[1]]),
#'   terra::values(JapTC[[2]]),
#'   terra::values(JapTC[[3]]), terra::values(JapTC[[4]])
#' ))
#'
#' out <- climPlot(xy,
#'   x.binSize = 5, y.binSize = 0.2, x.name = "Precipitation (mm)",
#'   y.name = "Temperature max (°C)"
#' )
#'
#' # output plots can be saved as:
#' ggplot2::ggsave(
#'   plot = out, filename = file.path(getwd(), "example_plot.pdf"),
#'   width = 17, height = 17, unit = "cm"
#' )
#'
climPlot <- function(xy, x.binSize, y.binSize, x.name = "V1", y.name = "V2") {
  xp <- xy[, 1]
  yp <- xy[, 3]
  xf <- xy[, 2]
  yf <- xy[, 4]

  # OPTIMIZATION: Pre-calculate ranges to avoid repeated min/max calls
  x_combined <- c(xp, xf)
  y_combined <- c(yp, yf)
  x_range <- range(x_combined)
  y_range <- range(y_combined)

  # bins per axis
  x.nbins <- floor(abs(x_range[2] - x_range[1]) / x.binSize)
  y.nbins <- floor(abs(y_range[2] - y_range[1]) / y.binSize)
  x.bin <- seq(floor(x_range[1]), ceiling(x_range[2]), length = x.nbins)
  y.bin <- seq(floor(y_range[1]), ceiling(y_range[2]), length = y.nbins)

  # define palette
  rf <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(11, "Spectral")))
  r <- rf(64)

  # present
  freq <- as.data.frame(table(findInterval(xp, x.bin), findInterval(yp, y.bin)))
  freq[, 1] <- as.numeric(freq[, 1])
  freq[, 2] <- as.numeric(freq[, 2])
  freq2Dp <- diag(nrow = x.nbins, ncol = y.nbins)
  freq2Dp[cbind(freq[, 1], freq[, 2])] <- freq[, 3]
  freq2Dp[freq2Dp[] == 0] <- NA

  # future
  freq <- as.data.frame(table(findInterval(xf, x.bin), findInterval(yf, y.bin)))
  freq[, 1] <- as.numeric(freq[, 1])
  freq[, 2] <- as.numeric(freq[, 2])
  freq2Df <- diag(nrow = x.nbins, ncol = y.nbins)
  freq2Df[cbind(freq[, 1], freq[, 2])] <- freq[, 3]
  freq2Df[freq2Df[] == 0] <- NA

  # Adjust values above the upper limit (minimum of the two maxima) to the upper limit so image.plot plot all on the same color scale
  UL <- min(max(freq2Dp, na.rm = TRUE), max(freq2Df, na.rm = TRUE))
  freq2Dp[freq2Dp > UL] <- UL
  freq2Df[freq2Df > UL] <- UL

  # OPTIMIZATION: Vectorized climate classification - eliminates nested loops
  freq2D <- matrix(NA, nrow = x.nbins, ncol = y.nbins)

  # Vectorized logical operations - much faster than nested loops
  present_na <- is.na(freq2Dp)
  future_na <- is.na(freq2Df)

  # Novel climates: in future but not present
  freq2D[present_na & !future_na] <- 1
  # Disappearing climates: in present but not future
  freq2D[!present_na & future_na] <- 2
  # Remnant climates: in both present and future
  freq2D[!present_na & !future_na] <- 0
  # NA remains NA (neither present nor future)

  # plot climate space
  Freq2Dpf <- rbind(
    data.frame(x = x.bin, y = rep(y.bin, each = length(x.bin)), freq = c(freq2Dp)),
    data.frame(x = x.bin, y = rep(y.bin, each = length(x.bin)), freq = c(freq2Df))
  )
  Freq2Dpf$clim <- factor(rep(c("Historical", "Present"), each = nrow(Freq2Dpf) / 2), levels = c("Historical", "Present"))
  Freq2Dpf <- Freq2Dpf[!is.na(Freq2Dpf$freq), ]

  Freq2D <- factor(c("Novel climate", "Remnant climate", "Disappearing climate")[c(freq2D) + 1],
    levels = c("Novel climate", "Remnant climate", "Disappearing climate")
  )
  Freq2D <- data.frame(x = x.bin, y = rep(y.bin, each = length(x.bin)), freq = Freq2D)
  Freq2D <- Freq2D[!is.na(Freq2D$freq), ]

  panelAB <- ggplot2::ggplot(Freq2Dpf, ggplot2::aes(x = .data$x, y = .data$y, fill = freq)) +
    ggplot2::geom_raster() +
    ggplot2::scale_fill_gradientn(
      colors = r,
      name = "Cell count",
      breaks = seq(0, UL, 20),
      guide = ggplot2::guide_colorbar(ticks.linewidth = 1.5)
    ) +
    ggplot2::facet_wrap(~clim, scale = "free_y", ncol = 2) +
    ggplot2::scale_x_continuous(limits = c(min(x.bin), max(x.bin))) +
    ggplot2::scale_y_continuous(limits = c(min(y.bin), max(y.bin))) +
    ggplot2::labs(x = x.name, y = y.name) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.justification = c(1, 0),
      legend.key.width = ggplot2::unit(2, "cm"),
      strip.text = ggplot2::element_blank()
    )

  panelC <- ggplot2::ggplot(Freq2D, ggplot2::aes(x = .data$x, y = .data$y, fill = freq)) +
    ggplot2::geom_raster() +
    ggplot2::scale_fill_manual(values = c("#56B4E9", "#009E73", "#D55E00"), name = "Climate type") +
    ggplot2::labs(x = x.name, y = y.name)

  panels <- cowplot::plot_grid(panelAB, panelC, nrow = 2, rel_heights = c(1.3, 1.0))

  return(panels)
}

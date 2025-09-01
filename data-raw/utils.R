#' Return data location and data if requested
#'
#' @param path Filename of data if known.
#'
#' @return The path as a character string, or a dataset
#' @export
#'
#' @examples
#' \dontrun{
#' VoCC_get_data()
#'
#' HSST <- VoCC_get_data(path = "HSST.tif")
#' }
VoCC_get_data <- function(path = NULL) {
  if (is.null(path)) {
    dir(system.file("extdata", package = "VoCC"))
  } else {
    if (substr(path, nchar(path) - 3 + 1, nchar(path)) == "tif") {
      terra::rast(system.file("extdata", path, package = "VoCC", mustWork = TRUE))
    } else {
      terra::vect(system.file("extdata", path, package = "VoCC", mustWork = TRUE))
    }
  }
}

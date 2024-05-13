## code to prepare `DATASET` dataset goes here

# We don't currently have the providence of these datasets.
# Use the existing ones for now. I will speak to Jorge.

load("data-raw/EEZ.RData")
EEZ <- terra::vect(EEZ)
terra::writeVector(EEZ, filename = "inst/extdata/EEZ.gpkg", overwrite = TRUE)

load("data-raw/HSST.RData")
HSST <- terra::rast(HSST)
terra::writeRaster(HSST, filename = "inst/extdata/HSST.tif", overwrite = TRUE)

load("data-raw/JapTC.rda")
raster::writeRaster(JapTC, filename = "data-raw/JapTC.tif", overwrite = TRUE) # Error converting directly to rast
JapTC <- terra::rast("data-raw/JapTC.tif")
terra::writeRaster(JapTC, filename = "inst/extdata/JapTC.tif")

# usethis::use_data(EEZ, HSST, JapTC, internal = FALSE, overwrite = TRUE, compress = "bzip2")

# Check suggested compression formats
# tools::checkRdaFiles("R")

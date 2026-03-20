## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
collapse = TRUE,
comment = "#>"
)

# Reduce memory usage and optimize temp file handling
terra::terraOptions(memfrac = 0.5, tempdir = tempdir())

# Force sequential processing to avoid parallel processing issues during package checking
Sys.setenv("R_PARALLELLY_AVAILABLECORES_FALLBACK" = "1")

## ----message=FALSE------------------------------------------------------------
library(VoCC)
library(dplyr)
library(tidyr)
library(data.table)
library(terra)
library(purrr)
library(future)

# For Plotting
library(ggplot2)
library(tidyterra)
library(patchwork)
library(mapplots)

## ----message=FALSE------------------------------------------------------------
library(VoCCdata)

## -----------------------------------------------------------------------------
str(marshift)

## -----------------------------------------------------------------------------
HadiSST <- terra::rast(system.file("extdata", "HadiSST.tif", package = "VoCCdata"))

# monthly to annual averages
r <- sumSeries(HadiSST, p = "1960-01/2009-12", yr0 = "1955-01-01", 
               l = terra::nlyr(HadiSST), 
               fun = function(x) colMeans(x, na.rm = TRUE), 
               freqin = "months", freqout = "years")

vt <- tempTrend(r, th = 10) # temporal trend
vg <- spatGrad(r, th = 0.0001, projected = FALSE) # spatial gradient
gv <- gVoCC(vt, vg) # climate velocity

# Now the distance-based velocities
# Take 1960-1970 as base period against 2000-2009
r2 <- c(terra::mean(r[[1:10]], na.rm = TRUE), terra::mean(r[[41:50]], na.rm = TRUE))

# prepare the data frame with the necessary variables
clim <- na.omit(data.frame(terra::values(r2), cid = 1:terra::ncell(r)))

clim[, c("x", "y")] <- terra::xyFromCell(r, clim$cid)


future::plan(future::multisession, workers = parallelly::availableCores(omit = 2))

# 1965-2004 (40 yr), 500 km search radius
v <- dVoCC(clim, n = 1, tdiff = 40, method = "Single", climTol = 0.1, geoTol = 500, distfun = "GreatCircle", trans = NA, lonlat = TRUE)

future::plan(future::sequential)


## -----------------------------------------------------------------------------
# Change sign as needed and create the distance-based velocity raster
# Change sign as needed - terra approach for value comparison
focal_vals <- terra::values(r2[[1]])[v$focal]
target_vals <- terra::values(r2[[2]])[v$target]
ind <- which(focal_vals > target_vals)
v$velBis <- v$vel
v$velBis[ind] <- v$vel[ind] * -1

# put output in raster format - create single layer empty template like raster(gv)
dv <- terra::rast(terra::ext(gv), resolution = terra::res(gv), crs = terra::crs(gv))
dv[v$focal] <- v$velBis

# Create point geometries and buffer them
coords <- terra::vect(cbind(marshift$long, marshift$lat), crs = "EPSG:4326")
buffer_size <- marshift$Shift * (marshift$timespan / 10) * 1000

# Get the mean velocity within the buffer for each data point.
# Match old raster::extract approach exactly
marshift$GV <- terra::extract(abs(gv[[1]]), coords,
                              buffer = buffer_size,
                              fun = mean, na.rm = TRUE,
                              weights = TRUE, exact = FALSE)[,2]

marshift$DV <- terra::extract(abs(dv), coords,
                              buffer = buffer_size,
                              fun = mean, na.rm = TRUE,
                              weights = TRUE, exact = FALSE)[,2]

## -----------------------------------------------------------------------------
# For points that didn't get values (NA), find nearest valid cells

missing_points <- coords[is.na(marshift$GV)] # Identify NAs
if(!is.empty(missing_points)){
  marine_cells <- terra::as.points(gv[[1]]) # vector of all valid marine cell locations.
  nearest_indices <- terra::nearest(missing_points, marine_cells) # Find the nearest marine cell
  nearest_values <- terra::extract(gv[[1]], marine_cells[nearest_indices]) # get the values from the nearest marine cells.
  marshift$GV[is.na(marshift$GV)] <- nearest_values[, 2] # Replace the NA values in `marshift$GV`
}

missing_points <- coords[is.na(marshift$DV)] # Identify NAs
if(!is.empty(missing_points)){
  marine_cells <- terra::as.points(dv[[1]]) # vector of all valid marine cell locations.
  nearest_cells <- terra::nearest(missing_points, marine_cells) # Find the nearest marine cell
  nearest_values <- terra::extract(dv[[1]], nearest_cells) # get the values from the nearest marine cells.
  marshift$DV[is.na(marshift$DV)] <- nearest_values[, 2] # Replace the NA values in `marshift$GV`
}


## -----------------------------------------------------------------------------
# fit the regression models
Mgv <- lm(Shift^(1 / 4) ~ I((GV * 10)^(1 / 4)), data = marshift, weights = years_data)
summary(Mgv)

## -----------------------------------------------------------------------------
Mdv <- lm(Shift^(1 / 4) ~ I((DV * 10)^(1 / 4)), data = marshift, weights = years_data)
summary(Mdv)

## -----------------------------------------------------------------------------
# first compare both velocities

p1 <- ggplot() +
  geom_spatraster(data = gv[[1]]) +
  scale_fill_distiller(palette = "RdBu", direction = -1, limits = c(-50, 50)) +
  ggtitle("Gradient-based vocc") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))

p2 <- ggplot() +
  geom_spatraster(data = dv[[1]]) +
  scale_fill_distiller(palette = "RdBu", direction = -1, limits = c(-20, 20)) +
  ggtitle("Distance-based vocc") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))

wrap_plots(p1, p2, ncol = 1)


## -----------------------------------------------------------------------------
# scatter plots with the resulting regression line
p1 <- ggplot(na.omit(marshift), aes(x = (GV * 10)^(1 / 4), y = Shift^(1 / 4))) +
  geom_point(color = "grey") +
  geom_smooth(method = lm, se = FALSE) +
  theme_classic() +
  scale_color_brewer(palette = "Accent") +
  labs(x = "Predicted shift (x^1/4; km/yr)", y = "Observed shift (y^1/4; km/yr)")

p2 <- ggplot(na.omit(marshift), aes(x = (DV * 10)^(1 / 4), y = Shift^(1 / 4))) +
  geom_point(color = "grey") +
  geom_smooth(method = lm, se = FALSE) +
  theme_classic() +
  scale_color_brewer(palette = "Accent") +
  labs(x = "Predicted shift (x^1/4; km/yr)", y = "Observed shift (y^1/4; km/yr)")

wrap_plots(p1, p2, nrow = 1)

## -----------------------------------------------------------------------------
# prepare raster layers
vel <- gv[[1]]
ang <- gv[[2]]
mn <- terra::app(r, mean, na.rm = T)

# generate a velocity layer centered and cropped to study region to extract the initial coordinates for the trajectories from
x1 <- terra::crop(gv[[1]], terra::ext(-180, 0, -90, 90))
x2 <- terra::crop(gv[[1]], terra::ext(0, 180, -90, 90))
terra::ext(x1) <- c(180, 360, -90, 90)
velc <- terra::merge(x1, x2)

# crop to the desired extent
# display restricted to +180 longitude to avoid plotting issues with date line crossing
velc <- terra::crop(velc, c(90, 180, -32, 33))

## -----------------------------------------------------------------------------
lonlat <- data.frame(terra::xyFromCell(velc, 1:ncell(velc)))
lonlat$vel <- terra::extract(vel, lonlat, ID = FALSE)
lonlat$ang <- terra::extract(ang, lonlat[, 1:2], ID = FALSE)
lonlat$mn <- terra::extract(mn, lonlat[, 1:2], ID = FALSE)
lonlat$lineID <- 1:nrow(lonlat)
lonlat <- drop_na(lonlat)


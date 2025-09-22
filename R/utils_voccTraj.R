# Internal helper functions ----------------------------------------------------

# Function for finding minimum distance
mfind <- function(rw) {
  X <- which.min(rw)
  return(X)
}

# Function to extract corner coordinates per row
mplace <- function(rw) {
  X <- rw[(rw[1] + 1)]
  Y <- rw[rw[1]]
  return(c(X, Y))
}

# Simple circular functions
deg2rad <- function(deg) {
  return(deg * pi / 180)
}
rad2deg <- function(rad) {
  return(180 * rad / pi)
}


#' Find distances from departure point to corners of ttcell
#' OPTIMIZED: Use direct Haversine formula instead of expensive sf operations
#'
#' @noRd
get_dist <- function(y1, x1, y2, x2) {
  # Convert degrees to radians
  lat1_rad <- y1 * pi / 180
  lon1_rad <- x1 * pi / 180
  lat2_rad <- y2 * pi / 180
  lon2_rad <- x2 * pi / 180
  
  # Haversine formula - much faster than sf::st_distance
  dlat <- lat2_rad - lat1_rad
  dlon <- lon2_rad - lon1_rad
  
  a <- sin(dlat/2)^2 + cos(lat1_rad) * cos(lat2_rad) * sin(dlon/2)^2
  c <- 2 * atan2(sqrt(a), sqrt(1-a))
  
  # Earth radius in meters (same units as sf::st_distance)
  R <- 6371000
  distance <- R * c
  
  return(distance)
}


#' Function for to find closest cooler/warmer cell within the union of two levels of 8-cell adjacency from the "departed" and "destination" cells
#' Use when working with the coarse grid!
#' JDE - Add bfr as an argument to voccTraj?
#' A from-to cell pair and coordinates of departure point; bfr is a search radius in km for nearest cooler cell, if you end on land
#'
#' @noRd
get_dest_cell_coarse <- function(rw, x_res, y_res, bfr, vel_values, mn_values, vel_raster, mn_raster) {
  
  # Add error handling for invalid inputs
  if (any(is.na(rw[1:4]))) {
    return(NA)
  }

  # Find clumps of ocean; the rule is that you can't jump from one clump to another,
  # because this would mean passing over unbroken land
  pos_depart <- tryCatch({
    data.frame(x = purrr::pluck(rw, 3),
               y = purrr::pluck(rw, 4)) %>% # Departure coordinates
      sf::st_as_sf(coords = c("x", "y"), crs = "EPSG:4326")
  }, error = function(e) {
    return(NULL)
  })
  
  if (is.null(pos_depart)) {
    return(NA)
  }

  clumped <- tryCatch({
    get_clumps(pos_depart, mn_raster, bfr, x_res, y_res)
  }, error = function(e) {
    return(NULL)
  })
  
  if (is.null(clumped)) {
    return(NA)
  }

  # Which clump did I start in?
  r1 <- purrr::pluck(rw, 1) # We use this cell a lot, so let's just make it an object
  from_clump <- terra::extract(clumped, terra::xyFromCell(mn_raster, r1)) %>%
    unlist()

  # What are the coordinates of cells within the search area that fall in the clump I came from?
  search_xy <- terra::xyFromCell(clumped, which(clumped[] == from_clump)) %>%
    as.data.frame()

  # Get the ssts in the cells to search
  or <- terra::extract(mn_raster, search_xy, cells = TRUE, xy = TRUE) %>%
    dplyr::rename(sst = 2, x = 4, y = 5)

  # MEMORY LEAK FIX: Use pre-extracted values instead of raster indexing
  if (vel_values[r1] > 0) {
    o <- or %>%
      dplyr::filter(.data$sst < mn_values[r1]) %>%
      stats::na.omit()
    if (nrow(o) == 0) {
      dest_cell <- NA
    } else {
      # OPTIMIZED: Use vectorized Haversine instead of sf::st_distance
      depart_coords <- sf::st_coordinates(pos_depart)
      distances <- get_dist(depart_coords[1,2], depart_coords[1,1], o$y, o$x)
      dest_cell <- o$cell[which.min(distances)]
    }
  } else { # Otherwise, find nearest warmer cell, if there is one
    o <- or %>%
      dplyr::filter(.data$sst > mn_values[r1]) %>%
      stats::na.omit()
    if (nrow(o) == 0) {
      dest_cell <- NA
    } else {
      # OPTIMIZED: Use vectorized Haversine instead of sf::st_distance
      depart_coords <- sf::st_coordinates(pos_depart)
      distances <- get_dist(depart_coords[1,2], depart_coords[1,1], o$y, o$x)
      dest_cell <- o$cell[which.min(distances)]
    }
  }
  return(dest_cell)
}





#' Function to find closest cooler/warmer cell within the union of two levels of 8-cell adjacency from the "departed" and "destination" cells
#'
#' @noRd
get_dest_cell_fine <- function(rw, x_res, y_res, bfr, vel_values, mn_values, vel_raster, mn_raster) {
  
  # Add error handling for invalid inputs
  if (any(is.na(rw[1:4]))) {
    return(NA)
  }

  # Find clumps of ocean; the rule is that you can't jump from one clump to another, because this would mean passing over unbroken land
  pos_depart <- tryCatch({
    data.frame(
      x = purrr::pluck(rw, 3),
      y = purrr::pluck(rw, 4)
    )
  }, error = function(e) {
    return(NULL)
  })
  
  if (is.null(pos_depart)) {
    return(NA)
  }

  xy <- tryCatch({
    terra::xyFromCell(mn_raster, as.vector(as.matrix(rw[1:2]))) %>% # Coordinates of cell centres for start cell and end cell
      as.data.frame() %>%
      sf::st_as_sf(coords = c("x", "y"), crs = "EPSG:4326")
  }, error = function(e) {
    return(NULL)
  })
  
  if (is.null(xy)) {
    return(NA)
  }

  clumped <- tryCatch({
    get_clumps(xy, mn_raster, bfr, x_res, y_res)
  }, error = function(e) {
    return(NULL)
  })
  
  if (is.null(clumped)) {
    return(NA)
  }

  # Which clump did I start in?
  r1 <- purrr::pluck(rw, 1) # We use this cell a lot, so let's just make it an object
  from_clump <- purrr::pluck(terra::extract(clumped, terra::xyFromCell(mn_raster, r1)), 1)

  # What are the coordinates of cells within the search area that fall in the clump I came from?
  search_xy <- terra::xyFromCell(clumped, which(clumped[] == from_clump)) %>%
    as.data.frame()
  search_cells <- terra::cellFromXY(mn_raster, search_xy) # Which cells are these

  # Get the ssts in the cells to search
  or <- terra::extract(mn_raster, search_xy, cells = TRUE, xy = TRUE) %>%
    dplyr::rename(sst = 2, x = 4, y = 5)

  # MEMORY LEAK FIX: Use pre-extracted values instead of raster indexing
  if (vel_values[r1] > 0) {
    # Find all cells in the search area that meet the sst criterion
    o <- or %>%
      dplyr::filter(.data$sst < mn_values[r1]) %>%
      stats::na.omit()
    if (nrow(o) == 0) {
      dest_cell <- NA # Set condition with which to ID stuck cells (i.e., if there are no suitable cells in the search area)
    } else {
      # OPTIMIZED: Use vectorized Haversine instead of sf::st_distance
      distances <- get_dist(pos_depart$y, pos_depart$x, o$y, o$x)
      potential_dest_cells <- o %>%
        dplyr::mutate(distances = distances) %>%
        dplyr::arrange(.data$distances)
      n_cells <- 25 # Number of "close" cells to pick among *** Can change this, or set as an argument (just remove the hard-code here)
      if (nrow(potential_dest_cells) >= n_cells) {
        dest_cell <- potential_dest_cells %>%
          dplyr::slice_head(n = n_cells) %>% # Get the closest cells
          dplyr::slice_sample(n = 1) %>% # Remove .preserve = TRUE which may cause seed issues
          dplyr::pull(.data$cell)
      } else { # If fewer than n_cells, sample 10% of what you have
        n_sample <- max(1, ceiling(nrow(potential_dest_cells) / 10))
        dest_cell <- potential_dest_cells %>%
          dplyr::slice_head(n = n_sample) %>%
          dplyr::slice_sample(n = 1) %>% # Remove .preserve = TRUE which may cause seed issues
          dplyr::pull(.data$cell)
      }
    }
  } else { # Otherwise, find nearest cells that is warmer, if there are any
    # Find all cells in the search area that meet the sst criterion
    o <- or %>%
      dplyr::filter(.data$sst > mn_values[r1]) %>%
      stats::na.omit()
    if (nrow(o) == 0) {
      dest_cell <- NA # Set condition with which to ID stuck cells (i.e., if there are no suitable cells in the search area)
    } else {
      # OPTIMIZED: Use vectorized Haversine instead of sf::st_distance
      distances <- get_dist(pos_depart$y, pos_depart$x, o$y, o$x)
      potential_dest_cells <- o %>%
        dplyr::mutate(distances = distances) %>%
        dplyr::arrange(.data$distances)

      n_cells <- 25 # Number of "close" cells to pick among *** Can change this, or set as an argument (just remove the hard-code here)
      if (nrow(potential_dest_cells) >= n_cells) {
        dest_cell <- potential_dest_cells %>%
          dplyr::slice_head(n = n_cells) %>% # Get the closest cells
          dplyr::slice_sample(n = 1) %>% # Remove .preserve = TRUE which may cause seed issues
          dplyr::pull(.data$cell)
      } else { # If fewer than n_cells, sample 10% of what you have
        n_sample <- max(1, ceiling(nrow(potential_dest_cells) / 10))
        dest_cell <- potential_dest_cells %>%
          dplyr::slice_head(n = n_sample) %>%
          dplyr::slice_sample(n = 1) %>% # Remove .preserve = TRUE which may cause seed issues
          dplyr::pull(.data$cell)
      }
    }
  }
  return(dest_cell)
}


# JDE - I think we should rearrange arguments to x, y to be consistent

# Find new destination, given velocity (km/yr), angle (º), time step (yr) and initial coordinates (ºlon, ºlat); max allowed jump is 1 deg
# vell = vel[fcells] %>% pull(1); angg = ang[fcells] %>% pull(1); timestep = tstep; ll = llold
#'
#' @noRd
destcoords <- function(vell, angg, timestep, ll, y_res, x_res) {
  # OPTIMIZED: Pre-compute trigonometric values to avoid repeated calculations
  angg_rad <- deg2rad(angg)
  cos_angg <- cos(angg_rad)
  sin_angg <- sin(angg_rad)
  
  latshift <- (abs(vell) * timestep * cos_angg) / 111.325 # Calculate shift in lat
  latnew <- ll[, 2] + latshift # Find new lat...first approximation
  
  # OPTIMIZED: Pre-compute cos(deg2rad(latnew)) for longitude calculation
  cos_latnew <- cos(deg2rad(latnew))
  lonshift <- (abs(vell) * timestep * sin_angg) / (111.325 * cos_latnew) # Shift in lon

  # Limit lonshift to at most 1 x resolution cell
  lonshift[lonshift > x_res] <- x_res
  lonshift[lonshift < -x_res] <- -x_res

  # Now on that basis, adjust latshift
  x_gt <- which(lonshift == x_res) # Indices for adjusted lon shifts
  if (length(x_gt) > 0) {
    cos_ll_gt <- cos(deg2rad(ll[x_gt, 2]))
    tan_angg_gt <- tan(angg_rad[x_gt])
    latshift[x_gt] <- (x_res * 111.325 * cos_ll_gt / tan_angg_gt) / 111.325
  }
  
  x_lt <- which(lonshift == -x_res) # Indices for adjusted lon shifts
  if (length(x_lt) > 0) {
    cos_ll_lt <- cos(deg2rad(ll[x_lt, 2]))
    tan_angg_lt <- tan(angg_rad[x_lt])
    latshift[x_lt] <- (x_res * 111.325 * cos_ll_lt / tan_angg_lt) / 111.325
  }
  
  latnew <- ll[, 2] + latshift # Find new lat by adding the adjusted lats

  # Stop new lat from jumping the poles
  latnew[latnew > 90] <- 90
  latnew[latnew < -90] <- -90

  # Adjust lon
  lonnew <- ll[, 1] + lonshift # Find new lon...first approximation

  # Adjust for dateline jumps
  lonnew <- lonnew - (360 * floor((lonnew + 180) / 360))

  # OPTIMIZED: Direct data.frame creation instead of pipe
  result <- data.frame(dlon = lonnew, dlat = latnew)
  return(result)
}


#'
#' @noRd
get_clumps <- function(xy, mn, bfr, x_res, y_res){
  
  # MEMORY LEAK FIX: Add error handling and explicit memory management
  tryCatch({
    # MEMORY LEAK FIX: Limit buffer size to prevent excessive memory usage
    max_bfr <- min(bfr, 200)  # Cap buffer at 200km to prevent memory issues
    sp_buffer <- sf::st_buffer(xy, max_bfr * 1000) # Buffer around departure point, remembering that buffer is in metres

    buffer_zone <- terra::extract(mn, sp_buffer, cells = TRUE, xy = TRUE) %>%
      dplyr::select(-"ID") %>%
      dplyr::distinct() %>%
      dplyr::rename(sst = 1) %>% #*** rename "climatology", if needed
      dplyr::select("x", "y", "sst", "cell") %>%
      tidyr::drop_na("sst")
    
    # Check if buffer_zone is empty
    if (nrow(buffer_zone) == 0) {
      return(NULL)
    }
    
    # MEMORY LEAK FIX: Limit buffer_zone size to prevent memory explosion
    if (nrow(buffer_zone) > 10000) {
      # Sample down to manageable size
      buffer_zone <- buffer_zone[sample(nrow(buffer_zone), 10000), ]
    }

    clumped_rast <- terra::rast(
      xmin = min(buffer_zone$x) - x_res/2,
      xmax = max(buffer_zone$x) + x_res/2,
      ymin = min(buffer_zone$y) - y_res/2,
      ymax = max(buffer_zone$y) + y_res/2,
      resolution = c(x_res, y_res),
      crs = "EPSG:4326"
    )

    # MEMORY LEAK FIX: Create SpatVector more efficiently and clean up immediately
    buffer_vect <- terra::vect(buffer_zone[, c("x", "y", "sst")], geom = c("x", "y"), crs = "EPSG:4326")
    
    clumped <- terra::rasterize(
      x = buffer_vect,
      y = clumped_rast, # The template to rasterize onto
      field = "sst") %>% # The data
      terra::patches(directions = 8, allowGaps = FALSE)
    
    # MEMORY LEAK FIX: Explicitly clean up temporary objects
    rm(buffer_vect, clumped_rast, sp_buffer, buffer_zone)
    
    return(clumped)
    
  }, error = function(e) {
    warning("Error in get_clumps: ", e$message)
    return(NULL)
  })
}


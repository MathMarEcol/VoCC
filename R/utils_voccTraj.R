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
#'
#' @noRd
get_dist <- function(y1, x1, y2, x2) {
  pt1 <- sf::st_as_sf(data.frame(x = x1, y = y1), coords = c("x", "y"), crs = sf::st_crs(terra::rast()))
  pt2 <- sf::st_as_sf(data.frame(x = x2, y = y2), coords = c("x", "y"), crs = sf::st_crs(terra::rast()))
  out <- sf::st_distance(pt1, pt2, by_element = TRUE) %>%
    as.vector()
  return(out)
}


#' Function for to find closest cooler/warmer cell within the union of two levels of 8-cell adjacency from the "departed" and "destination" cells
#' Use when working with the coarse grid!
#' JDE - Add bfr as an argument to voccTraj?
#' A from-to cell pair and coordinates of departure point; bfr is a search radius in km for nearest cooler cell, if you end on land
#'
#' @noRd
get_dest_cell_coarse <- function(rw, x_res, y_res, bfr, vel, mn) {
  
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
    get_clumps(pos_depart, mn, bfr, x_res, y_res)
  }, error = function(e) {
    return(NULL)
  })
  
  if (is.null(clumped)) {
    return(NA)
  }

  # Which clump did I start in?
  r1 <- purrr::pluck(rw, 1) # We use this cell a lot, so let's just make it an object
  from_clump <- terra::extract(clumped, terra::xyFromCell(mn, r1)) %>%
    unlist()

  # What are the coordinates of cells within the search area that fall in the clump I came from?
  search_xy <- terra::xyFromCell(clumped, which(clumped[] == from_clump)) %>%
    as.data.frame()

  # Get the ssts in the cells to search
  or <- terra::extract(mn, search_xy, cells = TRUE, xy = TRUE) %>%
    dplyr::rename(sst = 2, x = 4, y = 5)

  # If velocity is positive, find nearest cooler cell, if there is one
  if (vel[r1] > 0) {
    o <- or %>%
      dplyr::filter(.data$sst < unlist(mn[r1])) %>%
      stats::na.omit()
    # %>%
    #   bind_cols(., terra::xyFromCell(mn, .$cell)) #*** Do we need this line?
    if (nrow(o) == 0) {
      dest_cell <- NA
    } else {
      # pt <- terra::xyFromCell(mn, r1) %>%
      #   data.frame()

      dest_cell <- sf::st_distance(
        pos_depart,
        sf::st_as_sf(data.frame(o),
          coords = c("x", "y"),
          crs = "EPSG:4326"
        )
      ) %>%
        which.min() %>%
        {
          o$cell[.]
        } # The closest cell with appropriate sst
    }
  } else { # Otherwise, find nearest warmer cell, if there is one
    o <- or %>%
      dplyr::filter(.data$sst > unlist(mn[r1])) %>%
      stats::na.omit()
    # %>%
    #   bind_cols(., terra::xyFromCell(mn, .$sst)) #*** Do we need this line?
    if (nrow(o) == 0) {
      dest_cell <- NA
    } else {
      # pt <- terra::xyFromCell(mn, r1) %>%
      #   data.frame()

      dest_cell <- sf::st_distance(
        pos_depart,
        sf::st_as_sf(data.frame(o),
          coords = c("x", "y"),
          crs = "EPSG:4326"
        )
      ) %>%
        which.min() %>% # The closest cell with appropriate sst
        {
          o$cell[.]
        }
    }
  }
  return(dest_cell)
}





#' Function to find closest cooler/warmer cell within the union of two levels of 8-cell adjacency from the "departed" and "destination" cells
#'
#' @noRd
get_dest_cell_fine <- function(rw, x_res, y_res, bfr, vel, mn) { # A to-from cell pair and the buffer around a cell on land that can be searched for a suitable cell
  
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
    terra::xyFromCell(mn, as.vector(as.matrix(rw[1:2]))) %>% # Coordinates of cell centres for start cell and end cell
      as.data.frame() %>%
      sf::st_as_sf(coords = c("x", "y"), crs = "EPSG:4326")
  }, error = function(e) {
    return(NULL)
  })
  
  if (is.null(xy)) {
    return(NA)
  }

  clumped <- tryCatch({
    get_clumps(xy, mn, bfr, x_res, y_res)
  }, error = function(e) {
    return(NULL)
  })
  
  if (is.null(clumped)) {
    return(NA)
  }


  # Which clump did I start in?
  r1 <- purrr::pluck(rw, 1) # We use this cell a lot, so let's just make it an object
  from_clump <- purrr::pluck(terra::extract(clumped, terra::xyFromCell(mn, r1)), 1)

  # What are the coordinates of cells within the search area that fall in the clump I came from?
  search_xy <- terra::xyFromCell(clumped, which(clumped[] == from_clump)) %>%
    as.data.frame()
  search_cells <- terra::cellFromXY(mn, search_xy) # Which cells are these

  # Get the ssts in the cells to search
  or <- terra::extract(mn, search_xy, cells = TRUE, xy = TRUE) %>%
    dplyr::rename(sst = 2, x = 4, y = 5)

  # If velocity is positive, find nearest cells that are cooler, if there are any
  if (unlist(vel[r1]) > 0) {
    # Find all cells in the search area that meet the sst criterion
    o <- or %>%
      dplyr::filter(.data$sst < unlist(mn[r1])) %>%
      stats::na.omit()
    if (nrow(o) == 0) {
      dest_cell <- NA # Set condition with which to ID stuck cells (i.e., if there are no suitable cells in the search area)
    } else {
      potential_dest_cells <- sf::st_as_sf(o, coords = c("x", "y"), crs = sf::st_crs(terra::rast())) %>%
        dplyr::mutate(distances = sf::st_distance(
          x = sf::st_as_sf(pos_depart,
            coords = c("x", "y"),
            crs = "EPSG:4326"
          ),
          y = .data
        ) %>%
          as.vector()) %>%
        dplyr::arrange("distances")
      n_cells <- 25 # Number of "close" cells to pick among *** Can change this, or set as an argument (just remove the hard-code here)
      if (nrow(potential_dest_cells) >= n_cells) {
        dest_cell <- potential_dest_cells %>%
          dplyr::slice_head(n = n_cells) %>% # Get the closest cells
          dplyr::slice_sample(n = 1, .preserve = TRUE) %>% # Pick one at random with preserved grouping
          dplyr::pull(.data$cell)
      } else { # If fewer than n_cells, sample 10% of what you have
        n_sample <- max(1, ceiling(nrow(potential_dest_cells) / 10))
        dest_cell <- potential_dest_cells %>%
          dplyr::slice_head(n = n_sample) %>%
          dplyr::slice_sample(n = 1, .preserve = TRUE) %>%
          dplyr::pull(.data$cell)
      }
    }
  } else { # Otherwise, find nearest cells that is warmer, if there are any
    # Find all cells in the search area that meet the sst criterion
    o <- or %>%
      dplyr::filter(.data$sst > unlist(mn[r1])) %>%
      stats::na.omit()
    if (nrow(o) == 0) {
      dest_cell <- NA # Set condition with which to ID stuck cells (i.e., if there are no suitable cells in the search area)
    } else {
      potential_dest_cells <- sf::st_as_sf(o, coords = c("x", "y"), crs = sf::st_crs(terra::rast())) %>%
        dplyr::mutate(distances = sf::st_distance(
          x = sf::st_as_sf(pos_depart,
            coords = c("x", "y"),
            crs = "EPSG:4326"
          ),
          y = .data
        ) %>%
          base::as.vector()) %>%
        dplyr::arrange(.data$distances)

      n_cells <- 25 # Number of "close" cells to pick among *** Can change this, or set as an argument (just remove the hard-code here)
      if (nrow(potential_dest_cells) >= n_cells) {
        dest_cell <- potential_dest_cells %>%
          dplyr::slice_head(n = n_cells) %>% # Get the closest cells
          dplyr::slice_sample(n = 1, .preserve = TRUE) %>% # Pick one at random with preserved grouping
          dplyr::pull(.data$cell)
      } else { # If fewer than n_cells, sample 10% of what you have
        n_sample <- max(1, ceiling(nrow(potential_dest_cells) / 10))
        dest_cell <- potential_dest_cells %>%
          dplyr::slice_head(n = n_sample) %>%
          dplyr::slice_sample(n = 1, .preserve = TRUE) %>%
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
  latshift <- (abs(vell) * timestep * cos(deg2rad(angg))) / 111.325 # Calculate shift in lat
  latnew <- ll[, 2] + latshift # Find new lat...first approximation
  lonshift <- (abs(vell) * timestep * sin(deg2rad(angg))) / (111.325 * cos(deg2rad(latnew))) # Shift in lon

  # Limit lonshift to at most 1 x resolution cell
  lonshift[lonshift > x_res] <- x_res
  lonshift[lonshift < -x_res] <- -x_res

  # Now on that basis, adjust latshift
  x_gt <- which(lonshift == x_res) # Indices for adjusted lon shifts
  latshift[x_gt] <- ((x_res * 111.325 * cos(deg2rad(ll[x_gt, 2]))) / tan(deg2rad(angg[x_gt]))) / 111.325 # Using trig on distances
  x_lt <- which(lonshift == -x_res) # Indices for adjusted lon shifts
  latshift[x_lt] <- ((x_res * 111.325 * cos(deg2rad(ll[x_lt, 2]))) / tan(deg2rad(angg[x_lt]))) / 111.325 # Using trig on distances
  latnew <- ll[, 2] + latshift # Find new lat by adding the adjusted lats

  # Stop new lat from jumping the poles
  latnew[latnew > 90] <- 90
  latnew[latnew < -90] <- -90

  # Adjust lon
  lonnew <- ll[, 1] + lonshift # Find new lon...first approximation

  # Adjust for dateline jumps
  lonnew <- lonnew - (360 * floor((lonnew + 180) / 360))

  return(data.frame(lonnew, latnew) %>%
    stats::setNames(c("dlon", "dlat")))
}


#'
#' @noRd
get_clumps <- function(xy, mn, bfr, x_res, y_res){
  
  # Add error handling and memory management
  tryCatch({
    sp_buffer <- sf::st_buffer(xy, bfr * 1000) # Buffer around departure point, remembering that buffer is in metres

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

    clumped_rast <- terra::rast(
      xmin = min(buffer_zone$x) - x_res/2,
      xmax = max(buffer_zone$x) + x_res/2,
      ymin = min(buffer_zone$y) - y_res/2,
      ymax = max(buffer_zone$y) + y_res/2,
      resolution = c(x_res, y_res),
      crs = "EPSG:4326"
    )

    clumped <- terra::rasterize(
      x = buffer_zone %>% # Needs to be SpatVector to add the sst
        dplyr::select(-"cell") %>%
        terra::vect(geom = c("x", "y"), crs = "EPSG:4326"),
      y = clumped_rast, # The template to rasterize onto
      field = "sst") %>% # The data
      terra::patches(directions = 8, allowGaps = FALSE)
    
    return(clumped)
    
  }, error = function(e) {
    warning("Error in get_clumps: ", e$message)
    return(NULL)
  })
}


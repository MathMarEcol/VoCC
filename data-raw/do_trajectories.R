# Run trajectories for a file ---------------------------------------------

do_trajectories <- function(raster_stack, start_points, nyears) {
  # Get the rasters we need
  stacks <- build_coastal_mask(f_m = "OISST_Mask_data/mask.nc",
                               f_cm = "/Volumes/Jet_drive_too/Coastal_masks/coarseAusmask.nc",
                               f_fm = "/Volumes/Jet_drive_too/Coastal_masks/fineAusmask.nc",
                               rast_stk = raster_stack)
  fine_coast <- stacks[[3]]
  vel <- stacks[[1]][[1]]
  ang <- stacks[[1]][[2]]
  mn <- stacks[[1]][[3]]
  vel_c <- stacks[[2]][[1]]
  ang_c <- stacks[[2]][[2]]
  mn_c <- stacks[[2]][[3]]
  lk_up <- stacks$look_up

  lonlat <- start_points[,1:2] %>%
    mutate(ID = 1:nrow(.)) %>%
    as_tibble() %>%
    mutate(vel = terra::extract(vel, .[,1:2], ID = FALSE) %>%
             pull(voccMag),
           i_coast = terra::extract(fine_coast, .[,1:2], ID = FALSE) %>%
             pull(tos),
           to_drop = ifelse(is.na(vel) & is.na(i_coast), NA, 1)) %>%
    drop_na(to_drop) %>%
    dplyr::select(x, y, ID)
  out <- traject(lonlat,
                 vel = vel, ang = ang, mn = mn,
                 vel_c = vel_c, ang_c = ang_c, mn_c = mn_c,
                 fine_coast = fine_coast,
                 lk_up = lk_up,
                 tstep = 1/12, tyr = nyears)
  return(out)
}

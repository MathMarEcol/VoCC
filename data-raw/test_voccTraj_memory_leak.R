# Memory Leak Test Script for voccTraj
# This script runs voccTraj 10 times with identical arguments
# and monitors memory usage to detect potential memory leaks

devtools::load_all()
library(VoCCdata)
library(terra)
library(dplyr)
library(tidyr)

cat("=== voccTraj Memory Leak Test ===\n\n")

# Setup: Prepare data exactly as in the vignette (lines 208-233)
cat("Setting up data...\n")

# Load data
HadiSST <- terra::rast(system.file("extdata", "HadiSST.tif", package = "VoCCdata"))

# Monthly to annual averages
r <- sumSeries(HadiSST, p = "1960-01/2009-12", yr0 = "1955-01-01",
               l = terra::nlyr(HadiSST),
               fun = function(x) colMeans(x, na.rm = TRUE),
               freqin = "months", freqout = "years")

# Calculate climate velocity
vt <- tempTrend(r, th = 10)
vg <- spatGrad(r, th = 0.0001, projected = FALSE)
gv <- gVoCC(vt, vg)

# Prepare raster layers
vel <- gv[[1]]
ang <- gv[[2]]
mn <- terra::app(r, mean, na.rm = TRUE)

# Generate a velocity layer centered and cropped to study region
x1 <- terra::crop(gv[[1]], terra::ext(-180, 0, -90, 90))
x2 <- terra::crop(gv[[1]], terra::ext(0, 180, -90, 90))
terra::ext(x1) <- c(180, 360, -90, 90)
velc <- terra::merge(x1, x2)

# Crop to the desired extent
velc <- terra::crop(velc, c(90, 180, -32, 33))

# Populate the data frame with cell centroid coordinates
lonlat <- data.frame(terra::xyFromCell(velc, 1:ncell(velc)))
lonlat$vel <- terra::extract(vel, lonlat, ID = FALSE)
lonlat$ang <- terra::extract(ang, lonlat[, 1:2], ID = FALSE)
lonlat$mn <- terra::extract(mn, lonlat[, 1:2], ID = FALSE)
lonlat$lineID <- 1:nrow(lonlat)
lonlat <- drop_na(lonlat)

cat("Data setup complete.\n")
cat("Number of trajectories to calculate:", nrow(lonlat), "\n\n")

# Function to get memory usage in MB
get_memory_mb <- function() {
  gc_info <- gc(reset = TRUE)
  # Sum of used memory across all categories (Ncells * (Vcells size))
  used_mb <- sum(gc_info[, "used"] * c(8, 56)) / (1024^2)
  return(used_mb)
}

# Store results
results <- data.frame(
  iteration = integer(),
  memory_before_mb = numeric(),
  memory_after_mb = numeric(),
  memory_increase_mb = numeric(),
  time_seconds = numeric()
)

cat("Starting memory leak test (10 iterations)...\n")
cat(sprintf("%-10s | %-15s | %-15s | %-18s | %-12s\n",
            "Iteration", "Before (MB)", "After (MB)", "Increase (MB)", "Time (sec)"))
cat(paste(rep("-", 85), collapse = ""), "\n")

# Run voccTraj 10 times
for (i in 1:30) {

  # Force garbage collection and measure memory before
  gc(full = TRUE)
  mem_before <- get_memory_mb()

  # Time the execution
  start_time <- Sys.time()

  # Run voccTraj with same arguments as vignette line 239
  traj <- voccTraj(lonlat, vel, ang, mn, tyr = 10, tstep = 1/4, seed = 23)

  end_time <- Sys.time()
  elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))

  # Measure memory after (before cleanup)
  mem_after <- get_memory_mb()
  mem_increase <- mem_after - mem_before

  # Store results
  results <- rbind(results, data.frame(
    iteration = i,
    memory_before_mb = mem_before,
    memory_after_mb = mem_after,
    memory_increase_mb = mem_increase,
    time_seconds = elapsed
  ))

  # Print progress
  cat(sprintf("%-10d | %-15.2f | %-15.2f | %-18.2f | %-12.2f\n",
              i, mem_before, mem_after, mem_increase, elapsed))

  # Clean up explicitly
  rm(traj)
  gc(full = TRUE)

  # Small delay to let system stabilize
  Sys.sleep(0.5)
}

cat(paste(rep("-", 85), collapse = ""), "\n\n")

# Analysis
cat("=== ANALYSIS ===\n\n")
cat("Memory Usage Summary:\n")
cat(sprintf("  Mean memory before:   %.2f MB (SD: %.2f)\n",
            mean(results$memory_before_mb), sd(results$memory_before_mb)))
cat(sprintf("  Mean memory after:    %.2f MB (SD: %.2f)\n",
            mean(results$memory_after_mb), sd(results$memory_after_mb)))
cat(sprintf("  Mean memory increase: %.2f MB (SD: %.2f)\n",
            mean(results$memory_increase_mb), sd(results$memory_increase_mb)))
cat(sprintf("  Min increase:         %.2f MB\n", min(results$memory_increase_mb)))
cat(sprintf("  Max increase:         %.2f MB\n", max(results$memory_increase_mb)))
cat("\n")

cat("Execution Time Summary:\n")
cat(sprintf("  Mean execution time:  %.2f seconds (SD: %.2f)\n",
            mean(results$time_seconds), sd(results$time_seconds)))
cat(sprintf("  Min time:             %.2f seconds\n", min(results$time_seconds)))
cat(sprintf("  Max time:             %.2f seconds\n\n", max(results$time_seconds)))

# Check for memory leak patterns
# Calculate trend: if memory consistently increases across iterations, there's likely a leak
if (nrow(results) > 2) {
  # Simple linear regression of memory_before vs iteration
  lm_model <- lm(memory_before_mb ~ iteration, data = results)
  slope <- coef(lm_model)[2]
  p_value <- summary(lm_model)$coefficients[2, 4]

  cat("Memory Leak Detection:\n")
  cat(sprintf("  Trend slope:          %.4f MB per iteration\n", slope))
  cat(sprintf("  Statistical p-value:  %.4f\n", p_value))

  if (p_value < 0.05 && slope > 0.5) {
    cat("  ⚠️  WARNING: Significant upward trend detected - possible memory leak!\n")
  } else if (abs(slope) < 0.5) {
    cat("  ✓ No significant memory leak detected - baseline memory stable.\n")
  } else {
    cat("  ? Inconclusive - may need more iterations or longer monitoring.\n")
  }
}

cat("\n=== Test Complete ===\n")
cat("Results saved in 'results' data frame\n")
cat("To view full results: print(results)\n")

# Optionally save results to file
write.csv(results, "data-raw/voccTraj_memory_test_results.csv", row.names = FALSE)
cat("Results also saved to: voccTraj_memory_test_results.csv\n")

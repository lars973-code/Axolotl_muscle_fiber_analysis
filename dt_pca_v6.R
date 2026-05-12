library(tidyverse)
library(sf)
library(FNN)

# =============================================================================
# 1. LOAD CLEANED DATA FROM v6
# =============================================================================

df_clean <- read_csv("phalloidin_clean_with_hull_v6.csv")
optimal_ratios <- read_csv("sample_optimal_ratios_v6.csv")

df_hull <- df_clean %>% filter(within_hull == TRUE)
cat("Samples:", n_distinct(df_hull$Sample), "\n")
cat("Within-hull fibers:", nrow(df_hull), "\n")


# =============================================================================
# 2. REBUILD HULLS
# =============================================================================

cat("\nRebuilding hulls...\n")

hull_list <- map2(
  optimal_ratios$Sample,
  optimal_ratios$optimal_ratio,
  function(s, ratio) {
    xy <- df_hull %>% filter(Sample == s) %>%
      select(X, Y) %>% as.matrix()
    tryCatch(
      st_concave_hull(st_multipoint(xy), ratio = ratio),
      error = function(e) {
        warning("Hull failed for ", s, " -- falling back to convex hull.")
        st_convex_hull(st_multipoint(xy))
      }
    )
  }
)
names(hull_list) <- optimal_ratios$Sample
cat("Hulls rebuilt.\n")


# =============================================================================
# 3. DISTANCE TRANSFORM MATRIX + PCA
#
#    For each sample:
#      a) Extract hull boundary vertices
#      b) Lay a GRID_N x GRID_N regular grid over the bounding box
#      c) Determine which grid points fall inside the hull
#      d) For interior points: compute mean distance to all boundary vertices
#         For exterior points: set to 0
#      e) Treat the 2D grid as a weighted surface -- each interior grid cell
#         has coordinates (gx, gy) and a weight = mean_distance_value.
#         Run weighted PCA on (gx, gy) with weights = DT values.
#         This finds the axes along which the distance-weighted interior
#         mass is most spread out.
#
#    The mean distance is computed in chunks to manage memory.
# =============================================================================

GRID_N     <- 200       # grid resolution (200 x 200)
CHUNK_SIZE <- 5000      # rows per chunk for distance computation
BOUNDARY_RESAMPLE <- 300 # resample boundary to this many evenly-spaced points

cat("\nGrid resolution:", GRID_N, "x", GRID_N, "\n")
cat("Boundary resampled to:", BOUNDARY_RESAMPLE, "points per sample\n")

# --- Helper: resample polygon boundary to n evenly spaced points ---
resample_boundary <- function(hull_geom, n_pts) {
  # Convert hull to linestring for even sampling along perimeter
  ring <- st_cast(hull_geom, "MULTILINESTRING")
  total_len <- as.numeric(st_length(ring))
  if (total_len == 0 || is.na(total_len)) {
    # Fallback: just return raw coordinates
    coords <- st_coordinates(hull_geom)
    return(coords[, c("X", "Y")])
  }
  # Sample at evenly spaced distances along the boundary
  fracs <- seq(0, 1, length.out = n_pts + 1)[-(n_pts + 1)]
  pts <- st_line_sample(ring, sample = fracs)
  coords <- st_coordinates(pts)
  coords[, c("X", "Y")]
}

# --- Helper: compute DT-weighted PCA for one sample ---
compute_dt_pca <- function(s, hull_geom) {
  cat("  Processing", s, "...")
  
  # (a) Resample boundary to evenly spaced points
  bnd_xy <- tryCatch(
    resample_boundary(hull_geom, BOUNDARY_RESAMPLE),
    error = function(e) {
      coords <- st_coordinates(hull_geom)
      coords[, c("X", "Y")]
    }
  )
  if (is.null(bnd_xy) || nrow(bnd_xy) < 3) {
    warning("Skipping ", s, " -- insufficient boundary points.")
    return(NULL)
  }
  
  # (b) Build regular grid over bounding box (with small padding)
  bbox   <- st_bbox(hull_geom)
  pad_x  <- 0.02 * (bbox["xmax"] - bbox["xmin"])
  pad_y  <- 0.02 * (bbox["ymax"] - bbox["ymin"])
  gx_seq <- seq(bbox["xmin"] - pad_x, bbox["xmax"] + pad_x, length.out = GRID_N)
  gy_seq <- seq(bbox["ymin"] - pad_y, bbox["ymax"] + pad_y, length.out = GRID_N)
  grid   <- expand.grid(gx = gx_seq, gy = gy_seq)
  
  # (c) Determine which grid points are inside the hull
  grid_sf <- st_as_sf(grid, coords = c("gx", "gy"), crs = NA)
  inside  <- as.logical(st_intersects(grid_sf, hull_geom, sparse = FALSE)[, 1])
  grid$inside <- inside
  
  n_inside <- sum(inside)
  if (n_inside < 10) {
    warning("Skipping ", s, " -- only ", n_inside, " interior grid points.")
    return(NULL)
  }
  
  # (d) Compute mean distance to boundary for interior points (chunked)
  interior_xy <- as.matrix(grid[inside, c("gx", "gy")])
  bnd_mat     <- as.matrix(bnd_xy)
  
  mean_dists <- numeric(n_inside)
  n_chunks   <- ceiling(n_inside / CHUNK_SIZE)
  
  for (ch in seq_len(n_chunks)) {
    idx_start <- (ch - 1) * CHUNK_SIZE + 1
    idx_end   <- min(ch * CHUNK_SIZE, n_inside)
    chunk_xy  <- interior_xy[idx_start:idx_end, , drop = FALSE]
    
    # Distance matrix: each row is a grid point, each col is a boundary point
    # Using vectorized outer distance computation
    dx <- outer(chunk_xy[, 1], bnd_mat[, 1], "-")
    dy <- outer(chunk_xy[, 2], bnd_mat[, 2], "-")
    d  <- sqrt(dx^2 + dy^2)
    mean_dists[idx_start:idx_end] <- rowMeans(d)
  }
  
  # (e) Weighted PCA on interior grid coordinates
  #     Weights = mean distance values (higher weight for deep interior)
  weights <- mean_dists / sum(mean_dists)  # normalize to sum to 1
  
  # Weighted centroid
  cx <- sum(interior_xy[, 1] * weights)
  cy <- sum(interior_xy[, 2] * weights)
  
  # Weighted covariance matrix
  dx_c <- interior_xy[, 1] - cx
  dy_c <- interior_xy[, 2] - cy
  cov_xx <- sum(weights * dx_c^2)
  cov_yy <- sum(weights * dy_c^2)
  cov_xy <- sum(weights * dx_c * dy_c)
  cov_mat <- matrix(c(cov_xx, cov_xy, cov_xy, cov_yy), nrow = 2)
  
  # Eigen decomposition
  eig <- eigen(cov_mat, symmetric = TRUE)
  ev1 <- eig$values[1]
  ev2 <- eig$values[2]
  
  # PC directions (eigenvectors)
  pc1_dx <- eig$vectors[1, 1]
  pc1_dy <- eig$vectors[2, 1]
  pc2_dx <- eig$vectors[1, 2]
  pc2_dy <- eig$vectors[2, 2]
  
  pc1_angle <- atan2(pc1_dy, pc1_dx) * 180 / pi
  
  cat(" done. (", n_inside, "interior pts, angle=", round(pc1_angle, 1), "°)\n")
  
  tibble(
    Sample     = s,
    centroid_x = cx,
    centroid_y = cy,
    pc1_dx     = pc1_dx,
    pc1_dy     = pc1_dy,
    pc2_dx     = pc2_dx,
    pc2_dy     = pc2_dy,
    ev1        = ev1,
    ev2        = ev2,
    pct_var1   = round(100 * ev1 / (ev1 + ev2), 1),
    pc1_angle  = round(pc1_angle, 2),
    n_interior = n_inside
  )
}


# =============================================================================
# 4. RUN DT-PCA FOR ALL SAMPLES
# =============================================================================

cat("\nComputing distance-transform PCA...\n")

dt_pca_results <- map_dfr(names(hull_list), function(s) {
  compute_dt_pca(s, hull_list[[s]])
})

cat("\n--- DT-PCA results ---\n")
print(dt_pca_results %>% select(Sample, pc1_angle, pct_var1, n_interior), n = Inf)

write_csv(dt_pca_results, "dt_pca_results_v6.csv")
cat("Saved -> dt_pca_results_v6.csv\n")


# =============================================================================
# 5. COMPARE BOUNDARY PCA vs DT-PCA ANGLES
#    Load the boundary PCA results from the previous script for comparison.
# =============================================================================

if (file.exists("boundary_pca_results_v6.csv")) {
  bnd_pca <- read_csv("boundary_pca_results_v6.csv")
  
  comparison <- dt_pca_results %>%
    select(Sample, dt_angle = pc1_angle, dt_var1 = pct_var1) %>%
    left_join(
      bnd_pca %>% select(Sample, bnd_angle = pc1_angle, bnd_var1 = pct_var1),
      by = "Sample"
    ) %>%
    mutate(angle_diff = dt_angle - bnd_angle)
  
  cat("\n--- Boundary PCA vs DT-PCA comparison ---\n")
  print(comparison, n = Inf)
  
  cat("\nAngle difference summary:\n")
  cat("  Mean abs diff: ", round(mean(abs(comparison$angle_diff), na.rm = TRUE), 2), "°\n")
  cat("  Max abs diff:  ", round(max(abs(comparison$angle_diff), na.rm = TRUE), 2), "°\n")
  
  write_csv(comparison, "pca_method_comparison_v6.csv")
  cat("Saved -> pca_method_comparison_v6.csv\n")
}


# =============================================================================
# 6. BUILD AXIS LINE SEGMENTS FOR PLOTTING
# =============================================================================

AXIS_SCALE <- 0.4

axis_segments <- dt_pca_results %>%
  mutate(
    len1 = AXIS_SCALE * sqrt(ev1),
    len2 = AXIS_SCALE * sqrt(ev2),
    
    pc1_x_start = centroid_x - pc1_dx * len1,
    pc1_y_start = centroid_y - pc1_dy * len1,
    pc1_x_end   = centroid_x + pc1_dx * len1,
    pc1_y_end   = centroid_y + pc1_dy * len1,
    
    pc2_x_start = centroid_x - pc2_dx * len2,
    pc2_y_start = centroid_y - pc2_dy * len2,
    pc2_x_end   = centroid_x + pc2_dx * len2,
    pc2_y_end   = centroid_y + pc2_dy * len2
  )


# =============================================================================
# 7. HULL POLYGON COORDINATES
# =============================================================================

hull_polygons <- imap_dfr(hull_list, function(geom, s) {
  tryCatch({
    coords        <- as.data.frame(st_coordinates(geom))
    coords$Sample <- s
    coords
  }, error = function(e) NULL)
})


# =============================================================================
# 8. PLOT: DT-PCA AXES OVER HULL + FIBERS
# =============================================================================

cat("\nRendering DT-PCA axis overlay...\n")

df_thin <- df_hull %>%
  group_by(Sample) %>%
  slice_sample(n = 500) %>%
  ungroup()

p_dt_pca <- ggplot() +
  geom_point(data = df_thin,
             aes(x = X, y = Y),
             colour = "#88ccaa", size = 0.1, alpha = 0.3) +
  
  geom_polygon(data = hull_polygons,
               aes(x = X, y = Y),
               fill = NA, colour = "dodgerblue", linewidth = 0.4) +
  
  # PC1 axis (red)
  geom_segment(data = axis_segments,
               aes(x = pc1_x_start, y = pc1_y_start,
                   xend = pc1_x_end, yend = pc1_y_end),
               colour = "#ff4d4d", linewidth = 0.9,
               arrow = arrow(length = unit(0.08, "inches"),
                             ends = "both", type = "open")) +
  
  # PC2 axis (orange)
  geom_segment(data = axis_segments,
               aes(x = pc2_x_start, y = pc2_y_start,
                   xend = pc2_x_end, yend = pc2_y_end),
               colour = "#f4a261", linewidth = 0.7,
               arrow = arrow(length = unit(0.06, "inches"),
                             ends = "both", type = "open")) +
  
  # Centroid
  geom_point(data = dt_pca_results,
             aes(x = centroid_x, y = centroid_y),
             colour = "white", fill = "#333333",
             shape = 21, size = 1.5, stroke = 0.4) +
  
  # Angle label
  geom_text(data = dt_pca_results,
            aes(x = -Inf, y = Inf,
                label = paste0(pc1_angle, "° (", pct_var1, "%)")),
            hjust = -0.05, vjust = 1.5,
            colour = "#ff4d4d", size = 1.8) +
  
  facet_wrap(~Sample, scales = "free", ncol = 6) +
  theme_minimal(base_size = 7) +
  labs(
    title    = "Distance-transform PCA (200x200 grid) -- Old samples, outliers removed",
    subtitle = "Red = PC1  |  Orange = PC2  |  Label = PC1 angle & % var  |  Weights = mean dist to boundary",
    x = "X (pixels)", y = "Y (pixels)"
  ) +
  theme(
    strip.text      = element_text(size = 6),
    axis.text       = element_text(size = 4),
    panel.spacing   = unit(0.3, "lines"),
    legend.position = "none"
  )

png("dt_pca_overlay_v6.png", width = 5400, height = 4200, res = 200)
print(p_dt_pca)
dev.off()
cat("Saved -> dt_pca_overlay_v6.png\n")


# =============================================================================
# 9. ANGLE DISTRIBUTION SUMMARY
# =============================================================================

cat("\n--- DT-PCA PC1 angle summary ---\n")
cat("Mean:  ", round(mean(dt_pca_results$pc1_angle), 2), "°\n")
cat("SD:    ", round(sd(dt_pca_results$pc1_angle), 2), "°\n")
cat("Range: ", round(min(dt_pca_results$pc1_angle), 2), "° to",
    round(max(dt_pca_results$pc1_angle), 2), "°\n")

cat("\nAll outputs:\n")
cat("  dt_pca_results_v6.csv          -- per-sample DT-PCA results\n")
cat("  pca_method_comparison_v6.csv   -- boundary PCA vs DT-PCA angle comparison\n")
cat("  dt_pca_overlay_v6.png          -- hull overlay with DT-PCA axes\n")
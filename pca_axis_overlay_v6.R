library(tidyverse)
library(sf)

# =============================================================================
# 1. LOAD CLEANED DATA FROM v6
#    Uses the cleaned phalloidin data (outliers removed) and optimal ratios.
# =============================================================================

df_clean <- read_csv("phalloidin_clean_with_hull_v6.csv")
optimal_ratios <- read_csv("sample_optimal_ratios_v6.csv")

cat("Cleaned fibers loaded:", nrow(df_clean), "\n")
cat("Samples:", n_distinct(df_clean$Sample), "\n")

# Keep only within-hull fibers for PCA (these define the tissue boundary)
df_hull <- df_clean %>% filter(within_hull == TRUE)
cat("Within-hull fibers:", nrow(df_hull), "\n")


# =============================================================================
# 2. REBUILD HULLS (needed for boundary coordinates)
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
# 3. EXTRACT HULL BOUNDARY COORDINATES & RUN PCA
#    PCA is performed on the boundary vertices of each sample's concave hull.
#    This captures the dominant orientation of the tissue cross-section shape.
#
#    Outputs per sample:
#      - centroid (mean X, mean Y of boundary)
#      - PC1 and PC2 direction vectors (unit vectors)
#      - eigenvalues (variance explained along each axis)
#      - PC1 slope angle (degrees from horizontal)
# =============================================================================

cat("\nRunning PCA on hull boundary coordinates...\n")

pca_results <- map_dfr(names(hull_list), function(s) {
  # Extract boundary coordinates from hull geometry
  coords <- tryCatch({
    as.data.frame(st_coordinates(hull_list[[s]]))
  }, error = function(e) {
    warning("Could not extract coordinates for ", s)
    return(NULL)
  })
  
  if (is.null(coords) || nrow(coords) < 3) {
    warning("Too few boundary points for PCA on ", s)
    return(NULL)
  }
  
  xy <- as.matrix(coords[, c("X", "Y")])
  
  # Centroid
  cx <- mean(xy[, 1])
  cy <- mean(xy[, 2])
  
  # PCA
  pca <- prcomp(xy, center = TRUE, scale. = FALSE)
  
  # Eigenvalues (variances)
  ev1 <- pca$sdev[1]^2
  ev2 <- pca$sdev[2]^2
  
  # Unit direction vectors (rotation matrix columns)
  pc1_dx <- pca$rotation[1, 1]
  pc1_dy <- pca$rotation[2, 1]
  pc2_dx <- pca$rotation[1, 2]
  pc2_dy <- pca$rotation[2, 2]
  
  # Angle of PC1 from horizontal (degrees)
  pc1_angle <- atan2(pc1_dy, pc1_dx) * 180 / pi
  
  tibble(
    Sample    = s,
    centroid_x = cx,
    centroid_y = cy,
    pc1_dx    = pc1_dx,
    pc1_dy    = pc1_dy,
    pc2_dx    = pc2_dx,
    pc2_dy    = pc2_dy,
    ev1       = ev1,
    ev2       = ev2,
    pct_var1  = round(100 * ev1 / (ev1 + ev2), 1),
    pc1_angle = round(pc1_angle, 2)
  )
})

cat("\nPCA summary per sample:\n")
print(pca_results %>% select(Sample, pc1_angle, pct_var1, ev1, ev2), n = Inf)

write_csv(pca_results, "boundary_pca_results_v6.csv")
cat("Saved -> boundary_pca_results_v6.csv\n")


# =============================================================================
# 4. BUILD AXIS LINE SEGMENTS FOR PLOTTING
#    Scale each PC axis proportional to its eigenvalue (sqrt for length units)
#    so the line lengths visually represent variance along each direction.
# =============================================================================

AXIS_SCALE <- 0.4   # fraction of sqrt(eigenvalue) to use as half-length

axis_segments <- pca_results %>%
  mutate(
    # Half-lengths scaled by sqrt of eigenvalue
    len1 = AXIS_SCALE * sqrt(ev1),
    len2 = AXIS_SCALE * sqrt(ev2),
    
    # PC1 endpoints
    pc1_x_start = centroid_x - pc1_dx * len1,
    pc1_y_start = centroid_y - pc1_dy * len1,
    pc1_x_end   = centroid_x + pc1_dx * len1,
    pc1_y_end   = centroid_y + pc1_dy * len1,
    
    # PC2 endpoints
    pc2_x_start = centroid_x - pc2_dx * len2,
    pc2_y_start = centroid_y - pc2_dy * len2,
    pc2_x_end   = centroid_x + pc2_dx * len2,
    pc2_y_end   = centroid_y + pc2_dy * len2
  )


# =============================================================================
# 5. HULL POLYGON COORDINATES FOR OVERLAY
# =============================================================================

hull_polygons <- imap_dfr(hull_list, function(geom, s) {
  tryCatch({
    coords        <- as.data.frame(st_coordinates(geom))
    coords$Sample <- s
    coords
  }, error = function(e) NULL)
})


# =============================================================================
# 6. PLOT: FIBER SCATTER + HULL + PCA AXES
# =============================================================================

cat("\nRendering PCA axis overlay plot...\n")

# Thin fibers for rendering speed
df_thin <- df_hull %>%
  group_by(Sample) %>%
  slice_sample(n = 500) %>%
  ungroup()

p_pca <- ggplot() +
  # Fiber scatter
  geom_point(data = df_thin,
             aes(x = X, y = Y),
             colour = "#88ccaa", size = 0.1, alpha = 0.3) +
  
  # Hull boundary
  geom_polygon(data = hull_polygons,
               aes(x = X, y = Y),
               fill = NA, colour = "dodgerblue", linewidth = 0.4) +
  
  # PC1 axis (red, thicker -- dominant axis)
  geom_segment(data = axis_segments,
               aes(x = pc1_x_start, y = pc1_y_start,
                   xend = pc1_x_end, yend = pc1_y_end),
               colour = "#ff4d4d", linewidth = 0.9,
               arrow = arrow(length = unit(0.08, "inches"),
                             ends = "both", type = "open")) +
  
  # PC2 axis (orange, thinner -- minor axis)
  geom_segment(data = axis_segments,
               aes(x = pc2_x_start, y = pc2_y_start,
                   xend = pc2_x_end, yend = pc2_y_end),
               colour = "#f4a261", linewidth = 0.7,
               arrow = arrow(length = unit(0.06, "inches"),
                             ends = "both", type = "open")) +
  
  # Centroid marker
  geom_point(data = pca_results,
             aes(x = centroid_x, y = centroid_y),
             colour = "white", fill = "#333333",
             shape = 21, size = 1.5, stroke = 0.4) +
  
  # PC1 angle label
  geom_text(data = pca_results,
            aes(x = -Inf, y = Inf,
                label = paste0(pc1_angle, "°  (",  pct_var1, "%)")),
            hjust = -0.05, vjust = 1.5,
            colour = "#ff4d4d", size = 1.8) +
  
  facet_wrap(~Sample, scales = "free", ncol = 6) +
  theme_minimal(base_size = 7) +
  labs(
    title    = "PCA on hull boundary -- Old samples, outliers removed",
    subtitle = "Red = PC1 (dominant axis)  |  Orange = PC2  |  Label = PC1 angle & % variance",
    x = "X (pixels)", y = "Y (pixels)"
  ) +
  theme(
    strip.text      = element_text(size = 6),
    axis.text       = element_text(size = 4),
    panel.spacing   = unit(0.3, "lines"),
    legend.position = "none"
  )

png("pca_axis_overlay_v6.png", width = 5400, height = 4200, res = 200)
print(p_pca)
dev.off()
cat("Saved -> pca_axis_overlay_v6.png\n")


# =============================================================================
# 7. SUMMARY: ANGLE DISTRIBUTION ACROSS SAMPLES
#    Quick check of how much rotation variability exists.
# =============================================================================

cat("\n--- PC1 angle summary across samples ---\n")
cat("Mean:  ", round(mean(pca_results$pc1_angle), 2), "°\n")
cat("SD:    ", round(sd(pca_results$pc1_angle), 2), "°\n")
cat("Range: ", round(min(pca_results$pc1_angle), 2), "° to",
    round(max(pca_results$pc1_angle), 2), "°\n")

cat("\nAll outputs:\n")
cat("  boundary_pca_results_v6.csv  -- per-sample PCA: centroid, axes, angles, variance\n")
cat("  pca_axis_overlay_v6.png      -- hull + fiber scatter with PC1/PC2 axes overlaid\n")
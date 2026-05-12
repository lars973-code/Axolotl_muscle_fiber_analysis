library(tidyverse)
library(sf)

# =============================================================================
# 1. LOAD DATA
# =============================================================================

df_clean     <- read_csv("phalloidin_clean_with_hull_v6.csv")
dt_pca       <- read_csv("dt_pca_results_v6.csv")
optimal_ratios <- read_csv("sample_optimal_ratios_v6.csv")

df_hull <- df_clean %>% filter(within_hull == TRUE)
cat("Samples:", n_distinct(df_hull$Sample), "\n")
cat("Within-hull fibers:", nrow(df_hull), "\n")


# =============================================================================
# 2. NORMALIZE PC1 ANGLES TO 0–180°
#    PCA axes are bidirectional, so angles that differ by 180° are equivalent.
#    We map everything into [0, 180) for consistency, then compute the
#    rotation needed to bring PC1 to horizontal (0°).
# =============================================================================

dt_pca <- dt_pca %>%
  mutate(
    # Normalize to [0, 180): if angle is negative, add 180
    pc1_angle_norm = case_when(
      pc1_angle < 0   ~ pc1_angle + 180,
      pc1_angle >= 180 ~ pc1_angle - 180,
      TRUE             ~ pc1_angle
    ),
    # Rotation angle to bring PC1 to horizontal (0°)
    # Negative because we rotate the coordinates in the opposite direction
    rotation_deg = -pc1_angle_norm,
    rotation_rad = rotation_deg * pi / 180
  )

cat("\n--- Normalized PC1 angles ---\n")
print(dt_pca %>% select(Sample, pc1_angle, pc1_angle_norm, rotation_deg), n = Inf)

cat("\nNormalized angle summary:\n")
cat("  Mean: ", round(mean(dt_pca$pc1_angle_norm), 2), "°\n")
cat("  SD:   ", round(sd(dt_pca$pc1_angle_norm), 2), "°\n")
cat("  Range:", round(min(dt_pca$pc1_angle_norm), 2), "° to",
    round(max(dt_pca$pc1_angle_norm), 2), "°\n")


# =============================================================================
# 3. ROTATE & CENTER FIBER COORDINATES
#    For each sample:
#      a) Translate so DT-PCA centroid is at origin
#      b) Rotate by -pc1_angle_norm so PC1 aligns to horizontal
#
#    This produces new columns: X_aligned, Y_aligned
# =============================================================================

cat("\nAligning fiber coordinates...\n")

df_aligned <- df_hull %>%
  left_join(
    dt_pca %>% select(Sample, centroid_x, centroid_y, rotation_rad),
    by = "Sample"
  ) %>%
  mutate(
    # Center on DT-PCA centroid
    X_c = X - centroid_x,
    Y_c = Y - centroid_y,
    
    # Rotate to align PC1 with horizontal
    X_aligned = X_c * cos(rotation_rad) - Y_c * sin(rotation_rad),
    Y_aligned = X_c * sin(rotation_rad) + Y_c * cos(rotation_rad)
  )

cat("Alignment complete.\n")
cat("Aligned fibers:", nrow(df_aligned), "\n")


# =============================================================================
# 4. REBUILD & TRANSFORM HULLS
#    Apply the same translation + rotation to hull boundaries for plotting.
# =============================================================================

cat("\nRebuilding and transforming hulls...\n")

hull_list <- map2(
  optimal_ratios$Sample,
  optimal_ratios$optimal_ratio,
  function(s, ratio) {
    xy <- df_hull %>% filter(Sample == s) %>%
      select(X, Y) %>% as.matrix()
    tryCatch(
      st_concave_hull(st_multipoint(xy), ratio = ratio),
      error = function(e) st_convex_hull(st_multipoint(xy))
    )
  }
)
names(hull_list) <- optimal_ratios$Sample

# Extract hull polygon coords and apply the same transform
hull_polygons_aligned <- imap_dfr(hull_list, function(geom, s) {
  tryCatch({
    coords <- as.data.frame(st_coordinates(geom))
    
    params <- dt_pca %>% filter(Sample == s)
    cx   <- params$centroid_x
    cy   <- params$centroid_y
    rad  <- params$rotation_rad
    
    coords %>%
      mutate(
        X_c = X - cx,
        Y_c = Y - cy,
        X_aligned = X_c * cos(rad) - Y_c * sin(rad),
        Y_aligned = X_c * sin(rad) + Y_c * cos(rad),
        Sample = s
      )
  }, error = function(e) NULL)
})

cat("Hull transformation done.\n")


# =============================================================================
# 5. VERIFY ALIGNMENT -- PCA ON ALIGNED COORDINATES
#    After rotation, PC1 should be near 0° (horizontal) for all samples.
# =============================================================================

cat("\nVerifying alignment with PCA on aligned coordinates...\n")

verify_pca <- df_aligned %>%
  group_by(Sample) %>%
  summarise(
    n_fibers = n(),
    # Quick PCA on aligned coords
    pca_check = list({
      xy <- cbind(X_aligned, Y_aligned)
      pca <- prcomp(xy, center = TRUE, scale. = FALSE)
      ev <- pca$sdev^2
      angle <- atan2(pca$rotation[2, 1], pca$rotation[1, 1]) * 180 / pi
      # Normalize to [0, 180)
      if (angle < 0) angle <- angle + 180
      if (angle >= 180) angle <- angle - 180
      tibble(post_angle = round(angle, 2),
             post_var1  = round(100 * ev[1] / sum(ev), 1))
    }),
    .groups = "drop"
  ) %>%
  unnest(pca_check)

cat("\n--- Post-alignment PC1 angles (should be near 0° or 180°) ---\n")
print(verify_pca %>% select(Sample, post_angle, post_var1), n = Inf)

cat("\nPost-alignment angle summary:\n")
cat("  Mean: ", round(mean(verify_pca$post_angle), 2), "°\n")
cat("  SD:   ", round(sd(verify_pca$post_angle), 2), "°\n")


# =============================================================================
# 6. DIAGNOSTIC PLOTS
# =============================================================================

# --- 6a. BEFORE vs AFTER: side-by-side comparison ---
cat("\nRendering before/after comparison...\n")

# Thin for plotting
df_before_thin <- df_hull %>%
  group_by(Sample) %>%
  slice_sample(n = 500) %>%
  ungroup() %>%
  mutate(stage = "Before alignment")

df_after_thin <- df_aligned %>%
  group_by(Sample) %>%
  slice_sample(n = 500) %>%
  ungroup() %>%
  select(Sample, X = X_aligned, Y = Y_aligned) %>%
  mutate(stage = "After alignment")

# Hull polygons before
hull_polygons_before <- imap_dfr(hull_list, function(geom, s) {
  tryCatch({
    coords <- as.data.frame(st_coordinates(geom))
    coords$Sample <- s
    coords$stage  <- "Before alignment"
    coords
  }, error = function(e) NULL)
})

hull_polygons_after <- hull_polygons_aligned %>%
  select(X = X_aligned, Y = Y_aligned, Sample) %>%
  mutate(stage = "After alignment")


# --- 6b. Aligned overlay (all samples) ---
cat("Rendering aligned hull overlay...\n")

p_aligned <- ggplot() +
  geom_point(data = df_after_thin,
             aes(x = X, y = Y),
             colour = "#88ccaa", size = 0.1, alpha = 0.3) +
  geom_polygon(data = hull_polygons_after,
               aes(x = X, y = Y),
               fill = NA, colour = "dodgerblue", linewidth = 0.4) +
  # Horizontal reference line through origin (should align with PC1)
  geom_hline(yintercept = 0, colour = "#ff4d4d",
             linewidth = 0.3, linetype = "dashed") +
  geom_vline(xintercept = 0, colour = "#f4a261",
             linewidth = 0.3, linetype = "dashed") +
  facet_wrap(~Sample, scales = "free", ncol = 6) +
  theme_minimal(base_size = 7) +
  labs(
    title    = "Aligned samples -- PC1 rotated to horizontal, centered at origin",
    subtitle = "Red dashed = horizontal (PC1 target)  |  Orange dashed = vertical (PC2)",
    x = "X aligned (pixels)", y = "Y aligned (pixels)"
  ) +
  theme(
    strip.text      = element_text(size = 6),
    axis.text       = element_text(size = 4),
    panel.spacing   = unit(0.3, "lines"),
    legend.position = "none"
  )

png("aligned_overlay_v6.png", width = 5400, height = 4200, res = 200)
print(p_aligned)
dev.off()
cat("Saved -> aligned_overlay_v6.png\n")


# --- 6c. All samples superimposed (normalized scale) ---
cat("Rendering superimposed overlay...\n")

# Normalize each sample's coordinates to unit scale for superposition
df_superimposed <- df_aligned %>%
  group_by(Sample) %>%
  mutate(
    # Scale to [-1, 1] range based on max absolute extent
    max_extent = max(abs(c(X_aligned, Y_aligned))),
    X_norm = X_aligned / max_extent,
    Y_norm = Y_aligned / max_extent
  ) %>%
  ungroup()

hull_superimposed <- hull_polygons_aligned %>%
  left_join(
    df_superimposed %>%
      group_by(Sample) %>%
      summarise(max_extent = first(max_extent), .groups = "drop"),
    by = "Sample"
  ) %>%
  mutate(
    X_norm = X_aligned / max_extent,
    Y_norm = Y_aligned / max_extent
  )

# Thin for overlay
df_super_thin <- df_superimposed %>%
  group_by(Sample) %>%
  slice_sample(n = 200) %>%
  ungroup()

p_super <- ggplot() +
  geom_polygon(data = hull_superimposed,
               aes(x = X_norm, y = Y_norm, group = Sample),
               fill = NA, colour = "dodgerblue", linewidth = 0.15, alpha = 0.4) +
  geom_hline(yintercept = 0, colour = "#ff4d4d",
             linewidth = 0.4, linetype = "dashed") +
  geom_vline(xintercept = 0, colour = "#f4a261",
             linewidth = 0.4, linetype = "dashed") +
  theme_minimal(base_size = 10) +
  labs(
    title    = "All samples superimposed -- aligned & normalized to unit scale",
    subtitle = "Each blue outline = one sample hull  |  Dashed = aligned axes",
    x = "X (normalized)", y = "Y (normalized)"
  ) +
  theme(legend.position = "none")

png("superimposed_hulls_v6.png", width = 2400, height = 2400, res = 200)
print(p_super)
dev.off()
cat("Saved -> superimposed_hulls_v6.png\n")


# =============================================================================
# 7. EXPORT
# =============================================================================

# Export all columns -- avoids hardcoding column names that may vary
write_csv(df_aligned, "phalloidin_aligned_v6.csv")

write_csv(dt_pca %>%
            select(Sample, centroid_x, centroid_y,
                   pc1_angle, pc1_angle_norm, rotation_deg,
                   pct_var1, ev1, ev2),
          "alignment_parameters_v6.csv")

cat("\n=== Alignment outputs ===\n")
cat("  phalloidin_aligned_v6.csv    -- fiber coords with X_aligned, Y_aligned\n")
cat("  alignment_parameters_v6.csv  -- per-sample rotation & centroid params\n")
cat("  aligned_overlay_v6.png       -- per-sample aligned hull overlay\n")
cat("  superimposed_hulls_v6.png    -- all hulls superimposed at unit scale\n")
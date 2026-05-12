library(tidyverse)
library(sf)

# =============================================================================
# SIGN DISAMBIGUATION -- "more mass to the right"
#
# Run this AFTER alignment_v6.R has completed.
# It loads the aligned data, checks whether each sample has more fiber mass
# on the positive-X (right) or negative-X (left) side, and flips (180°
# rotation = negate both X and Y) any sample where mass is on the left.
# =============================================================================

# --- 1. Load aligned data ---
df_aligned <- read_csv("phalloidin_aligned_v6.csv")
optimal_ratios <- read_csv("sample_optimal_ratios_v6.csv")
dt_pca <- read_csv("dt_pca_results_v6.csv") %>%
  mutate(
    pc1_angle_norm = if_else(pc1_angle < 0, pc1_angle + 180, pc1_angle),
    rotation_deg   = -pc1_angle_norm,
    rotation_rad   = rotation_deg * pi / 180
  )

cat("Loaded", nrow(df_aligned), "aligned fibers across",
    n_distinct(df_aligned$Sample), "samples\n")


# --- 2. Compute mass asymmetry per sample ---
#     Use mean of X_aligned: if negative, more mass is to the left.
#     Also use a weighted version counting fibers on each side.

mass_check <- df_aligned %>%
  group_by(Sample) %>%
  summarise(
    mean_x     = mean(X_aligned),
    n_right    = sum(X_aligned > 0),
    n_left     = sum(X_aligned <= 0),
    pct_right  = round(100 * n_right / n(), 1),
    needs_flip = mean_x < 0,
    .groups    = "drop"
  )

cat("\n--- Mass asymmetry check ---\n")
print(mass_check %>% select(Sample, mean_x, pct_right, needs_flip), n = Inf)

n_flip <- sum(mass_check$needs_flip)
cat("\nSamples to flip:", n_flip, "of", nrow(mass_check), "\n")


# --- 3. Apply flip (180° rotation = negate X and Y) ---
flip_samples <- mass_check %>% filter(needs_flip) %>% pull(Sample)

df_aligned <- df_aligned %>%
  mutate(
    X_aligned = if_else(Sample %in% flip_samples, -X_aligned, X_aligned),
    Y_aligned = if_else(Sample %in% flip_samples, -Y_aligned, Y_aligned)
  )

cat("Flip applied.\n")


# --- 4. Verify: all samples should now have mean_x > 0 ---
verify <- df_aligned %>%
  group_by(Sample) %>%
  summarise(mean_x = round(mean(X_aligned), 1),
            .groups = "drop")

cat("\n--- Post-flip mean X (should all be positive or near zero) ---\n")
print(verify, n = Inf)


# --- 5. Rebuild aligned hulls for plotting ---
cat("\nRebuilding hulls for plotting...\n")

df_hull_orig <- read_csv("phalloidin_clean_with_hull_v6.csv") %>%
  filter(within_hull == TRUE)

hull_list <- map2(
  optimal_ratios$Sample,
  optimal_ratios$optimal_ratio,
  function(s, ratio) {
    xy <- df_hull_orig %>% filter(Sample == s) %>%
      select(X, Y) %>% as.matrix()
    tryCatch(
      st_concave_hull(st_multipoint(xy), ratio = ratio),
      error = function(e) st_convex_hull(st_multipoint(xy))
    )
  }
)
names(hull_list) <- optimal_ratios$Sample

# Transform hull coordinates using original alignment params + flip
hull_polygons_aligned <- imap_dfr(hull_list, function(geom, s) {
  coords <- as.data.frame(st_coordinates(geom))
  params <- dt_pca %>% filter(Sample == s)
  cx  <- as.numeric(params$centroid_x[1])
  cy  <- as.numeric(params$centroid_y[1])
  rad <- as.numeric(params$rotation_deg[1]) * pi / 180
  
  flip <- s %in% flip_samples
  flip_sign <- if (flip) -1 else 1
  
  X_c <- coords$X - cx
  Y_c <- coords$Y - cy
  
  data.frame(
    X_aligned = flip_sign * (X_c * cos(rad) - Y_c * sin(rad)),
    Y_aligned = flip_sign * (X_c * sin(rad) + Y_c * cos(rad)),
    Sample    = s
  )
})


# --- 6. Plot: aligned overlay with sign fixed ---
cat("Rendering sign-fixed aligned overlay...\n")

df_thin <- df_aligned %>%
  group_by(Sample) %>%
  slice_sample(n = 500) %>%
  ungroup()

p_aligned <- ggplot() +
  geom_point(data = df_thin,
             aes(x = X_aligned, y = Y_aligned),
             colour = "#88ccaa", size = 0.1, alpha = 0.3) +
  geom_polygon(data = hull_polygons_aligned,
               aes(x = X_aligned, y = Y_aligned),
               fill = NA, colour = "dodgerblue", linewidth = 0.4) +
  geom_hline(yintercept = 0, colour = "#ff4d4d",
             linewidth = 0.3, linetype = "dashed") +
  geom_vline(xintercept = 0, colour = "#f4a261",
             linewidth = 0.3, linetype = "dashed") +
  facet_wrap(~Sample, scales = "free", ncol = 6) +
  theme_minimal(base_size = 7) +
  labs(
    title    = "Aligned samples -- sign-fixed (more mass to the right)",
    subtitle = "Red = horizontal (PC1)  |  Orange = vertical (PC2)  |  Mass biased rightward",
    x = "X aligned (pixels)", y = "Y aligned (pixels)"
  ) +
  theme(
    strip.text      = element_text(size = 6),
    axis.text       = element_text(size = 4),
    panel.spacing   = unit(0.3, "lines"),
    legend.position = "none"
  )

png("aligned_signfix_overlay_v6.png", width = 5400, height = 4200, res = 200)
print(p_aligned)
dev.off()
cat("Saved -> aligned_signfix_overlay_v6.png\n")


# --- 7. Superimposed plot ---
cat("Rendering superimposed plot...\n")

df_superimposed <- df_aligned %>%
  group_by(Sample) %>%
  mutate(
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
    title    = "All samples superimposed -- aligned, sign-fixed, normalized",
    subtitle = "Each blue outline = one sample  |  Mass biased rightward",
    x = "X (normalized)", y = "Y (normalized)"
  ) +
  theme(legend.position = "none")

png("superimposed_signfix_v6.png", width = 2400, height = 2400, res = 200)
print(p_super)
dev.off()
cat("Saved -> superimposed_signfix_v6.png\n")


# --- 8. Export ---
write_csv(df_aligned, "phalloidin_aligned_signfix_v6.csv")

write_csv(mass_check %>% select(Sample, mean_x, pct_right, needs_flip),
          "sign_flip_log_v6.csv")

cat("\n=== Sign-fix outputs ===\n")
cat("  phalloidin_aligned_signfix_v6.csv  -- aligned fibers with sign corrected\n")
cat("  sign_flip_log_v6.csv               -- which samples were flipped\n")
cat("  aligned_signfix_overlay_v6.png     -- per-sample overlay\n")
cat("  superimposed_signfix_v6.png        -- all hulls superimposed\n")
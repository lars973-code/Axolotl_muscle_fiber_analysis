library(tidyverse)
library(sf)
library(FNN)       # fast k-nearest-neighbor search


# =============================================================================
# 1. LOAD DATA & FILTER TO OLD + PHALLOIDIN
#    Identical to v5 -- only Old tissue, Phalloidin channel.
# =============================================================================

data_full <- read_csv("COMPLETE_COMBINED_MUSCLE_FIBER_DATA.csv")
data_full <- data_full %>% filter(Age == "Old")
cat("After Age=='Old' filter:", nrow(data_full), "rows\n")

df_phalloidin <- data_full %>%
  filter(Channel_Code == "C1") %>%
  rename(fiber_index = `...1`)

cat("Phalloidin rows (Old only):", nrow(df_phalloidin), "\n")
cat("Samples:", n_distinct(df_phalloidin$Sample), "\n")


# =============================================================================
# 2. OUTLIER METRICS
#    Three complementary metrics computed per fiber, per sample:
#
#    (A) mean_dist       -- mean Euclidean distance to ALL other fibers.
#                           Catches globally isolated fibers.
#
#    (B) local_linearity -- ratio of 1st to 2nd eigenvalue from PCA on the
#                           k nearest neighbors.  High ratio = neighbors are
#                           arranged in a line (peninsulas / spurs).
#
#    (C) angular_coverage -- fraction of angular sectors (around each fiber)
#                            that contain at least one neighbor within a
#                            radius.  Low coverage = neighbors only on one
#                            side (boundary peaks).
#
#    All three are z-scored within each sample so thresholds are comparable
#    across samples of different sizes/densities.
# =============================================================================

K_NEIGHBORS      <- 15     # neighbors for linearity & angular metrics
N_SECTORS        <- 8      # angular bins (45° each) for coverage metric
RADIUS_QUANTILE  <- 0.10   # angular coverage radius = this quantile of
#   all pairwise nearest-neighbor distances

# --- Helper: compute all three metrics for one sample ---
compute_outlier_metrics <- function(df_s, key) {
  xy <- as.matrix(df_s[, c("X", "Y")])
  n  <- nrow(xy)
  
  if (n < K_NEIGHBORS + 1) {
    warning("Sample ", key$Sample, " has only ", n,
            " fibers -- skipping outlier metrics.")
    return(df_s %>% mutate(mean_dist = NA_real_,
                           local_linearity = NA_real_,
                           angular_coverage = NA_real_))
  }
  
  # ----------------------------------------------------------
  # (A) Mean distance to all other fibers
  # ----------------------------------------------------------
  # Use dist matrix for exact computation (feasible for <~10k fibers)
  d_mat     <- as.matrix(dist(xy))
  mean_dist <- rowMeans(d_mat)
  
  # ----------------------------------------------------------
  # (B) Local linearity via PCA on k nearest neighbors
  # ----------------------------------------------------------
  knn_res <- get.knn(xy, k = K_NEIGHBORS)  # excludes self
  knn_idx <- knn_res$nn.index
  
  local_linearity <- map_dbl(seq_len(n), function(i) {
    nb_xy <- xy[knn_idx[i, ], , drop = FALSE]
    # Include focal fiber in the neighborhood for PCA
    pts   <- rbind(xy[i, ], nb_xy)
    ev    <- prcomp(pts, center = TRUE, scale. = FALSE)$sdev^2
    # Ratio of 1st to 2nd eigenvalue (high = linear)
    if (length(ev) >= 2 && ev[2] > 0) ev[1] / ev[2] else Inf
  })
  
  # ----------------------------------------------------------
  # (C) Angular coverage
  #     For each fiber, look at its k nearest neighbors.
  #     Divide the circle into N_SECTORS equal slices.
  #     Score = fraction of sectors that contain >= 1 neighbor
  #             within a local radius.
  # ----------------------------------------------------------
  # Radius: use a generous local distance (quantile of knn distances)
  knn_dists   <- knn_res$nn.dist
  radius      <- quantile(as.vector(knn_dists[, 1]), RADIUS_QUANTILE)
  # Use a slightly larger radius: median of k-th neighbor distance
  
  radius <- median(knn_dists[, K_NEIGHBORS])
  
  sector_breaks <- seq(-pi, pi, length.out = N_SECTORS + 1)
  
  angular_coverage <- map_dbl(seq_len(n), function(i) {
    nb_xy  <- xy[knn_idx[i, ], , drop = FALSE]
    dx     <- nb_xy[, 1] - xy[i, 1]
    dy     <- nb_xy[, 2] - xy[i, 2]
    dists  <- sqrt(dx^2 + dy^2)
    angles <- atan2(dy, dx)
    
    # Only consider neighbors within the radius
    in_range <- dists <= radius
    if (sum(in_range) == 0) return(0)
    
    angles_in <- angles[in_range]
    # Count occupied sectors
    sector_ids <- cut(angles_in, breaks = sector_breaks,
                      include.lowest = TRUE, labels = FALSE)
    n_occupied <- n_distinct(na.omit(sector_ids))
    n_occupied / N_SECTORS
  })
  
  df_s %>%
    mutate(mean_dist        = mean_dist,
           local_linearity  = local_linearity,
           angular_coverage = angular_coverage)
}

cat("\nComputing outlier metrics per sample...\n")

df_phalloidin <- df_phalloidin %>%
  group_by(Sample) %>%
  group_modify(compute_outlier_metrics) %>%
  ungroup()

cat("Outlier metrics computed.\n")


# =============================================================================
# 3. Z-SCORE WITHIN SAMPLE & COMPOSITE OUTLIER FLAG
#    A fiber is flagged for removal if it meets ANY TWO of three criteria:
#      - z_mean_dist        > THRESH_DIST       (globally far away)
#      - z_local_linearity  > THRESH_LINEARITY   (on a peninsula)
#      - z_angular_coverage < THRESH_COVERAGE_Z  (neighbors on one side)
#
#    Additionally, angular_coverage raw < HARD_COVERAGE removes fibers
#    with very poor surround regardless of other metrics.
# =============================================================================

THRESH_DIST       <- 2.0   # z-score: mean distance
THRESH_LINEARITY  <- 2.0   # z-score: linearity
THRESH_COVERAGE_Z <- -1.5  # z-score: angular coverage (low = bad)
HARD_COVERAGE     <- 0.25  # raw: fewer than 25% of sectors occupied

df_phalloidin <- df_phalloidin %>%
  group_by(Sample) %>%
  mutate(
    z_mean_dist        = (mean_dist        - mean(mean_dist, na.rm = TRUE))        / sd(mean_dist, na.rm = TRUE),
    z_local_linearity  = (local_linearity  - mean(local_linearity, na.rm = TRUE))  / sd(local_linearity, na.rm = TRUE),
    z_angular_coverage = (angular_coverage - mean(angular_coverage, na.rm = TRUE)) / sd(angular_coverage, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    flag_dist     = z_mean_dist        > THRESH_DIST,
    flag_linear   = z_local_linearity  > THRESH_LINEARITY,
    flag_coverage = z_angular_coverage < THRESH_COVERAGE_Z,
    flag_hard_cov = angular_coverage   < HARD_COVERAGE,
    
    # Composite: 2-of-3 soft flags OR hard coverage cutoff
    n_flags       = flag_dist + flag_linear + flag_coverage,
    outlier       = (n_flags >= 2) | flag_hard_cov
  )

cat("\nOutlier summary:\n")
df_phalloidin %>%
  group_by(Sample) %>%
  summarise(total    = n(),
            outliers = sum(outlier, na.rm = TRUE),
            pct_out  = round(100 * outliers / total, 1),
            .groups  = "drop") %>%
  print(n = Inf)

cat("\nGlobal outlier rate:",
    sum(df_phalloidin$outlier, na.rm = TRUE), "of", nrow(df_phalloidin),
    sprintf("(%.1f%%)\n",
            100 * mean(df_phalloidin$outlier, na.rm = TRUE)))


# =============================================================================
# 4. CREATE CLEANED DATA (outliers removed)
# =============================================================================

df_clean <- df_phalloidin %>% filter(!outlier)
cat("\nCleaned fibers:", nrow(df_clean), "\n")


# =============================================================================
# 5. CONCAVE HULL PERIMETER SWEEP ON CLEANED DATA
# =============================================================================

RATIO_MIN <- 0.5
RATIO_MAX <- 1.0
N_STEPS   <- 50
ratio_vec <- seq(RATIO_MIN, RATIO_MAX, length.out = N_STEPS)

hull_perimeter <- function(xy_mat, ratio) {
  tryCatch({
    h <- st_concave_hull(st_multipoint(xy_mat), ratio = ratio)
    st_length(st_cast(h, "MULTILINESTRING"))
  }, error = function(e) NA_real_)
}

cat("\nSweeping concave hull ratios on cleaned data...\n")
sweep_results <- df_clean %>%
  group_by(Sample) %>%
  group_modify(function(df_s, key) {
    xy <- as.matrix(df_s[, c("X", "Y")])
    tibble(ratio = ratio_vec,
           perimeter = map_dbl(ratio_vec, ~hull_perimeter(xy, .x)))
  }) %>%
  ungroup()
write_csv(sweep_results, "hull_perimeter_sweep_v6.csv")
cat("Sweep complete.\n")


# =============================================================================
# 6. SECOND DERIVATIVE & OPTIMAL RATIO
# =============================================================================

optimal_ratios <- sweep_results %>%
  group_by(Sample) %>%
  arrange(ratio, .by_group = TRUE) %>%
  mutate(
    d1 = (lead(perimeter) - lag(perimeter)) / (lead(ratio) - lag(ratio)),
    d2 = (lead(d1)        - lag(d1))        / (lead(ratio) - lag(ratio))
  ) %>%
  filter(!is.na(d2)) %>%
  slice_min(d2, n = 1, with_ties = FALSE) %>%
  select(Sample,
         optimal_ratio    = ratio,
         perimeter_at_opt = perimeter,
         d2_at_opt        = d2) %>%
  ungroup()


# =============================================================================
# 7. MANUAL OVERRIDES
#    Carry forward the v5 overrides.  These may no longer be needed after
#    outlier removal -- inspect diagnostics and update as necessary.
# =============================================================================

MANUAL_OVERRIDES <- tribble(
  ~Sample,      ~optimal_ratio,  ~note,
  "a4_0001",    0.85,            "carried from v5 -- verify after outlier cleanup",
  "a4_0005",    0.85,            "carried from v5 -- verify after outlier cleanup",
  "a4_0007",    0.85,            "carried from v5 -- verify after outlier cleanup"
)

if (nrow(MANUAL_OVERRIDES) > 0) {
  # Only override samples that exist in the data
  MANUAL_OVERRIDES <- MANUAL_OVERRIDES %>%
    filter(Sample %in% optimal_ratios$Sample)
  
  if (nrow(MANUAL_OVERRIDES) > 0) {
    optimal_ratios <- optimal_ratios %>%
      rows_update(
        MANUAL_OVERRIDES %>% select(Sample, optimal_ratio),
        by = "Sample"
      ) %>%
      left_join(MANUAL_OVERRIDES %>% select(Sample, note), by = "Sample") %>%
      mutate(overridden = !is.na(note))
  } else {
    optimal_ratios <- optimal_ratios %>%
      mutate(note = NA_character_, overridden = FALSE)
  }
} else {
  optimal_ratios <- optimal_ratios %>%
    mutate(note = NA_character_, overridden = FALSE)
}

cat("\nOptimal ratios per sample:\n")
print(optimal_ratios %>% select(Sample, optimal_ratio, overridden, note), n = Inf)


# =============================================================================
# 8. BUILD FINAL HULLS ON CLEANED DATA
# =============================================================================

cat("\nBuilding hulls on cleaned data...\n")

hull_list <- map2(
  optimal_ratios$Sample,
  optimal_ratios$optimal_ratio,
  function(s, ratio) {
    xy <- df_clean %>% filter(Sample == s) %>%
      select(X, Y) %>% as.matrix()
    tryCatch(
      st_concave_hull(st_multipoint(xy), ratio = ratio),
      error = function(e) {
        warning("Hull failed for ", s, " at ratio=", ratio,
                " -- falling back to convex hull.")
        st_convex_hull(st_multipoint(xy))
      }
    )
  }
)
names(hull_list) <- optimal_ratios$Sample
cat("Hulls built.\n")


# =============================================================================
# 9. WITHIN-HULL FLAGGING (on cleaned fibers)
# =============================================================================

cat("Flagging fibers within hull...\n")

hulls_sf <- st_sf(
  Sample   = names(hull_list),
  geometry = st_sfc(hull_list),
  crs      = NA
)

within_flags <- map_dfr(names(hull_list), function(s) {
  df_s  <- df_clean %>% filter(Sample == s)
  pts_s <- st_as_sf(df_s %>% select(fiber_index, X, Y),
                    coords = c("X", "Y"), crs = NA)
  hull_s <- hulls_sf %>% filter(Sample == s)
  hits   <- st_intersects(pts_s, hull_s, sparse = TRUE)
  tibble(fiber_index = df_s$fiber_index,
         within_hull = lengths(hits) > 0)
})

df_clean <- df_clean %>%
  left_join(within_flags, by = "fiber_index")

cat("Flagging done.\n")
cat("\nFibers outside hull (post-cleaning):",
    sum(!df_clean$within_hull, na.rm = TRUE),
    sprintf("(%.1f%%)\n",
            100 * mean(!df_clean$within_hull, na.rm = TRUE)))


# =============================================================================
# 10. DIAGNOSTIC PLOTS
# =============================================================================

# --- 10a. Outlier metric distributions per sample ---
cat("\nRendering outlier metric diagnostic...\n")

metric_long <- df_phalloidin %>%
  select(Sample, fiber_index, outlier,
         z_mean_dist, z_local_linearity, z_angular_coverage) %>%
  pivot_longer(cols = starts_with("z_"),
               names_to = "metric", values_to = "z_value")

p_metrics <- ggplot(metric_long, aes(x = z_value, fill = outlier)) +
  geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
  scale_fill_manual(values = c("FALSE" = "#00c875", "TRUE" = "#ff4d4d"),
                    labels = c("FALSE" = "Kept", "TRUE" = "Outlier")) +
  facet_wrap(~metric, scales = "free", ncol = 3) +
  theme_minimal(base_size = 9) +
  labs(title    = "Outlier metric z-score distributions (all Old samples pooled)",
       subtitle = "Green = kept  |  Red = flagged outlier",
       x = "Z-score", y = "Count", fill = "Status") +
  theme(legend.position = "bottom")

png("outlier_metrics_hist_v6.png", width = 3600, height = 1800, res = 200)
print(p_metrics)
dev.off()
cat("Saved -> outlier_metrics_hist_v6.png\n")


# --- 10b. Per-sample scatter: outliers highlighted ---
cat("Rendering per-sample outlier scatter...\n")

p_outlier_scatter <- ggplot(df_phalloidin,
                            aes(x = X, y = Y, colour = outlier)) +
  geom_point(size = 0.15, alpha = 0.5) +
  scale_colour_manual(values = c("FALSE" = "#00c875", "TRUE" = "#ff4d4d"),
                      labels = c("FALSE" = "Kept", "TRUE" = "Outlier")) +
  facet_wrap(~Sample, scales = "free", ncol = 6) +
  theme_minimal(base_size = 7) +
  labs(title    = "Outlier detection -- Old samples (* fibers flagged for removal)",
       subtitle = "Green = kept  |  Red = outlier",
       x = "X (pixels)", y = "Y (pixels)", colour = "Status") +
  theme(strip.text      = element_text(size = 6),
        axis.text       = element_text(size = 4),
        panel.spacing   = unit(0.3, "lines"),
        legend.position = "bottom")

png("outlier_scatter_v6.png", width = 5400, height = 3600, res = 200)
print(p_outlier_scatter)
dev.off()
cat("Saved -> outlier_scatter_v6.png\n")


# --- 10c. Hull overlay on cleaned data ---
cat("Rendering hull overlay on cleaned data...\n")

hull_polygons <- imap_dfr(hull_list, function(geom, s) {
  tryCatch({
    coords        <- as.data.frame(st_coordinates(geom))
    coords$Sample <- s
    coords
  }, error = function(e) NULL)
})

df_thin <- df_clean %>%
  group_by(Sample) %>%
  slice_sample(n = 400) %>%
  ungroup()

p_hull <- ggplot() +
  geom_point(data = df_thin,
             aes(x = X, y = Y, colour = within_hull),
             size = 0.15, alpha = 0.35) +
  geom_polygon(data        = hull_polygons,
               aes(x = X, y = Y),
               fill        = NA, colour = "dodgerblue", linewidth = 0.45,
               inherit.aes = FALSE) +
  geom_text(data      = optimal_ratios,
            aes(label = paste0("r=", round(optimal_ratio, 3),
                               if_else(overridden, "*", "")),
                x     = -Inf, y = Inf),
            hjust = -0.1, vjust = 1.5,
            colour = "#00c875", size = 1.8, inherit.aes = FALSE) +
  scale_colour_manual(values = c("TRUE" = "#00c875", "FALSE" = "#ff4d4d")) +
  facet_wrap(~Sample, scales = "free") +
  theme_minimal(base_size = 7) +
  labs(
    title    = "v6 hull overlay -- Old samples, outliers removed (* = manual override)",
    subtitle = "Green = within hull  |  Red = outside hull  |  Blue = hull boundary",
    x = "X (pixels)", y = "Y (pixels)", colour = "Within Hull"
  ) +
  theme(strip.text      = element_text(size = 6),
        axis.text       = element_text(size = 4),
        panel.spacing   = unit(0.3, "lines"),
        legend.position = "bottom")

png("hull_overlay_v6.png", width = 4800, height = 3600, res = 200)
print(p_hull)
dev.off()
cat("Saved -> hull_overlay_v6.png\n")


# --- 10d. Sweep d2 curves ---
cat("Rendering sweep diagnostic...\n")

d2_plot_data <- sweep_results %>%
  group_by(Sample) %>%
  arrange(ratio, .by_group = TRUE) %>%
  mutate(
    d1 = (lead(perimeter) - lag(perimeter)) / (lead(ratio) - lag(ratio)),
    d2 = (lead(d1) - lag(d1))               / (lead(ratio) - lag(ratio))
  ) %>%
  filter(!is.na(d2)) %>%
  ungroup()

p_sweep <- ggplot(d2_plot_data) +
  geom_line(aes(x = ratio, y = scale(perimeter)),
            colour = "grey60", linewidth = 0.4) +
  geom_line(aes(x = ratio, y = scale(d2)),
            colour = "#f4a261", linewidth = 0.4, linetype = "dashed") +
  geom_vline(data = optimal_ratios,
             aes(xintercept = optimal_ratio,
                 colour     = overridden),
             linewidth = 0.55) +
  scale_colour_manual(values = c("FALSE" = "#00c875", "TRUE" = "#ff6b6b"),
                      labels = c("FALSE" = "Elbow (auto)", "TRUE" = "Manual override"),
                      name   = "Ratio source") +
  facet_wrap(~Sample, scales = "free_y", ncol = 6) +
  theme_minimal(base_size = 7) +
  labs(
    title    = "v6 perimeter sweep & d2 -- cleaned Old samples",
    subtitle = "Grey = scaled perimeter  |  Orange = scaled d2  |  Green = auto  |  Red = override",
    x = "Hull ratio", y = "Scaled value"
  ) +
  theme(strip.text      = element_text(size = 6),
        axis.text       = element_text(size = 4),
        panel.spacing   = unit(0.3, "lines"),
        legend.position = "bottom")

png("hull_sweep_v6.png", width = 5400, height = 3600, res = 200)
print(p_sweep)
dev.off()
cat("Saved -> hull_sweep_v6.png\n")


# =============================================================================
# 11. EXPORT
# =============================================================================

# Full data with outlier metrics and flags (all fibers, including outliers)
write_csv(df_phalloidin, "phalloidin_outlier_metrics_v6.csv")

# Cleaned data with hull flags (outliers removed)
write_csv(df_clean, "phalloidin_clean_with_hull_v6.csv")

# Optimal ratios
write_csv(optimal_ratios %>% select(Sample, optimal_ratio, overridden, note),
          "sample_optimal_ratios_v6.csv")

cat("\n=== v6 outputs ===\n")
cat("  phalloidin_outlier_metrics_v6.csv  -- all Old fibers + 3 outlier metrics + flags\n")
cat("  phalloidin_clean_with_hull_v6.csv  -- cleaned fibers (outliers removed) + hull flag\n")
cat("  sample_optimal_ratios_v6.csv       -- per-sample optimal ratios\n")
cat("  hull_perimeter_sweep_v6.csv        -- full sweep data on cleaned fibers\n")
cat("  outlier_metrics_hist_v6.png        -- z-score distributions for 3 metrics\n")
cat("  outlier_scatter_v6.png             -- per-sample scatter with outliers highlighted\n")
cat("  hull_overlay_v6.png                -- hull overlay on cleaned data\n")
cat("  hull_sweep_v6.png                  -- d2 sweep curves on cleaned data\n")
library(tidyverse)
library(sf)


# =============================================================================
# RESUME SCRIPT — picks up from saved hull_perimeter_sweep.csv
# No re-sweep needed. All prior fixes applied:
#   - fast st_intersects flagging (not st_join)
#   - thinned overlay plot (400 pts/sample max)
#   - boundary sample extension for a4_0005 and a6_0003
# =============================================================================


# =============================================================================
# 1. LOAD DATA
# =============================================================================

cat("Loading data...\n")
sweep_base    <- read_csv("hull_perimeter_sweep.csv")
data_full     <- read_csv("COMPLETE_COMBINED_MUSCLE_FIBER_DATA.csv")

df_phalloidin <- data_full %>%
  filter(Channel_Code == "C1") %>%
  rename(fiber_index = `...1`)

cat("Phalloidin fibers:", nrow(df_phalloidin), "\n")
cat("Sweep rows loaded:", nrow(sweep_base), "\n")


# =============================================================================
# 2. EXTEND BOUNDARY SAMPLES (a4_0005 and a6_0003 were near sweep edges)
# =============================================================================

RATIO_MIN        <- 0.6
RATIO_MAX        <- 1.0
BOUNDARY_SAMPLES <- c("a4_0005", "a6_0003")

hull_perimeter <- function(xy_mat, ratio) {
  tryCatch({
    h <- st_concave_hull(st_multipoint(xy_mat), ratio = ratio)
    st_length(st_cast(h, "MULTILINESTRING"))
  }, error = function(e) NA_real_)
}

extended_sweeps <- map_dfr(BOUNDARY_SAMPLES, function(s) {
  existing <- sweep_base %>% filter(Sample == s) %>% arrange(ratio)
  p  <- existing$perimeter
  r  <- existing$ratio
  d1 <- (lead(p) - lag(p)) / (lead(r) - lag(r))
  d2 <- (lead(d1) - lag(d1)) / (lead(r) - lag(r))
  opt_ratio <- existing$ratio[which.min(d2)]
  
  extend_low  <- opt_ratio < (RATIO_MIN + 0.1)
  extend_high <- opt_ratio > (RATIO_MAX - 0.1)
  
  new_ratios <- c(
    if (extend_low)  seq(max(0.01, RATIO_MIN - 0.15), RATIO_MIN, length.out = 15) else NULL,
    if (extend_high) seq(RATIO_MAX, min(0.85, RATIO_MAX + 0.15), length.out = 15) else NULL
  )
  new_ratios <- setdiff(round(new_ratios, 4), round(existing$ratio, 4))
  if (length(new_ratios) == 0) return(NULL)
  
  xy <- df_phalloidin %>% filter(Sample == s) %>% select(X, Y) %>% as.matrix()
  cat("Extending sweep for", s, ":", length(new_ratios), "extra steps\n")
  tibble(Sample = s, ratio = new_ratios,
         perimeter = map_dbl(new_ratios, ~hull_perimeter(xy, .x)))
})

sweep_full <- bind_rows(sweep_base, extended_sweeps) %>%
  arrange(Sample, ratio)

cat("Sweep extended. Total rows:", nrow(sweep_full), "\n")


# =============================================================================
# 3. SECOND DERIVATIVES & BOTH OPTIMAL RATIOS
# =============================================================================

sweep_d2 <- sweep_full %>%
  group_by(Sample) %>%
  arrange(ratio, .by_group = TRUE) %>%
  mutate(
    d1 = (lead(perimeter) - lag(perimeter)) / (lead(ratio) - lag(ratio)),
    d2 = (lead(d1) - lag(d1))               / (lead(ratio) - lag(ratio))
  ) %>%
  filter(!is.na(d2)) %>%
  ungroup()

optimal_elbow <- sweep_d2 %>%
  group_by(Sample) %>%
  slice_min(d2, n = 1, with_ties = FALSE) %>%
  transmute(Sample, ratio_elbow = ratio, perim_elbow = perimeter, d2_elbow = d2)

optimal_plateau <- sweep_d2 %>%
  group_by(Sample) %>%
  slice_min(abs(d2), n = 1, with_ties = FALSE) %>%
  transmute(Sample, ratio_plateau = ratio, perim_plateau = perimeter, d2_plateau = d2)

comparison_tbl <- left_join(optimal_elbow, optimal_plateau, by = "Sample") %>%
  mutate(
    ratio_diff     = round(ratio_plateau - ratio_elbow, 3),
    pct_perim_diff = round(100 * (perim_plateau - perim_elbow) / perim_elbow, 2)
  )

cat("\n=== Method comparison ===\n")
print(comparison_tbl %>% select(Sample, ratio_elbow, ratio_plateau, ratio_diff, pct_perim_diff), n = Inf)

cat("\nSummary of ratio differences (plateau - elbow):\n")
cat(sprintf("  Mean:   %+.3f\n", mean(comparison_tbl$ratio_diff)))
cat(sprintf("  Median: %+.3f\n", median(comparison_tbl$ratio_diff)))
cat(sprintf("  Range:  [%+.3f, %+.3f]\n", min(comparison_tbl$ratio_diff), max(comparison_tbl$ratio_diff)))

large_diff <- comparison_tbl %>% filter(abs(ratio_diff) > 0.15)
if (nrow(large_diff) > 0) {
  cat("\nSamples with large disagreement (|diff| > 0.15):\n")
  print(large_diff %>% select(Sample, ratio_elbow, ratio_plateau, ratio_diff))
}


# =============================================================================
# 4. BUILD HULLS FOR BOTH METHODS
# =============================================================================

cat("\nBuilding hulls...\n")

build_hull <- function(sample_name, ratio) {
  xy <- df_phalloidin %>% filter(Sample == sample_name) %>%
    select(X, Y) %>% as.matrix()
  tryCatch(
    st_concave_hull(st_multipoint(xy), ratio = ratio),
    error = function(e) {
      warning("Falling back to convex hull for ", sample_name)
      st_convex_hull(st_multipoint(xy))
    }
  )
}

hulls_elbow   <- map2(comparison_tbl$Sample, comparison_tbl$ratio_elbow,   build_hull)
hulls_plateau <- map2(comparison_tbl$Sample, comparison_tbl$ratio_plateau,  build_hull)
names(hulls_elbow)   <- comparison_tbl$Sample
names(hulls_plateau) <- comparison_tbl$Sample

cat("Hulls built.\n")


# =============================================================================
# 5. FAST WITHIN-HULL FLAGGING (st_intersects per sample — not global st_join)
# =============================================================================

flag_within_hull <- function(hull_named_list, label_col) {
  hulls_sf <- st_sf(
    Sample   = names(hull_named_list),
    geometry = st_sfc(hull_named_list),
    crs      = NA
  )
  map_dfr(names(hull_named_list), function(s) {
    df_s  <- df_phalloidin %>% filter(Sample == s)
    pts_s <- st_as_sf(df_s %>% select(fiber_index, X, Y),
                      coords = c("X", "Y"), crs = NA)
    hull_s <- hulls_sf %>% filter(Sample == s)
    hits   <- st_intersects(pts_s, hull_s, sparse = TRUE)
    tibble(fiber_index = df_s$fiber_index,
           !!label_col := lengths(hits) > 0)
  })
}

cat("Flagging elbow hulls...\n")
flags_elbow   <- flag_within_hull(hulls_elbow,   "within_elbow")
cat("Flagging plateau hulls...\n")
flags_plateau <- flag_within_hull(hulls_plateau, "within_plateau")
cat("Flagging done.\n")

df_compare <- df_phalloidin %>%
  left_join(flags_elbow,   by = "fiber_index") %>%
  left_join(flags_plateau, by = "fiber_index") %>%
  mutate(
    agreement = case_when(
      within_elbow &  within_plateau ~ "both_inside",
      !within_elbow & !within_plateau ~ "both_outside",
      within_elbow & !within_plateau ~ "elbow_only",
      !within_elbow &  within_plateau ~ "plateau_only"
    )
  )

cat("\n=== Agreement between methods ===\n")
print(table(df_compare$agreement))

cat("\nPer-sample disagreement:\n")
df_compare %>%
  group_by(Sample) %>%
  summarise(total        = n(),
            disagreement = sum(agreement %in% c("elbow_only", "plateau_only")),
            pct_disagree = round(100 * disagreement / total, 1)) %>%
  arrange(desc(pct_disagree)) %>%
  print(n = Inf)


# =============================================================================
# 6. DIAGNOSTIC PLOTS
# =============================================================================

# --- 6a. d² curves — both picks marked per sample ---
cat("\nRendering d² curve plot...\n")

p_d2 <- ggplot(sweep_d2) +
  geom_line(aes(x = ratio, y = scale(perimeter)),
            colour = "white", linewidth = 0.35) +
  geom_line(aes(x = ratio, y = scale(d2)),
            colour = "#f4a261", linewidth = 0.35, linetype = "dashed") +
  geom_vline(data = comparison_tbl, aes(xintercept = ratio_elbow),
             colour = "#ff6b6b", linewidth = 0.55) +
  geom_vline(data = comparison_tbl, aes(xintercept = ratio_plateau),
             colour = "#00c875", linewidth = 0.55, linetype = "dotted") +
  facet_wrap(~Sample, scales = "free_y", ncol = 6) +
  theme_dark(base_size = 7) +
  labs(title    = "Perimeter sweep & d² — elbow (red) vs plateau (green) per sample",
       subtitle = "White = scaled perimeter  |  Orange = scaled d²  |  Red = min(d²)  |  Green = min(|d²|)",
       x = "Hull ratio", y = "Scaled value") +
  theme(strip.text    = element_text(size = 6, colour = "white"),
        axis.text     = element_text(size = 4),
        panel.spacing = unit(0.25, "lines"))

png("comparison_d2_curves.png", width = 5400, height = 3600, res = 200)
print(p_d2)
dev.off()
cat("Saved → comparison_d2_curves.png\n")


# --- 6b. Hull overlay — thinned to 400 pts/sample for rendering speed ---
cat("Rendering hull overlay plot (thinned)...\n")

get_hull_coords <- function(hull_named_list, method_label) {
  imap_dfr(hull_named_list, function(geom, s) {
    tryCatch({
      coords        <- as.data.frame(st_coordinates(geom))
      coords$Sample <- s
      coords$method <- method_label
      coords
    }, error = function(e) NULL)
  })
}

hull_coords <- bind_rows(
  get_hull_coords(hulls_elbow,   "Elbow — min(d2)"),
  get_hull_coords(hulls_plateau, "Plateau — min(abs(d2))")
)

# Thin points for rendering — 400 per sample max
df_compare_thin <- df_compare %>%
  group_by(Sample) %>%
  slice_sample(n = 400) %>%
  ungroup()

agree_colours <- c(
  "both_inside"  = "#aaaaaa",
  "both_outside" = "#333333",
  "elbow_only"   = "#ff6b6b",
  "plateau_only" = "#00c875"
)

p_overlay <- ggplot() +
  geom_point(data  = df_compare_thin,
             aes(x = X, y = Y, colour = agreement),
             size = 0.15, alpha = 0.4) +
  geom_polygon(data = hull_coords,
               aes(x = X, y = Y, colour = method, linetype = method),
               fill = NA, linewidth = 0.5, inherit.aes = FALSE) +
  scale_colour_manual(values = c(
    agree_colours,
    "Elbow — min(d2)"       = "#ff6b6b",
    "Plateau — min(abs(d2))" = "#00c875"
  )) +
  scale_linetype_manual(values = c(
    "Elbow — min(d2)"        = "solid",
    "Plateau — min(abs(d2))" = "dotted"
  )) +
  facet_wrap(~Sample, scales = "free", ncol = 6) +
  theme_dark(base_size = 7) +
  labs(title    = "Hull comparison — Elbow (red solid) vs Plateau (green dotted)",
       subtitle = "Grey=both inside  |  Dark=both outside  |  Red=elbow only  |  Green=plateau only",
       x = "X (pixels)", y = "Y (pixels)") +
  theme(strip.text      = element_text(size = 6, colour = "white"),
        axis.text       = element_text(size = 4),
        panel.spacing   = unit(0.25, "lines"),
        legend.position = "bottom",
        legend.text     = element_text(size = 6))

png("comparison_hull_overlay.png", width = 5400, height = 4200, res = 200)
print(p_overlay)
dev.off()
cat("Saved → comparison_hull_overlay.png\n")


# --- 6c. Ratio scatter ---
cat("Rendering ratio scatter...\n")

p_scatter <- ggplot(comparison_tbl,
                    aes(x = ratio_elbow, y = ratio_plateau, label = Sample)) +
  geom_abline(slope = 1, intercept = 0,
              colour = "white", linetype = "dashed", linewidth = 0.4) +
  geom_point(aes(colour = abs(ratio_diff)), size = 3) +
  geom_text(aes(label = Sample), size = 2, colour = "white",
            hjust = -0.15, vjust = 0.5) +
  scale_colour_gradient(low = "#00c875", high = "#ff6b6b", name = "|ratio diff|") +
  theme_dark(base_size = 9) +
  labs(title    = "Elbow vs Plateau optimal ratio — per sample",
       subtitle = "Points on dashed line = perfect agreement; red = large disagreement",
       x = "Elbow ratio — min(d2)",
       y = "Plateau ratio — min(abs(d2))")

png("comparison_ratio_scatter.png", width = 2000, height = 1800, res = 200)
print(p_scatter)
dev.off()
cat("Saved → comparison_ratio_scatter.png\n")


# =============================================================================
# 7. EXPORT
# =============================================================================

write_csv(comparison_tbl, "sample_ratio_comparison.csv")
write_csv(
  df_compare %>% select(fiber_index, Sample, X, Y, Area, Perim., Mean,
                        within_elbow, within_plateau, agreement),
  "phalloidin_hull_comparison.csv"
)

cat("\nAll outputs written:\n")
cat("  sample_ratio_comparison.csv     — per-sample elbow vs plateau ratios\n")
cat("  phalloidin_hull_comparison.csv  — fibers with both hull flags + agreement\n")
cat("  comparison_d2_curves.png        — d² curves with both picks\n")
cat("  comparison_hull_overlay.png     — hull overlay (thinned scatter)\n")
cat("  comparison_ratio_scatter.png    — scatter of elbow vs plateau ratios\n")
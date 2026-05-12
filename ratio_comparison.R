install.packages("patchwork")
library(tidyverse)
library(sf)
library(patchwork)   # for side-by-side plot panels

# =============================================================================
# PREREQUISITE: This script reads outputs already produced by
# alignment_boundary_detect_v3.R — specifically:
#   - hull_perimeter_sweep.csv   (full ratio x perimeter table)
#   - COMPLETE_COMBINED_MUSCLE_FIBER_DATA.csv  (fiber positions)
# Run v3 through Section 2 (sweep) first if those files don't exist.
# =============================================================================

BOUNDARY_SAMPLES <- c("a4_0005", "a6_0003")  # optima near sweep edges in v3
RATIO_MIN        <- 0.05
RATIO_MAX        <- 0.60
N_STEPS          <- 50


# =============================================================================
# 1. LOAD SWEEP DATA & EXTEND BOUNDARY SAMPLES IF NEEDED
# =============================================================================

sweep_base <- read_csv("hull_perimeter_sweep.csv")

data_full    <- read_csv("COMPLETE_COMBINED_MUSCLE_FIBER_DATA.csv")
df_phalloidin <- data_full %>%
  filter(Channel_Code == "C1") %>%
  rename(fiber_index = `...1`)

hull_perimeter <- function(xy_mat, ratio) {
  tryCatch({
    h <- st_concave_hull(st_multipoint(xy_mat), ratio = ratio)
    st_length(st_cast(h, "MULTILINESTRING"))
  }, error = function(e) NA_real_)
}

# Extend sweep for boundary samples: push range out by 20% on the flagged side
extended_sweeps <- map_dfr(BOUNDARY_SAMPLES, function(s) {
  existing   <- sweep_base %>% filter(Sample == s) %>% arrange(ratio)
  opt_ratio  <- existing %>% slice_min(d2 <- {
    p  <- existing$perimeter
    r  <- existing$ratio
    d1 <- (lead(p) - lag(p)) / (lead(r) - lag(r))
    (lead(d1) - lag(d1)) / (lead(r) - lag(r))
  }, n = 1) %>% pull(ratio)
  
  # Determine which boundary is at risk and extend that side
  extend_low  <- opt_ratio < (RATIO_MIN + 0.05)
  extend_high <- opt_ratio > (RATIO_MAX - 0.05)
  
  new_ratios <- c(
    if (extend_low)  seq(max(0.01, RATIO_MIN - 0.15), RATIO_MIN, length.out = 15) else NULL,
    if (extend_high) seq(RATIO_MAX, min(0.85, RATIO_MAX + 0.15), length.out = 15) else NULL
  )
  new_ratios <- setdiff(round(new_ratios, 4), round(existing$ratio, 4))
  
  if (length(new_ratios) == 0) return(NULL)
  
  xy <- df_phalloidin %>% filter(Sample == s) %>% select(X, Y) %>% as.matrix()
  cat("Extending sweep for", s, ":", length(new_ratios), "extra steps\n")
  
  tibble(
    Sample    = s,
    ratio     = new_ratios,
    perimeter = map_dbl(new_ratios, ~hull_perimeter(xy, .x))
  )
})

sweep_full <- bind_rows(sweep_base, extended_sweeps) %>%
  arrange(Sample, ratio)


# =============================================================================
# 2. COMPUTE SECOND DERIVATIVES & DERIVE BOTH OPTIMAL RATIOS PER SAMPLE
# =============================================================================

sweep_d2 <- sweep_full %>%
  group_by(Sample) %>%
  arrange(ratio, .by_group = TRUE) %>%
  mutate(
    d1 = (lead(perimeter) - lag(perimeter)) / (lead(ratio) - lag(ratio)),
    d2 = (lead(d1)        - lag(d1))        / (lead(ratio) - lag(ratio))
  ) %>%
  filter(!is.na(d2)) %>%
  ungroup()

optimal_elbow <- sweep_d2 %>%
  group_by(Sample) %>%
  slice_min(d2, n = 1, with_ties = FALSE) %>%        # most negative d2
  transmute(Sample,
            ratio_elbow   = ratio,
            perim_elbow   = perimeter,
            d2_elbow      = d2)

optimal_plateau <- sweep_d2 %>%
  group_by(Sample) %>%
  slice_min(abs(d2), n = 1, with_ties = FALSE) %>%   # d2 closest to zero
  transmute(Sample,
            ratio_plateau = ratio,
            perim_plateau = perimeter,
            d2_plateau    = d2)

comparison_tbl <- left_join(optimal_elbow, optimal_plateau, by = "Sample") %>%
  mutate(
    ratio_diff     = round(ratio_plateau - ratio_elbow, 3),
    pct_perim_diff = round(100 * (perim_plateau - perim_elbow) / perim_elbow, 2)
  )

cat("\n=== Method comparison ===\n")
print(comparison_tbl %>%
        select(Sample, ratio_elbow, ratio_plateau, ratio_diff, pct_perim_diff),
      n = Inf)

cat("\nSummary of ratio differences (plateau - elbow):\n")
cat(sprintf("  Mean:   %+.3f\n", mean(comparison_tbl$ratio_diff)))
cat(sprintf("  Median: %+.3f\n", median(comparison_tbl$ratio_diff)))
cat(sprintf("  Range:  [%+.3f, %+.3f]\n",
            min(comparison_tbl$ratio_diff), max(comparison_tbl$ratio_diff)))

large_diff <- comparison_tbl %>% filter(abs(ratio_diff) > 0.15)
if (nrow(large_diff) > 0) {
  cat("\nSamples with large disagreement (|diff| > 0.15):\n")
  print(large_diff %>% select(Sample, ratio_elbow, ratio_plateau, ratio_diff))
}


# =============================================================================
# 3. BUILD HULLS FOR BOTH METHODS
# =============================================================================

build_hull <- function(sample_name, ratio) {
  xy <- df_phalloidin %>%
    filter(Sample == sample_name) %>%
    select(X, Y) %>% as.matrix()
  tryCatch(
    st_concave_hull(st_multipoint(xy), ratio = ratio),
    error = function(e) st_convex_hull(st_multipoint(xy))
  )
}

hulls_elbow   <- map2(comparison_tbl$Sample, comparison_tbl$ratio_elbow,   build_hull)
hulls_plateau <- map2(comparison_tbl$Sample, comparison_tbl$ratio_plateau,  build_hull)
names(hulls_elbow)   <- comparison_tbl$Sample
names(hulls_plateau) <- comparison_tbl$Sample


# =============================================================================
# 4. WITHIN-HULL FLAGS FOR BOTH METHODS
# =============================================================================

flag_within_hull <- function(hull_named_list, label_col) {
  
  # Build hull sf object once
  hulls_sf <- st_sf(
    Sample   = names(hull_named_list),
    geometry = st_sfc(hull_named_list),
    crs      = NA
  )
  
  # Process sample-by-sample — avoids building a 62k-point global sf object
  # and instead does 37 small intersects (~500-4000 pts each), which is 
  # dramatically faster due to smaller spatial index per call
  results <- map_dfr(names(hull_named_list), function(s) {
    
    df_s <- df_phalloidin %>% filter(Sample == s)
    
    pts_s <- st_as_sf(
      df_s %>% select(fiber_index, X, Y),
      coords = c("X", "Y"),
      crs    = NA
    )
    
    hull_s <- hulls_sf %>% filter(Sample == s)
    
    # st_intersects returns a sparse list — much faster than st_join
    hits <- st_intersects(pts_s, hull_s, sparse = TRUE)
    
    tibble(
      fiber_index    = df_s$fiber_index,
      !!label_col   := lengths(hits) > 0
    )
  })
  
  results
}


flags_elbow   <- flag_within_hull(hulls_elbow,   "within_elbow")
flags_plateau <- flag_within_hull(hulls_plateau, "within_plateau")

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

cat("\nPer-sample disagreement (fibers classified differently):\n")
df_compare %>%
  group_by(Sample) %>%
  summarise(
    total        = n(),
    disagreement = sum(agreement %in% c("elbow_only", "plateau_only")),
    pct_disagree = round(100 * disagreement / total, 1)
  ) %>%
  arrange(desc(pct_disagree)) %>%
  print(n = Inf)


# =============================================================================
# 5. DIAGNOSTIC PLOTS
# =============================================================================

# --- 5a. Per-sample d² curves with both method picks ---

p_d2 <- ggplot(sweep_d2) +
  geom_line(aes(x = ratio, y = scale(perimeter)),
            colour = "white", linewidth = 0.35) +
  geom_line(aes(x = ratio, y = scale(d2)),
            colour = "#f4a261", linewidth = 0.35, linetype = "dashed") +
  geom_vline(data = comparison_tbl,
             aes(xintercept = ratio_elbow),
             colour = "#ff6b6b", linewidth = 0.55, linetype = "solid") +
  geom_vline(data = comparison_tbl,
             aes(xintercept = ratio_plateau),
             colour = "#00c875", linewidth = 0.55, linetype = "dotted") +
  facet_wrap(~Sample, scales = "free_y", ncol = 6) +
  theme_dark(base_size = 7) +
  labs(
    title    = "Perimeter sweep & d² — elbow (red) vs plateau (green) per sample",
    subtitle = "White = scaled perimeter  |  Orange dashed = scaled d²  |  Red solid = min(d²)  |  Green dotted = min(|d²|)",
    x = "Hull ratio", y = "Scaled value"
  ) +
  theme(strip.text = element_text(size = 6, colour = "white"),
        axis.text  = element_text(size = 4),
        panel.spacing = unit(0.25, "lines"))

png("comparison_d2_curves.png", width = 5400, height = 3600, res = 200)
print(p_d2)
dev.off()
cat("\nPlot saved → comparison_d2_curves.png\n")


# --- 5b. Hull overlay: elbow vs plateau, side by side per sample ---

get_hull_coords <- function(hull_named_list, method_label) {
  imap_dfr(hull_named_list, function(geom, s) {
    tryCatch({
      coords         <- as.data.frame(st_coordinates(geom))
      coords$Sample  <- s
      coords$method  <- method_label
      coords
    }, error = function(e) NULL)
  })
}

hull_coords <- bind_rows(
  get_hull_coords(hulls_elbow,   "Elbow — min(d²)"),
  get_hull_coords(hulls_plateau, "Plateau — min(|d²|)")
)

# Colour fibers by agreement category
agree_colours <- c(
  "both_inside"  = "#aaaaaa",
  "both_outside" = "#444444",
  "elbow_only"   = "#ff6b6b",
  "plateau_only" = "#00c875"
)

p_overlay <- ggplot() +
  geom_point(
    data  = df_compare,
    aes(x = X, y = Y, colour = agreement),
    size  = 0.07, alpha = 0.3
  ) +
  geom_polygon(
    data = hull_coords,
    aes(x = X, y = Y, colour = method, linetype = method),
    fill        = NA, linewidth = 0.45,
    inherit.aes = FALSE
  ) +
  scale_colour_manual(values = c(
    agree_colours,
    "Elbow — min(d²)"     = "#ff6b6b",
    "Plateau — min(|d²|)" = "#00c875"
  )) +
  scale_linetype_manual(values = c(
    "Elbow — min(d²)"     = "solid",
    "Plateau — min(|d²|)" = "dotted"
  )) +
  facet_wrap(~Sample, scales = "free", ncol = 6) +
  theme_dark(base_size = 7) +
  labs(
    title    = "Hull comparison — Elbow (red solid) vs Plateau (green dotted)",
    subtitle = "Fibers: grey=both inside  |  dark=both outside  |  red=elbow only  |  green=plateau only",
    x = "X (pixels)", y = "Y (pixels)"
  ) +
  theme(strip.text      = element_text(size = 6, colour = "white"),
        axis.text       = element_text(size = 4),
        panel.spacing   = unit(0.25, "lines"),
        legend.position = "bottom",
        legend.text     = element_text(size = 6))

png("comparison_hull_overlay.png", width = 5400, height = 4200, res = 200)
print(p_overlay)
dev.off()
cat("Plot saved → comparison_hull_overlay.png\n")


# --- 5c. Ratio scatter: how much do the two methods disagree per sample? ---

p_scatter <- ggplot(comparison_tbl,
                    aes(x = ratio_elbow, y = ratio_plateau, label = Sample)) +
  geom_abline(slope = 1, intercept = 0,
              colour = "white", linetype = "dashed", linewidth = 0.4) +
  geom_point(aes(colour = abs(ratio_diff)), size = 3) +
  ggrepel::geom_text_repel(size = 2.5, colour = "white", max.overlaps = 20) +
  scale_colour_gradient(low = "#00c875", high = "#ff6b6b",
                        name = "|ratio diff|") +
  theme_dark(base_size = 9) +
  labs(
    title    = "Elbow vs Plateau optimal ratio — per sample",
    subtitle = "Points on dashed line = perfect agreement; red = large disagreement",
    x = "Elbow ratio — min(d²)",
    y = "Plateau ratio — min(|d²|)"
  )

png("comparison_ratio_scatter.png", width = 2000, height = 1800, res = 200)
print(p_scatter)
dev.off()
cat("Plot saved → comparison_ratio_scatter.png\n")


# =============================================================================
# 6. EXPORT
# =============================================================================

write_csv(comparison_tbl, "sample_ratio_comparison.csv")
write_csv(df_compare %>%
            select(fiber_index, Sample, X, Y, Area, Perim., Mean,
                   within_elbow, within_plateau, agreement),
          "phalloidin_hull_comparison.csv")

cat("\nOutputs:\n")
cat("  sample_ratio_comparison.csv     — per-sample elbow vs plateau ratios\n")
cat("  phalloidin_hull_comparison.csv  — fibers with both hull flags + agreement\n")
cat("  comparison_d2_curves.png        — d² sweep curves, both picks marked\n")
cat("  comparison_hull_overlay.png     — hull shapes overlaid on fibers\n")
cat("  comparison_ratio_scatter.png    — scatter of elbow vs plateau ratios\n")
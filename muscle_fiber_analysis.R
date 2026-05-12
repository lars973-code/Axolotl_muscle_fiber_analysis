# =============================================================================
# MUSCLE FIBER ANALYSIS PIPELINE
# =============================================================================
# Analyses: mCherry positivity, fiber size profiles, radial density,
#           spatial asymmetry, DAPI nuclei, area vs intensity,
#           and per-sample summaries.
#
# Usage:
#   1. Edit CONFIG below for each new dataset
#   2. source("muscle_fiber_analysis.R")
#   3. All plots and tables are saved to output_dir
#
# Required input columns:
#   Area, Mean, X, Y, Perim., IntDen, Age, Animal, Image,
#   Channel, Sample, Positive, X_scaled, Y_scaled
#
# Channel values (case-insensitive):
#   Magenta -> Phalloidin  (all muscle fibers, geometry)
#   Red     -> mCherry     (labeled cells, intensity + Positive flag)
#   Blue    -> DAPI        (nuclei)
# =============================================================================

library(tidyverse)
library(ggplot2)
library(patchwork)

# =============================================================================
# CONFIG  <-- edit this section for each new dataset
# =============================================================================

CONFIG <- list(

  # Path to input CSV
  input_file    = "LFNG_scaled_fibers.csv",

  # Tissue type: "limb" or "tail"
  tissue_type   = "limb",

  # Radial zone method: "equal_width" | "quantile"
  #   equal_width -> 6 equal rings from 0 to r_max (good default)
  #   quantile    -> 6 zones each containing equal numbers of fibers
  zone_method   = "equal_width",
  n_zones       = 6,

  # Filter to specific age groups: NULL keeps all, e.g. c("Old") or c("Young")
  age_filter    = NULL,

  # Filter to specific samples: NULL keeps all
  sample_filter = NULL,

  # Min fiber area to keep (removes segmentation noise)
  min_area      = 5,

  # Output directory (created if it doesn't exist)
  output_dir    = "fiber_analysis_output",

  # Plot dimensions (inches)
  plot_width    = 10,
  plot_height   = 7,
  plot_dpi      = 200
)

# =============================================================================
# 0. SETUP
# =============================================================================

dir.create(CONFIG$output_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== Muscle Fiber Analysis Pipeline ===\n")
cat("Input :", CONFIG$input_file, "\n")
cat("Tissue:", CONFIG$tissue_type, "\n")
cat("Zones :", CONFIG$zone_method, "(n =", CONFIG$n_zones, ")\n\n")

# --- Colour palettes ---
RING_COLORS <- c("#185FA5","#378ADD","#85B7EB","#B5D4F4","#9FE1CB",
                 "#1D9E75","#5DCAA5","#0F6E56")

QUAD_COLORS <- c(Q1_NE = "#378ADD", Q2_NW = "#1D9E75",
                 Q3_SW = "#D85A30", Q4_SE = "#D4537E")

POS_COLORS  <- c(Positive = "#3C3489", Negative = "#888780")

# Age colours are built from the actual groups present in the data (section 1)
AGE_PALETTE <- c("#3C3489","#1D9E75","#D85A30","#D4537E","#888780")

# --- Theme & save helper ---
theme_fiber <- function() {
  theme_minimal(base_size = 11) +
    theme(
      plot.title    = element_text(size = 13, face = "bold", margin = margin(b = 6)),
      plot.subtitle = element_text(size = 10, colour = "grey45", margin = margin(b = 10)),
      axis.title    = element_text(size = 10, colour = "grey30"),
      axis.text     = element_text(size = 9),
      strip.text    = element_text(size = 9, face = "bold"),
      legend.title  = element_text(size = 9),
      legend.text   = element_text(size = 9),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey92"),
      plot.caption  = element_text(size = 8, colour = "grey55", margin = margin(t = 8))
    )
}

save_plot <- function(p, filename, w = CONFIG$plot_width, h = CONFIG$plot_height) {
  path <- file.path(CONFIG$output_dir, filename)
  ggsave(path, p, width = w, height = h, dpi = CONFIG$plot_dpi)
  cat("  Saved:", path, "\n")
}

# =============================================================================
# 1. LOAD & VALIDATE
# =============================================================================

cat("Loading data...\n")
raw <- read_csv(CONFIG$input_file, show_col_types = FALSE)
cat("  Raw rows:", nrow(raw), "\n")

required_cols <- c("Area", "Mean", "IntDen", "X_scaled", "Y_scaled",
                   "Channel", "Age", "Animal", "Sample", "Positive")
missing_cols <- setdiff(required_cols, names(raw))
if (length(missing_cols) > 0) {
  stop("Input CSV is missing required columns: ", paste(missing_cols, collapse = ", "))
}

# Apply optional filters
if (!is.null(CONFIG$age_filter))    raw <- raw %>% filter(Age %in% CONFIG$age_filter)
if (!is.null(CONFIG$sample_filter)) raw <- raw %>% filter(Sample %in% CONFIG$sample_filter)
raw <- raw %>% filter(Area >= CONFIG$min_area)

cat("  After filters:", nrow(raw), "rows |",
    n_distinct(raw$Sample), "samples |",
    n_distinct(raw$Age), "age group(s)\n\n")

# Build age colour map from whatever groups are present
ages       <- sort(unique(raw$Age))
AGE_COLORS <- setNames(AGE_PALETTE[seq_along(ages)], ages)

# =============================================================================
# 2. SPLIT CHANNELS
# =============================================================================

phal    <- raw %>% filter(str_to_lower(Channel) == "magenta")
mcherry <- raw %>% filter(str_to_lower(Channel) == "red")
dapi    <- raw %>% filter(str_to_lower(Channel) == "blue")

has_phal    <- nrow(phal)    > 0
has_mcherry <- nrow(mcherry) > 0
has_dapi    <- nrow(dapi)    > 0

cat("Channel counts:\n")
cat("  Phalloidin (Magenta):", nrow(phal),    "\n")
cat("  mCherry    (Red)    :", nrow(mcherry), "\n")
cat("  DAPI       (Blue)   :", nrow(dapi),    "\n\n")

if (!has_phal)    warning("No Phalloidin (Magenta) data — geometry analyses will be skipped.")
if (!has_mcherry) warning("No mCherry (Red) data — positivity analyses will be skipped.")
if (!has_dapi)    warning("No DAPI (Blue) data — nuclei analyses will be skipped.")

# Ensure Positive is logical (handles TRUE/FALSE strings, 0/1 integers, NA)
if (has_mcherry) {
  mcherry <- mcherry %>% mutate(Positive = as.logical(Positive))
}

# =============================================================================
# 3. SPATIAL FEATURES
# =============================================================================

cat("Computing spatial features...\n")

# --- Quadrant assignment ---
assign_quadrant <- function(df) {
  df %>% mutate(
    quadrant = case_when(
      X_scaled >= 0 & Y_scaled >= 0 ~ "Q1_NE",
      X_scaled <  0 & Y_scaled >= 0 ~ "Q2_NW",
      X_scaled <  0 & Y_scaled <  0 ~ "Q3_SW",
      TRUE                           ~ "Q4_SE"
    )
  )
}

# --- Radial zone assignment ---
# Zone breaks are computed ONCE from the phalloidin distribution so all
# channels share identical boundaries. Falls back to the full dataset if
# phalloidin is absent.

zone_labels <- switch(as.character(CONFIG$n_zones),
  "5" = c("Core", "Inner", "Mid", "Outer", "Edge"),
  "6" = c("Core", "Inner", "Mid", "Outer", "Edge", "Far"),
  paste0("Z", seq_len(CONFIG$n_zones))
)

compute_r_breaks <- function(df, method, n_zones) {
  r <- sqrt(df$X_scaled^2 + df$Y_scaled^2)
  if (method == "equal_width") {
    seq(0, max(r, na.rm = TRUE) * 1.001, length.out = n_zones + 1)
  } else {  # quantile
    breaks <- quantile(r, probs = seq(0, 1, length.out = n_zones + 1), na.rm = TRUE)
    breaks[length(breaks)] <- breaks[length(breaks)] * 1.001
    as.numeric(breaks)
  }
}

assign_radial_zones <- function(df, r_breaks) {
  df %>% mutate(
    r          = sqrt(X_scaled^2 + Y_scaled^2),
    radial_bin = cut(r, breaks = r_breaks, labels = zone_labels,
                     include.lowest = TRUE, right = TRUE)
  )
}

ref_data <- if (has_phal) phal else raw
r_breaks <- compute_r_breaks(ref_data, CONFIG$zone_method, CONFIG$n_zones)

if (has_phal)    phal    <- phal    %>% assign_quadrant() %>% assign_radial_zones(r_breaks)
if (has_mcherry) mcherry <- mcherry %>% assign_quadrant() %>% assign_radial_zones(r_breaks)
if (has_dapi)    dapi    <- dapi    %>% assign_quadrant() %>% assign_radial_zones(r_breaks)

# Enforce consistent factor ordering so plots always go Core -> Far
ring_pal <- setNames(RING_COLORS[seq_along(zone_labels)], zone_labels)

if (has_phal)    phal$radial_bin    <- factor(phal$radial_bin,    levels = zone_labels)
if (has_mcherry) mcherry$radial_bin <- factor(mcherry$radial_bin, levels = zone_labels)
if (has_dapi)    dapi$radial_bin    <- factor(dapi$radial_bin,    levels = zone_labels)

cat("  Zone levels:", paste(zone_labels, collapse = " > "), "\n")
cat("  Age groups :", paste(ages, collapse = ", "), "\n\n")

# =============================================================================
# 4. NORMALISE PER SAMPLE
# =============================================================================

cat("Normalising measurements per sample...\n")

minmax <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (diff(rng) == 0) return(rep(0, length(x)))
  (x - rng[1]) / diff(rng)
}

if (has_phal) {
  phal <- phal %>%
    group_by(Sample) %>%
    mutate(Area_norm  = minmax(Area),
           Perim_norm = minmax(`Perim.`)) %>%
    ungroup()
}

if (has_mcherry) {
  mcherry <- mcherry %>%
    group_by(Sample) %>%
    mutate(Area_norm   = minmax(Area),
           Mean_norm   = minmax(Mean),
           IntDen_norm = minmax(IntDen)) %>%
    ungroup() %>%
    mutate(Pos_label = if_else(Positive, "Positive", "Negative"))
}

# =============================================================================
# 5. SUMMARY TABLES
# =============================================================================

cat("Computing summary tables...\n")

# 5a. Overall channel counts
summary_counts <- bind_rows(
  if (has_phal)    phal    %>% summarise(channel = "Phalloidin", n = n(), samples = n_distinct(Sample), ages = paste(sort(unique(Age)), collapse = "/")),
  if (has_mcherry) mcherry %>% summarise(channel = "mCherry",   n = n(), samples = n_distinct(Sample), ages = paste(sort(unique(Age)), collapse = "/")),
  if (has_dapi)    dapi    %>% summarise(channel = "DAPI",      n = n(), samples = n_distinct(Sample), ages = paste(sort(unique(Age)), collapse = "/"))
)
write_csv(summary_counts, file.path(CONFIG$output_dir, "00_summary_counts.csv"))

# 5b. Phalloidin area by age × radial zone
if (has_phal) {
  phal_area_summary <- phal %>%
    group_by(Age, radial_bin) %>%
    summarise(n = n(), mean_area = mean(Area), median_area = median(Area),
              sd_area = sd(Area), se_area = sd(Area) / sqrt(n()), .groups = "drop")
  write_csv(phal_area_summary, file.path(CONFIG$output_dir, "01_phal_area_by_age_zone.csv"))
}

# 5c–5g. mCherry summaries
if (has_mcherry) {

  pos_rate_summary <- mcherry %>%
    group_by(Age, radial_bin) %>%
    summarise(n = n(), n_positive = sum(Positive, na.rm = TRUE),
              pos_rate = mean(Positive, na.rm = TRUE) * 100, .groups = "drop")
  write_csv(pos_rate_summary, file.path(CONFIG$output_dir, "02_mcherry_pos_rate_by_age_zone.csv"))

  mc_area_summary <- mcherry %>%
    group_by(Age, Pos_label) %>%
    summarise(n = n(), mean_area = mean(Area), median_area = median(Area),
              sd_area = sd(Area), se_area = sd(Area) / sqrt(n()), .groups = "drop")
  write_csv(mc_area_summary, file.path(CONFIG$output_dir, "03_mcherry_area_pos_vs_neg.csv"))

  mc_int_summary <- mcherry %>%
    filter(Positive) %>%
    group_by(radial_bin) %>%
    summarise(n = n(), mean_intensity = mean(Mean), median_intensity = median(Mean),
              sd_intensity = sd(Mean), .groups = "drop")
  write_csv(mc_int_summary, file.path(CONFIG$output_dir, "04_mcherry_intensity_by_zone.csv"))

  quad_summary <- mcherry %>%
    group_by(quadrant) %>%
    summarise(n = n(), pos_rate = mean(Positive, na.rm = TRUE) * 100, .groups = "drop")
  write_csv(quad_summary, file.path(CONFIG$output_dir, "05_positivity_by_quadrant.csv"))

  sample_summary <- mcherry %>%
    group_by(Sample, Age, Animal) %>%
    summarise(n = n(), n_positive = sum(Positive, na.rm = TRUE),
              pos_rate = mean(Positive, na.rm = TRUE) * 100,
              mean_area = mean(Area), .groups = "drop") %>%
    arrange(Age, pos_rate)
  write_csv(sample_summary, file.path(CONFIG$output_dir, "06_per_sample_summary.csv"))
}

# 5h. DAPI nuclei by zone
if (has_dapi) {
  dapi_summary <- dapi %>%
    group_by(radial_bin) %>%
    summarise(n = n(), mean_nucleus_area = mean(Area), sd_nucleus_area = sd(Area),
              se_nucleus_area = sd(Area) / sqrt(n()), .groups = "drop")
  write_csv(dapi_summary, file.path(CONFIG$output_dir, "07_dapi_by_zone.csv"))
}

cat("  Summary tables written.\n\n")

# =============================================================================
# 6. PLOTS
# =============================================================================

cat("Generating plots...\n")

# ── Plot 1: Fiber area — mCherry positive vs negative by age ─────────────────
if (has_mcherry) {
  p1 <- mc_area_summary %>%
    ggplot(aes(x = Age, y = mean_area, fill = Pos_label)) +
    geom_col(position = position_dodge(0.7), width = 0.6) +
    geom_errorbar(aes(ymin = mean_area - se_area, ymax = mean_area + se_area),
                  position = position_dodge(0.7), width = 0.2, linewidth = 0.5) +
    scale_fill_manual(values = POS_COLORS, name = "mCherry") +
    labs(title    = "Fiber area: mCherry positive vs negative",
         subtitle = "Mean ± SE | Phalloidin-measured area",
         x = NULL, y = "Mean fiber area (µm²)") +
    theme_fiber()
  save_plot(p1, "01_area_pos_vs_neg.png", h = 5)
}

# ── Plot 2: Fiber area by radial zone × pos/neg ──────────────────────────────
if (has_mcherry) {
  area_radial <- mcherry %>%
    group_by(radial_bin, Pos_label) %>%
    summarise(mean_area = mean(Area), n = n(),
              se = sd(Area) / sqrt(n()), .groups = "drop")

  p2 <- area_radial %>%
    ggplot(aes(x = radial_bin, y = mean_area, fill = Pos_label)) +
    geom_col(position = position_dodge(0.7), width = 0.6) +
    geom_errorbar(aes(ymin = mean_area - se, ymax = mean_area + se),
                  position = position_dodge(0.7), width = 0.2, linewidth = 0.5) +
    scale_fill_manual(values = POS_COLORS, name = "mCherry") +
    labs(title = "Fiber area by radial zone — mCherry positive vs negative",
         x = "Radial zone (center → periphery)", y = "Mean fiber area (µm²)") +
    theme_fiber()
  save_plot(p2, "02_area_by_zone_pos_neg.png")
}

# ── Plot 3: mCherry intensity profile by zone (positive fibers only) ─────────
if (has_mcherry && nrow(mc_int_summary) > 0) {
  p3 <- mc_int_summary %>%
    ggplot(aes(x = radial_bin, y = mean_intensity, group = 1)) +
    geom_line(colour = "#D85A30", linewidth = 1.2) +
    geom_point(aes(size = n), colour = "#D85A30", fill = "#FAECE7",
               shape = 21, stroke = 1.2) +
    scale_size_continuous(range = c(2, 7), name = "Fiber count") +
    labs(title    = "mCherry fluorescence intensity by radial zone",
         subtitle = "Positive fibers only | Point size = fiber count in zone",
         x = "Radial zone (center → periphery)", y = "Mean mCherry intensity") +
    theme_fiber()
  save_plot(p3, "03_mcherry_intensity_profile.png", h = 5)
}

# ── Plot 4: mCherry positivity rate — age × radial zone ──────────────────────
if (has_mcherry) {
  p4 <- pos_rate_summary %>%
    ggplot(aes(x = radial_bin, y = pos_rate, fill = Age)) +
    geom_col(position = position_dodge(0.7), width = 0.6) +
    geom_hline(yintercept = 50, linetype = "dashed", colour = "grey60", linewidth = 0.5) +
    scale_fill_manual(values = AGE_COLORS) +
    scale_y_continuous(limits = c(0, 105), labels = scales::percent_format(scale = 1)) +
    labs(title    = "mCherry positivity rate by radial zone and age",
         subtitle = "% of fibers flagged mCherry-positive in each zone | Dashed = 50%",
         x = "Radial zone (center → periphery)", y = "Positivity rate (%)") +
    theme_fiber()
  save_plot(p4, "04_positivity_rate_by_zone_age.png")
}

# ── Plot 5: Phalloidin fiber area profile by age ──────────────────────────────
if (has_phal) {
  p5 <- phal_area_summary %>%
    ggplot(aes(x = radial_bin, y = mean_area, colour = Age, group = Age)) +
    geom_ribbon(aes(ymin = mean_area - se_area, ymax = mean_area + se_area,
                    fill = Age), alpha = 0.12, colour = NA) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) +
    scale_colour_manual(values = AGE_COLORS) +
    scale_fill_manual(values   = AGE_COLORS) +
    labs(title    = "Phalloidin fiber area profile by radial zone",
         subtitle = "All fibers | Ribbon = ±1 SE",
         x = "Radial zone (center → periphery)", y = "Mean fiber area (µm²)") +
    theme_fiber()
  save_plot(p5, "05_phal_area_profile.png")
}

# ── Plot 6: Spatial asymmetry — positivity by quadrant ───────────────────────
if (has_mcherry) {
  p6 <- quad_summary %>%
    mutate(quadrant = factor(quadrant, levels = names(QUAD_COLORS))) %>%
    ggplot(aes(x = quadrant, y = pos_rate, fill = quadrant)) +
    geom_col(width = 0.6) +
    geom_hline(yintercept = mean(quad_summary$pos_rate),
               linetype = "dashed", colour = "grey50", linewidth = 0.6) +
    scale_fill_manual(values = QUAD_COLORS, guide = "none") +
    scale_y_continuous(labels = scales::percent_format(scale = 1)) +
    labs(title    = "Spatial asymmetry: mCherry positivity by quadrant",
         subtitle = "Dashed line = overall mean positivity rate",
         x = NULL, y = "Positivity rate (%)") +
    theme_fiber()
  save_plot(p6, "06_spatial_asymmetry.png", h = 5)
}

# ── Plot 7: Per-sample positivity rates ──────────────────────────────────────
if (has_mcherry && nrow(sample_summary) > 0) {
  p7 <- sample_summary %>%
    mutate(Sample = fct_reorder(Sample, pos_rate)) %>%
    ggplot(aes(x = pos_rate, y = Sample, fill = Age)) +
    geom_col() +
    geom_vline(xintercept = 50, linetype = "dashed", colour = "grey50", linewidth = 0.5) +
    scale_fill_manual(values = AGE_COLORS) +
    scale_x_continuous(labels = scales::percent_format(scale = 1)) +
    facet_wrap(~Age, scales = "free_y", ncol = 2) +
    labs(title    = "Per-sample mCherry positivity rate",
         subtitle = "Sorted low to high within each age group | Dashed = 50%",
         x = "Positivity rate (%)", y = NULL) +
    theme_fiber() +
    theme(legend.position = "none")
  save_plot(p7, "07_per_sample_positivity.png",
            h = max(5, nrow(sample_summary) * 0.25 + 2))
}

# ── Plot 8: Area vs intensity scatter (mCherry positive, normalised) ─────────
if (has_mcherry) {
  pos_fibers <- mcherry %>% filter(Positive)
  if (nrow(pos_fibers) > 0) {
    n_scatter    <- min(5000L, nrow(pos_fibers))
    scatter_data <- pos_fibers %>% slice_sample(n = n_scatter)
    p8 <- scatter_data %>%
      ggplot(aes(x = Area_norm, y = Mean_norm, colour = Age)) +
      geom_point(alpha = 0.25, size = 0.9) +
      geom_smooth(method = "lm", se = TRUE, linewidth = 1.1, alpha = 0.2) +
      scale_colour_manual(values = AGE_COLORS) +
      facet_wrap(~Age) +
      labs(title    = "Normalised fiber area vs mCherry intensity",
           subtitle = "Positive fibers only | Per-sample min-max normalisation | Line = linear fit",
           x = "Normalised fiber area (0–1)", y = "Normalised mCherry intensity (0–1)",
           caption  = paste0("n = ", nrow(scatter_data), " sampled positive fibers")) +
      theme_fiber() +
      theme(legend.position = "none")
    save_plot(p8, "08_area_vs_intensity_scatter.png")
  }
}

# ── Plot 9: DAPI nuclei count and size by zone ────────────────────────────────
if (has_dapi) {
  p9a <- dapi_summary %>%
    ggplot(aes(x = radial_bin, y = n, fill = radial_bin)) +
    geom_col() +
    scale_fill_manual(values = ring_pal, guide = "none") +
    scale_y_continuous(labels = scales::comma) +
    labs(title = "DAPI nuclei count by radial zone",
         x = "Radial zone", y = "Nucleus count") +
    theme_fiber()

  p9b <- dapi_summary %>%
    ggplot(aes(x = radial_bin, y = mean_nucleus_area, fill = radial_bin)) +
    geom_col() +
    geom_errorbar(aes(ymin = mean_nucleus_area - se_nucleus_area,
                      ymax = mean_nucleus_area + se_nucleus_area),
                  width = 0.25, linewidth = 0.5) +
    scale_fill_manual(values = ring_pal, guide = "none") +
    labs(title    = "Mean nucleus area by radial zone",
         subtitle = "Error bars = ±1 SE",
         x = "Radial zone", y = "Mean nucleus area (µm²)") +
    theme_fiber()

  p9 <- p9a + p9b +
    plot_annotation(title = "DAPI nuclear distribution",
                    theme = theme(plot.title = element_text(size = 13, face = "bold")))
  save_plot(p9, "09_dapi_nuclei.png")
}

# ── Plot 10: Phalloidin spatial scatter by sample ────────────────────────────
if (has_phal) {
  scatter_phal <- phal %>%
    group_by(Sample) %>%
    slice_sample(prop = 1) %>%   # shuffle within group
    slice_head(n = 300) %>%      # then cap at 300 (keeps all if group < 300)
    ungroup()

  n_samples <- n_distinct(scatter_phal$Sample)
  p10 <- scatter_phal %>%
    ggplot(aes(x = X_scaled, y = Y_scaled, colour = radial_bin)) +
    geom_point(size = 0.5, alpha = 0.4) +
    scale_colour_manual(values = ring_pal, name = "Zone") +
    coord_fixed() +
    facet_wrap(~Sample, ncol = 6) +
    labs(title    = "Phalloidin fiber positions — spatial scatter by sample",
         subtitle = "Coloured by radial zone | Scaled coordinates",
         x = "X scaled", y = "Y scaled") +
    theme_fiber() +
    theme(axis.text  = element_text(size = 6),
          strip.text = element_text(size = 7))
  save_plot(p10, "10_spatial_scatter_by_sample.png",
            w = 14, h = max(6, ceiling(n_samples / 6) * 2.5))
}

# ── Plot 11: Growth type diagnostic (fiber size vs positivity by zone) ────────
if (has_phal && has_mcherry) {
  diag_data <- pos_rate_summary %>%
    left_join(phal_area_summary %>% select(Age, radial_bin, mean_area),
              by = c("Age", "radial_bin"))

  if (CONFIG$tissue_type == "tail") {
    edge_zones <- tail(zone_labels, 2)
    diag_data  <- diag_data %>%
      mutate(zone_type = if_else(radial_bin %in% edge_zones,
                                 "Peripheral (SHP candidate)", "Interior (HT/MHP)"))
    p11 <- diag_data %>%
      ggplot(aes(x = mean_area, y = pos_rate, colour = Age,
                 shape = zone_type, label = radial_bin)) +
      geom_path(aes(group = Age), linewidth = 0.6, linetype = "dashed", alpha = 0.5) +
      geom_point(size = 4) +
      geom_text(size = 3, vjust = -0.8, show.legend = FALSE) +
      scale_shape_manual(values = c("Peripheral (SHP candidate)" = 17,
                                    "Interior (HT/MHP)"          = 16),
                         name = "Zone type")
  } else {
    p11 <- diag_data %>%
      ggplot(aes(x = mean_area, y = pos_rate, colour = Age, label = radial_bin)) +
      geom_path(aes(group = Age), linewidth = 0.6, linetype = "dashed", alpha = 0.5) +
      geom_point(size = 4) +
      geom_text(size = 3, vjust = -0.8, show.legend = FALSE)
  }

  p11 <- p11 +
    scale_colour_manual(values = AGE_COLORS) +
    scale_y_continuous(labels = scales::percent_format(scale = 1)) +
    labs(title    = paste(tools::toTitleCase(CONFIG$tissue_type),
                          "tissue: fiber size vs mCherry positivity by zone"),
         subtitle = "Each point = one radial zone | Path connects zones Core → Far",
         x = "Mean phalloidin fiber area (µm²)", y = "mCherry positivity rate (%)") +
    theme_fiber()
  save_plot(p11, "11_growth_type_diagnostic.png")
}

# =============================================================================
# 7. SESSION SUMMARY
# =============================================================================

cat("\n=== Analysis complete ===\n")
cat("Output directory:", CONFIG$output_dir, "\n\n")
cat("Files written:\n")
list.files(CONFIG$output_dir) %>% paste0("  ", .) %>% cat(sep = "\n")
cat("\n")

cat("Key counts:\n")
if (has_phal)    cat("  Phalloidin fibers:", nrow(phal),    "\n")
if (has_mcherry) cat("  mCherry fibers   :", nrow(mcherry), "\n")
if (has_dapi)    cat("  DAPI nuclei      :", nrow(dapi),    "\n")
if (has_mcherry) cat("  Samples          :", n_distinct(mcherry$Sample), "\n")

if (has_mcherry && exists("pos_rate_summary") && nrow(pos_rate_summary) > 0) {
  for (ag in ages) {
    ag_rates <- pos_rate_summary %>% filter(Age == ag) %>% pull(pos_rate)
    if (length(ag_rates) > 0) {
      cat(sprintf("  %s positivity: %.1f%% (mean across zones)\n", ag, mean(ag_rates, na.rm = TRUE)))
    }
  }
}

if (has_mcherry && exists("quad_summary") && nrow(quad_summary) > 0) {
  quad_range <- diff(range(quad_summary$pos_rate))
  cat(sprintf("  Quadrant asymmetry: %.1f%% range across quadrants\n", quad_range))
  if (quad_range > 10) cat("  ** Asymmetry > 10% — consider anisotropy analysis **\n")
}

# Axolotl Muscle Fiber Growth Analysis

**Quantifying hypertrophy and hyperplasia in axolotl skeletal muscle using fluorescence microscopy cross-sections**

Collaborative research project with the [Mount Desert Island Biological Laboratory (MDIBL)](https://mdibl.org/) Computational & Genomic Data Science (CGDS) Core and the Murawala Research Group.

---

## Project Overview

This project investigates whether muscle fiber growth in axolotls (*Ambystoma mexicanum*) more closely resembles the mammalian or zebrafish model:

| Species | Growth Mechanisms |
|---------|-------------------|
| **Mammals (mice)** | Hypertrophy only |
| **Zebrafish** | Hypertrophy + mosaic hyperplasia + stratified hyperplasia |
| **Axolotl (tail)** | All three mechanisms — under investigation |
| **Axolotl (limb)** | Predominantly hypertrophy, minor mosaic hyperplasia |

The goal is to **systematically differentiate and quantify three types of muscle growth** from fluorescence microscopy cross-sections:

- **HT — Hypertrophy**: growth of existing fibers via fusion into multi-nucleated cells
- **MHP — Mosaic Hyperplasia**: new fibers forming and integrating *within* a muscle bundle
- **SHP — Stratified Hyperplasia**: new fibers forming along the *periphery* of a muscle bundle

---

## Biological Background

```
Progenitor cells → Satellite cells → Mature muscle fiber
                ↘                  ↗
                 (direct conversion)
```

Progenitor cells can signal satellite cells to fuse with existing fibers (hypertrophy), or convert directly into satellite cells that eventually become new fibers (hyperplasia). Tail samples contain all three cell types; limb samples contain only satellite cells and mature fibers.

---

## Data

Samples are fluorescence microscopy cross-sections of axolotl **tail** or **limb** tissue, collected at day 0 labeling and after 60 days of muscle growth. Four imaging channels are used:

| Channel | Label | Purpose |
|---------|-------|---------|
| mCherry | Cre-Lox targeted cells | Identifies labeled cell populations at day 0 |
| Phalloidin (C1) | All muscle fibers | Fiber geometry: position, area, perimeter |
| DAPI | All nuclei | Counts nuclei per fiber cross-section |
| GFP | — | Not used in current analysis |

**Dataset:** 164,185 Phalloidin-segmented fibers and 160,435 mCherry-channel measurements across 46 samples (Old and Young animals), segmented by the MDIBL microscopy core and provided as CSV files.

---

## Analysis Pipeline

The pipeline runs in eight stages, each handled by a dedicated script. Intermediate outputs are passed between stages via CSV files.

```
1. initial_visualization.R          Raw data QC plots
         ↓
2. optimal_ratio_search.R           Hull ratio sweep & elbow detection (development)
   ratio_comparison.R               Compare elbow vs plateau methods (development)
         ↓
3. alignment_boundary_detect_v6.R   Outlier removal + final concave hull fitting
         ↓
4. dt_pca_v6.R                      Distance-transform PCA for orientation axes
   pca_axis_overlay_v6.R            PCA axis visualization
         ↓
5. alignment_v6.R                   Rotate & center samples to common frame
         ↓
6. sign_fix_v6.R                    Disambiguate PCA axis direction (mass-right convention)
         ↓
7. muscle_fiber_analysis.R          Full downstream analysis: zoning, positivity,
                                    area profiles, DAPI, spatial asymmetry, diagnostics
```

---

### Stage 1 — Initial Visualization (`initial_visualization.R`)

Quick QC plot of raw Phalloidin fiber positions across all samples to check data completeness and identify outlier samples before any processing.

---

### Stage 2 — Hull Ratio Development (`optimal_ratio_search.R`, `ratio_comparison.R`)

Development scripts used to determine the best method for selecting concave hull ratio parameters. `optimal_ratio_search.R` implements a perimeter sweep (ratio 0.6–1.0) and second-derivative elbow detection, with boundary-sample extension for samples whose optima fell near the sweep edges. `ratio_comparison.R` compares the elbow method against a plateau method side-by-side, producing diagnostic plots that informed the final choice of elbow detection for Stage 3.

---

### Stage 3 — Boundary Detection & Outlier Removal (`alignment_boundary_detect_v6.R`)

The core tissue boundary pipeline. For each sample, three complementary per-fiber outlier metrics are computed before hull fitting:

- **`mean_dist`** — mean Euclidean distance to all other fibers; catches globally isolated fibers
- **`local_linearity`** — ratio of 1st to 2nd eigenvalue from PCA on k nearest neighbors; high ratio indicates peninsula/spur arrangements
- **`angular_coverage`** — fraction of angular sectors around each fiber containing at least one neighbor; low coverage indicates edge-only neighbors

All three metrics are z-scored within each sample and combined into a composite outlier flag. After removing flagged fibers, a concave hull is fit to the cleaned point cloud using a **second-derivative elbow method** on a perimeter sweep (ratio 0.5–1.0, 50 steps). Fibers are then classified as within- or outside-hull using `sf::st_intersects()`.

![Hull Overlay](figures/hull_overlay_final.png)
*Concave hull boundaries across Old-tissue samples. Blue = boundary, green = within-hull fibers.*

**Outputs:** `phalloidin_clean_with_hull_v6.csv`, `sample_optimal_ratios_v6.csv`

---

### Stage 4 — PCA Orientation (`dt_pca_v6.R`, `pca_axis_overlay_v6.R`)

Determines the principal orientation axis of each tissue cross-section so that all samples can be rotated to a common frame. Rather than running PCA directly on fiber positions (which can be biased by fiber density), this stage uses a **distance-transform PCA**: the hull boundary is resampled to evenly spaced points, and PCA is run on those boundary coordinates. This gives an orientation axis driven by tissue shape rather than internal fiber distribution. `pca_axis_overlay_v6.R` overlays the resulting axes on hull plots for visual verification.

**Output:** `dt_pca_results_v6.csv`

---

### Stage 5 — Alignment (`alignment_v6.R`)

Rotates and centers each sample's fiber coordinates using the PC1 angle from Stage 4. PC1 angles are normalized to [0°, 180°) to resolve the bidirectionality of PCA axes, then a rotation matrix is applied to bring PC1 to horizontal. Hull boundaries are transformed alongside the fiber coordinates. A verification PCA on aligned coordinates confirms successful rotation.

**Output:** `phalloidin_aligned_v6.csv`

---

### Stage 6 — Sign Disambiguation (`sign_fix_v6.R`)

After rotation, PCA axes remain 180°-ambiguous (pointing left or right). This script resolves that by applying a **"more mass to the right" convention**: for each sample, if the mean aligned X coordinate is negative (more fibers on the left), the sample is flipped via a 180° rotation (negating both X and Y). This ensures consistent left-right orientation across all samples before spatial analysis.

**Output:** `phalloidin_signed_v6.csv`

---

### Stage 7 — Downstream Analysis (`muscle_fiber_analysis.R`)

A configurable analysis script that takes the aligned, normalized fiber data and produces all summary tables and figures. Key CONFIG parameters are tissue type (`"limb"` or `"tail"`) and radial zone method (`equal_width` or `quantile`). The script:

- Splits data by channel (Phalloidin, mCherry, DAPI)
- Assigns fibers to **radial zones** (Core → Inner → Mid → Outer → Edge → Far) and **quadrants** (NE/NW/SW/SE)
- Normalizes coordinates and intensity values per sample
- Produces all summary CSVs and figures (see Results below)

---

## Results

### Fiber Area Profiles

**Overall fiber area by age and mCherry status:**

![Area: Positive vs Negative](figures/01_area_pos_vs_neg.png)

mCherry-positive fibers in Young animals are substantially larger than all other groups, suggesting pronounced hypertrophy of labeled fibers in younger animals.

**Fiber area by radial zone:**

![Area by Zone](figures/02_area_by_zone_pos_neg.png)

Fiber area peaks in the Mid zone and declines toward the periphery in both age groups, consistent with larger, more mature fibers residing in the tissue interior.

![Area Profile](figures/05_phal_area_profile.png)

The radial area profile differs between Old (purple) and Young (green) animals, with Young animals showing a sharper central peak in the Core zone.

---

### mCherry Positivity & Spatial Distribution

**Positivity rate by radial zone:**

![Positivity Rate by Zone](figures/04_positivity_rate_by_zone_age.png)

Old animals (purple) show high mCherry positivity across all zones (63–83%), while Young animals (green) show near-zero rates throughout. This age-dependent pattern is a key signature distinguishing growth mechanisms between the two groups.

**mCherry intensity profile:**

![mCherry Intensity](figures/03_mcherry_intensity_profile.png)

Among positive fibers only, mCherry intensity increases steeply from Core toward the periphery. Peripheral fibers carry stronger original labeling signal — consistent with stratified hyperplasia adding new labeled fibers at the tissue margins.

**Spatial asymmetry by quadrant:**

![Spatial Asymmetry](figures/06_spatial_asymmetry.png)

NW (Q2) and SW (Q3) quadrants show slightly elevated positivity rates, pointing to mild directional asymmetry in growth activity — consistent with expected dorsal-ventral differences in tail samples.

**Per-sample positivity:**

![Per-sample Positivity](figures/07_per_sample_positivity.png)

Old samples span a wide range of positivity rates (~10–100%); Young samples are consistently low (<30%).

---

### Fiber Area vs. mCherry Intensity

![Area vs. Intensity](figures/08_area_vs_intensity_scatter.png)

In Young animals, larger fibers tend to have lower mCherry intensity (negative correlation), consistent with label dilution through hypertrophic growth. In Old animals, intensity is relatively flat across fiber sizes.

---

### DAPI Nuclear Distribution

![DAPI Nuclei](figures/09_dapi_nuclei.png)

Nucleus count increases from Core to Edge/Far zones, and nucleus area is largest at the Edge zone (~60 µm²). Peripheral nuclear accumulation supports the presence of stratified hyperplasia at tissue margins.

---

### Spatial Scatter by Sample

![Spatial Scatter](figures/10_spatial_scatter_by_sample.png)

PCA-aligned, scaled fiber coordinates across all 46 samples, coloured by radial zone. Consistent spatial structure across samples validates the alignment and normalization approach.

---

### Growth Type Diagnostic

![Growth Type Diagnostic](figures/11_growth_type_diagnostic.png)

Each point is a radial zone; the path connects Core → Far. Old animals (purple) follow a non-linear trajectory — high positivity at both small (Core) and large (Mid) fibers — consistent with co-occurring hypertrophy and mosaic hyperplasia. Young animals (green) show uniformly low positivity regardless of zone or fiber size.

---

## Repository Structure

```
axolotl-muscle-fiber-analysis/
├── README.md
├── scripts/
│   ├── initial_visualization.R            # Stage 1: raw data QC
│   ├── optimal_ratio_search.R             # Stage 2: hull ratio sweep (development)
│   ├── ratio_comparison.R                 # Stage 2: elbow vs plateau comparison
│   ├── alignment_boundary_detect_v6.R     # Stage 3: outlier removal + hull fitting
│   ├── dt_pca_v6.R                        # Stage 4: distance-transform PCA
│   ├── pca_axis_overlay_v6.R              # Stage 4: PCA axis visualization
│   ├── alignment_v6.R                     # Stage 5: rotation & centering
│   ├── sign_fix_v6.R                      # Stage 6: sign disambiguation
│   └── muscle_fiber_analysis.R            # Stage 7: full downstream analysis
├── docs/
│   ├── project_description.md             # Biological background and aims
│   └── methods_notes.md                   # Algorithm design decisions
├── figures/
│   ├── hull_overlay_final.png
│   ├── 01_area_pos_vs_neg.png
│   ├── 02_area_by_zone_pos_neg.png
│   ├── 03_mcherry_intensity_profile.png
│   ├── 04_positivity_rate_by_zone_age.png
│   ├── 05_phal_area_profile.png
│   ├── 06_spatial_asymmetry.png
│   ├── 07_per_sample_positivity.png
│   ├── 08_area_vs_intensity_scatter.png
│   ├── 09_dapi_nuclei.png
│   ├── 10_spatial_scatter_by_sample.png
│   └── 11_growth_type_diagnostic.png
└── results/
    ├── 00_summary_counts.csv
    ├── 01_phal_area_by_age_zone.csv
    ├── 02_mcherry_pos_rate_by_age_zone.csv
    ├── 03_mcherry_area_pos_vs_neg.csv
    ├── 04_mcherry_intensity_by_zone.csv
    ├── 05_positivity_by_quadrant.csv
    ├── 06_per_sample_summary.csv
    └── 07_dapi_by_zone.csv
```

---

## Dependencies

**R packages:** `tidyverse`, `sf`, `FNN`, `patchwork`

**Input data:** `COMPLETE_COMBINED_MUSCLE_FIBER_DATA.csv` — per-fiber measurements with `Sample`, `Age`, `Channel_Code`, `X`, `Y`, area, perimeter, and intensity columns. Not included in this repository (available from MDIBL CGDS Core).

---

## Key Findings

- **Age strongly predicts mCherry positivity**: Old animals retain high labeling rates (63–83% across zones); Young animals are consistently near-zero
- **Peripheral enrichment of new fibers**: Rising mCherry intensity and DAPI accumulation toward outer zones are consistent with stratified hyperplasia at tissue margins
- **Fiber size–intensity anticorrelation in Young animals**: Consistent with label dilution through hypertrophic growth
- **Mild spatial asymmetry**: NW/SW quadrants show elevated positivity, consistent with directional growth bias expected in tail samples
- **Growth type diagnostic**: Zone-level area vs. positivity trajectories distinguish Old and Young tissue and are consistent with co-occurring hypertrophy and mosaic hyperplasia in Old animals

---

## Collaborators

- **MDIBL CGDS Core** — image segmentation, data provision, analysis collaboration
- **Murawala Research Group** — biological design, experimental execution, interpretation

---

## Status

Analysis complete. Manuscript in preparation.

# ============================================================
# 03_4_4_CD8_subtypes_inside_outside_distance_v1_2.R
#
# Goal
# ----
# User supplies ONE distance cutoff (um), e.g. 30.
# For every FOV:
#   Group 1 = distance <= cutoff
#   Group 2 = distance > cutoff
#
# Analyses are performed for:
#   1) All FOVs together (ID + IE + IN)
#   2) ID only
#   3) IE only
#   4) IN only
#
# Outputs
# -------
# A. CD8 subtype composition pie charts
#    - two pies per scope: <= cutoff vs > cutoff
#    - default pie summary is FOV-balanced mean composition
#    - pooled-cell composition is also exported to CSV
#
# B. Violin plots for each CD8 subtype, comparing the two distance groups
#    Set 1 denominator:
#       subtype cells / all CD8 cells in the same FOV x distance group
#    Set 2 denominator:
#       subtype cells / all non-tumor cells in the same FOV x distance group
#
# C. Paired statistics
#    Because the two distance groups come from the SAME FOV/patient,
#    the default test is a two-sided paired Wilcoxon signed-rank test.
#    Raw p values are shown on the figures.
#    BH-adjusted p values across the 5 CD8 subtypes within each
#    scope x denominator are exported to CSV.
#
# Important project conventions
# -----------------------------
# - disttonearesttumor is already in um in the current compatible metadata.
#   DO NOT multiply by 1000.
# - One patient = one FOV.
# - immunotype_BW_top_40 is the current ID / IE / IN classification.
# - Fine CD8 labels come from celltype_new.
# - Tumor cells are excluded from the "all non-tumor" denominator.
#
# Immunotype colors remain fixed project-wide:
#   ID = #000000
#   IE = #0000FF
#   IN = #FF0000
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ============================================================
# 01 paths
# ============================================================
find_project_root <- function(start = getwd()) {
  current <- tryCatch(
    normalizePath(start, winslash = "/", mustWork = TRUE),
    error = function(e) NA_character_
  )
  if (is.na(current)) stop("Cannot resolve the current working directory.")
  
  repeat {
    if (file.exists(file.path(current, "TNBC_TMA_6k_CosMx.Rproj"))) return(current)
    parent <- dirname(current)
    if (identical(parent, current)) break
    current <- parent
  }
  
  stop(
    "Cannot find TNBC_TMA_6k_CosMx.Rproj. ",
    "Open the R project or run the script from the repository root/subdirectory."
  )
}

project_dir <- find_project_root()

meta_file <- file.path(
  project_dir,
  "results/03_2/data/03_2_meta_with_immunotype_BW_top_40.rds"
)

out_root <- file.path(project_dir, "results/03_4_4")
data_dir <- file.path(out_root, "data")
fig_dir  <- file.path(out_root, "figures")

dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir,  recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 02 USER PARAMETERS
# ============================================================
# Main parameter to change.
distance_cutoff_um <- 40

# For subtype / all-CD8 proportions.
min_cd8_cells_per_FOV_group <- 3L

# For subtype / all-nontumor proportions.
min_nontumor_cells_per_FOV_group <- 1L

# Pie summary:
#   "FOV_balanced" = calculate composition in each FOV first,
#                    then average proportions across FOVs.
#   "pooled_cells" = pool all CD8 cells across FOVs.
pie_summary_mode <- "FOV_balanced"

# Paired comparison between <=cutoff and >cutoff within the same FOV.
# Options: "wilcoxon", "paired_t"
paired_test <- "wilcoxon"

show_paired_lines <- FALSE

# Relatively tall figures.
pie_pdf_width <- 9.5
pie_pdf_height <- 6.5
violin_pdf_width <- 9.5
violin_pdf_height <- 11.0
violin_ncol <- 2L

# Hide labels for very small pie wedges, matching 03_4_2.
label_min_pct <- 4

# ============================================================
# 03 labels and colors
# ============================================================
immunotype_order <- c("ID", "IE", "IN")

immunotype_cols <- c(
  "ID" = "#000000",
  "IE" = "#0000FF",
  "IN" = "#FF0000"
)

cd8_subtypes <- c(
  "cd8_naive",
  "cd8_tcm",
  "cd8_tem",
  "cd8_cytotoxic",
  "cd8_exhausted"
)

cd8_subtype_labels <- c(
  "cd8_naive"     = "Naive",
  "cd8_tcm"       = "TCM",
  "cd8_tem"       = "TEM",
  "cd8_cytotoxic" = "Cytotoxic",
  "cd8_exhausted" = "Exhausted"
)

cd8_subtype_cols <- c(
  "cd8_naive"     = "#2166AC",
  "cd8_tcm"       = "#1B9E77",
  "cd8_tem"       = "#E6AB02",
  "cd8_cytotoxic" = "#D73027",
  "cd8_exhausted" = "#7B3294"
)

scope_order <- c("All_FOVs", "ID", "IE", "IN")

scope_labels <- c(
  "All_FOVs" = "All FOVs (ID + IE + IN)",
  "ID" = "ID",
  "IE" = "IE",
  "IN" = "IN"
)

distance_cols <- c(
  "Within" = "#4D4D4D",
  "Beyond" = "#D9D9D9"
)

# ============================================================
# 04 helpers
# ============================================================
safe_tag <- function(x) gsub("[^A-Za-z0-9_]+", "_", x)

save_csv <- function(x, filename) {
  out <- file.path(data_dir, paste0(filename, ".csv"))
  fwrite(as.data.table(x), out)
  message("Saved: ", out)
  invisible(out)
}

save_rds <- function(x, filename) {
  out <- file.path(data_dir, paste0(filename, ".rds"))
  saveRDS(x, out)
  message("Saved: ", out)
  invisible(out)
}

format_p <- function(p) {
  if (length(p) == 0 || is.na(p) || !is.finite(p)) return("NA")
  if (p < 1e-4) return(formatC(p, format = "e", digits = 1))
  if (p < 0.001) return(formatC(p, format = "e", digits = 2))
  formatC(p, format = "f", digits = 3)
}

p_to_star <- function(p) {
  if (length(p) == 0 || is.na(p) || !is.finite(p)) return("NA")
  if (p < 1e-4) return("****")
  if (p < 1e-3) return("***")
  if (p < 1e-2) return("**")
  if (p < 0.05) return("*")
  "ns"
}

run_paired_test <- function(x, y, test = c("wilcoxon", "paired_t")) {
  test <- match.arg(test)
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  
  if (length(x) < 3L) return(NA_real_)
  
  if (test == "wilcoxon") {
    out <- tryCatch(
      suppressWarnings(
        stats::wilcox.test(
          x = x,
          y = y,
          paired = TRUE,
          alternative = "two.sided",
          exact = FALSE
        )
      ),
      error = function(e) NULL
    )
  } else {
    out <- tryCatch(
      stats::t.test(
        x = x,
        y = y,
        paired = TRUE,
        alternative = "two.sided"
      ),
      error = function(e) NULL
    )
  }
  
  if (is.null(out)) return(NA_real_)
  as.numeric(out$p.value)
}

# ============================================================
# 05 input checks
# ============================================================
if (!file.exists(meta_file)) stop("Metadata file not found: ", meta_file)

if (!is.numeric(distance_cutoff_um) ||
    length(distance_cutoff_um) != 1L ||
    !is.finite(distance_cutoff_um) ||
    distance_cutoff_um < 0) {
  stop("distance_cutoff_um must be one finite non-negative number.")
}

if (!pie_summary_mode %in% c("FOV_balanced", "pooled_cells")) {
  stop('pie_summary_mode must be "FOV_balanced" or "pooled_cells".')
}

if (!paired_test %in% c("wilcoxon", "paired_t")) {
  stop('paired_test must be "wilcoxon" or "paired_t".')
}

# ============================================================
# 06 read metadata
# ============================================================
meta <- readRDS(meta_file)
meta_dt <- as.data.table(meta)

required_cols <- c(
  "FOV",
  "Run_Tissue_name",
  "BC6k_refined",
  "celltype_new",
  "disttonearesttumor",
  "immunotype_BW_top_40"
)

missing_cols <- setdiff(required_cols, names(meta_dt))
if (length(missing_cols) > 0) {
  stop("Missing required metadata columns: ", paste(missing_cols, collapse = ", "))
}

meta_dt[, FOV := as.character(FOV)]
meta_dt[, Run_Tissue_name := as.character(Run_Tissue_name)]
meta_dt[, BC6k_refined := as.character(BC6k_refined)]
meta_dt[, celltype_new := as.character(celltype_new)]
meta_dt[, immunotype_BW_top_40 := as.character(immunotype_BW_top_40)]

# Current project convention: already in um.
meta_dt[, dist_to_tumor_um := as.numeric(disttonearesttumor)]

meta_dt <- meta_dt[
  immunotype_BW_top_40 %in% immunotype_order &
    is.finite(dist_to_tumor_um)
]

meta_dt[, immunotype_BW_top_40 := factor(
  immunotype_BW_top_40,
  levels = immunotype_order
)]

message("Cells retained with ID/IE/IN and finite tumor distance: ", nrow(meta_dt))
message("FOVs retained: ", uniqueN(meta_dt$FOV))

# ============================================================
# 07 distance groups
# ============================================================
within_label <- paste0("\u2264", distance_cutoff_um, " \u00b5m")
beyond_label <- paste0(">", distance_cutoff_um, " \u00b5m")

meta_dt[, distance_group := fifelse(
  dist_to_tumor_um <= distance_cutoff_um,
  "Within",
  "Beyond"
)]

meta_dt[, distance_group := factor(
  distance_group,
  levels = c("Within", "Beyond")
)]

# ============================================================
# 08 define non-tumor and CD8 cells
# ============================================================
meta_dt[, is_nontumor := BC6k_refined != "tumor"]
meta_dt[, is_cd8 := BC6k_refined == "T.cell.CD8"]

meta_dt[, cd8_subtype := fifelse(
  is_cd8 & celltype_new %in% cd8_subtypes,
  celltype_new,
  NA_character_
)]

cd8_qc <- meta_dt[
  is_cd8 == TRUE,
  .(
    total_broad_CD8 = .N,
    classified_into_5_subtypes = sum(!is.na(cd8_subtype)),
    unclassified_broad_CD8 = sum(is.na(cd8_subtype))
  )
]

save_csv(cd8_qc, "03_4_4_CD8_classification_QC")

if (cd8_qc$unclassified_broad_CD8 > 0) {
  warning(
    "There are ",
    cd8_qc$unclassified_broad_CD8,
    " broad T.cell.CD8 cells not assigned to one of the five requested ",
    "celltype_new CD8 subtypes. They remain in the total-CD8 denominator."
  )
}

# ============================================================
# 09 FOV annotation
# ============================================================
fov_anno <- unique(meta_dt[, .(
  FOV,
  Run_Tissue_name,
  immunotype_BW_top_40
)])

dup_fov <- fov_anno[, .N, by = FOV][N > 1]
if (nrow(dup_fov) > 0) {
  stop(
    "At least one FOV has more than one immunotype/slide annotation. ",
    "Please inspect metadata before running 03_4_4."
  )
}

save_csv(fov_anno, "03_4_4_FOV_annotation")

# ============================================================
# 10 count cells by FOV x distance group
# ============================================================
base_counts <- meta_dt[
  is_nontumor == TRUE,
  .(
    n_nontumor = .N,
    n_cd8_total = sum(is_cd8, na.rm = TRUE)
  ),
  by = .(
    FOV,
    Run_Tissue_name,
    immunotype_BW_top_40,
    distance_group
  )
]

subtype_counts <- meta_dt[
  is_nontumor == TRUE & !is.na(cd8_subtype),
  .N,
  by = .(
    FOV,
    Run_Tissue_name,
    immunotype_BW_top_40,
    distance_group,
    cd8_subtype
  )
]
setnames(subtype_counts, "N", "n_subtype")

# Build complete FOV x distance x subtype grid only for observed FOV-distance bins.
observed_bins <- unique(base_counts[, .(
  FOV,
  distance_group,
  Run_Tissue_name,
  immunotype_BW_top_40,
  n_nontumor,
  n_cd8_total
)])

subtype_grid <- observed_bins[
  ,
  .(cd8_subtype = cd8_subtypes),
  by = .(
    FOV,
    distance_group,
    Run_Tissue_name,
    immunotype_BW_top_40,
    n_nontumor,
    n_cd8_total
  )
]

subtype_long <- merge(
  subtype_grid,
  subtype_counts[, .(
    FOV,
    distance_group,
    cd8_subtype,
    n_subtype
  )],
  by = c("FOV", "distance_group", "cd8_subtype"),
  all.x = TRUE
)

subtype_long[is.na(n_subtype), n_subtype := 0L]

subtype_long[, cd8_subtype := factor(
  cd8_subtype,
  levels = cd8_subtypes
)]

subtype_long[, valid_cd8_denom :=
               n_cd8_total >= min_cd8_cells_per_FOV_group]

subtype_long[, valid_nontumor_denom :=
               n_nontumor >= min_nontumor_cells_per_FOV_group]

subtype_long[
  valid_cd8_denom == TRUE,
  prop_within_CD8 := n_subtype / n_cd8_total
]

subtype_long[
  valid_nontumor_denom == TRUE,
  prop_within_nontumor := n_subtype / n_nontumor
]

subtype_long[, pct_within_CD8 := 100 * prop_within_CD8]
subtype_long[, pct_within_nontumor := 100 * prop_within_nontumor]

save_csv(
  subtype_long,
  "03_4_4_FOV_distance_CD8_subtype_counts_and_proportions"
)

save_rds(
  subtype_long,
  "03_4_4_FOV_distance_CD8_subtype_counts_and_proportions"
)

# ============================================================
# 11 QC by immunotype x distance group
# ============================================================
distance_qc <- base_counts[
  ,
  .(
    n_FOV = uniqueN(FOV),
    median_n_nontumor = as.numeric(median(as.numeric(n_nontumor))),
    median_n_CD8 = as.numeric(median(as.numeric(n_cd8_total))),
    n_FOV_CD8_ge_min = uniqueN(
      FOV[n_cd8_total >= min_cd8_cells_per_FOV_group]
    ),
    n_FOV_nontumor_ge_min = uniqueN(
      FOV[n_nontumor >= min_nontumor_cells_per_FOV_group]
    )
  ),
  by = .(
    immunotype_BW_top_40,
    distance_group
  )
]

save_csv(distance_qc, "03_4_4_distance_group_QC_by_immunotype")

# ============================================================
# 12 add analysis scopes
# ============================================================
all_scope <- copy(subtype_long)
all_scope[, scope := "All_FOVs"]

type_scope <- copy(subtype_long)
type_scope[, scope := as.character(immunotype_BW_top_40)]

analysis_long <- rbindlist(
  list(all_scope, type_scope),
  use.names = TRUE,
  fill = TRUE
)

analysis_long <- analysis_long[scope %in% scope_order]
analysis_long[, scope := factor(scope, levels = scope_order)]

# ============================================================
# 13 pie summaries
# ============================================================
# FOV-balanced composition.
pie_fov_balanced <- analysis_long[
  valid_cd8_denom == TRUE & is.finite(prop_within_CD8),
  .(
    mean_prop = mean(prop_within_CD8, na.rm = TRUE),
    median_prop = as.numeric(median(as.numeric(prop_within_CD8), na.rm = TRUE)),
    n_FOV = uniqueN(FOV)
  ),
  by = .(
    scope,
    distance_group,
    cd8_subtype
  )
]

pie_fov_balanced[
  ,
  pie_prop := mean_prop / sum(mean_prop, na.rm = TRUE),
  by = .(
    scope,
    distance_group
  )
]
pie_fov_balanced[, pie_pct := 100 * pie_prop]

# Pooled-cell composition.
# Collapse to one row per FOV-distance first so denominator is not repeated 5 times.
pooled_denoms <- unique(
  analysis_long[, .(
    scope,
    FOV,
    distance_group,
    n_cd8_total
  )]
)

pooled_num <- analysis_long[
  ,
  .(
    pooled_subtype_count = sum(n_subtype, na.rm = TRUE)
  ),
  by = .(
    scope,
    distance_group,
    cd8_subtype
  )
]

pooled_denom <- pooled_denoms[
  ,
  .(
    pooled_total_CD8 = sum(n_cd8_total, na.rm = TRUE)
  ),
  by = .(
    scope,
    distance_group
  )
]

pie_pooled <- merge(
  pooled_num,
  pooled_denom,
  by = c("scope", "distance_group"),
  all.x = TRUE
)

pie_pooled[, pooled_prop := pooled_subtype_count / pooled_total_CD8]
pie_pooled[
  ,
  pie_prop := pooled_prop / sum(pooled_prop, na.rm = TRUE),
  by = .(
    scope,
    distance_group
  )
]
pie_pooled[, pie_pct := 100 * pie_prop]

save_csv(pie_fov_balanced, "03_4_4_pie_composition_FOV_balanced")
save_csv(pie_pooled, "03_4_4_pie_composition_pooled_cells")

# ============================================================
# 14 pie plot function
# ============================================================
make_pie_page <- function(scope_i) {
  if (pie_summary_mode == "FOV_balanced") {
    d <- copy(pie_fov_balanced[scope == scope_i])
    subtitle_txt <- paste0(
      "FOV-balanced mean CD8 subtype composition; ",
      "FOV-distance groups require >= ",
      min_cd8_cells_per_FOV_group,
      " total CD8 cells"
    )
  } else {
    d <- copy(pie_pooled[scope == scope_i])
    subtitle_txt <- "Pooled CD8-cell composition across FOVs"
  }
  
  if (nrow(d) == 0) {
    return(
      ggplot() +
        theme_void() +
        labs(title = paste0(scope_labels[[scope_i]], ": no valid CD8 data"))
    )
  }
  
  d[, distance_plot := factor(
    as.character(distance_group),
    levels = c("Within", "Beyond"),
    labels = c(within_label, beyond_label)
  )]
  
  d[, subtype_label := factor(
    as.character(cd8_subtype),
    levels = cd8_subtypes,
    labels = unname(cd8_subtype_labels[cd8_subtypes])
  )]
  
  # Match the proven 03_4_2 pie-label logic exactly:
  # do NOT manually calculate cumulative y positions because ggplot's
  # stacking order follows the fill-factor stack and can differ from
  # the current data.table row order. position_stack(vjust = 0.5)
  # guarantees that every percentage is centered in its own wedge.
  d[, pie_label := ifelse(
    pie_pct >= label_min_pct,
    sprintf("%.1f%%", pie_pct),
    ""
  )]
  
  ggplot(
    d,
    aes(x = "", y = pie_prop, fill = cd8_subtype)
  ) +
    geom_col(width = 1, color = "white", linewidth = 0.7) +
    coord_polar(theta = "y") +
    geom_text(
      aes(label = pie_label),
      position = position_stack(vjust = 0.5),
      size = 3.4,
      fontface = "bold",
      color = "white"
    ) +
    facet_wrap(~ distance_plot, nrow = 1) +
    scale_fill_manual(
      values = cd8_subtype_cols,
      breaks = cd8_subtypes,
      labels = cd8_subtype_labels[cd8_subtypes],
      drop = FALSE
    ) +
    labs(
      title = paste0("CD8 subtype composition — ", scope_labels[[scope_i]]),
      subtitle = subtitle_txt,
      x = NULL,
      y = NULL,
      fill = "CD8 subtype"
    ) +
    theme_void(base_size = 12) +
    theme(
      strip.background = element_rect(fill = "grey93", color = NA),
      strip.text = element_text(face = "bold", size = 11),
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 9.5),
      legend.position = "right"
    )
}

# ============================================================
# 15 paired violin statistics
# ============================================================
run_violin_stats <- function(d, value_col, scope_i, denominator_name) {
  pair_wide <- dcast(
    d[
      is.finite(get(value_col)),
      .(
        FOV,
        cd8_subtype,
        distance_group,
        value = get(value_col)
      )
    ],
    FOV + cd8_subtype ~ distance_group,
    value.var = "value"
  )
  
  if (!all(c("Within", "Beyond") %in% names(pair_wide))) {
    return(list(stats = data.table(), paired_FOVs = character()))
  }
  
  pair_wide <- pair_wide[is.finite(Within) & is.finite(Beyond)]
  
  if (nrow(pair_wide) == 0) {
    return(list(stats = data.table(), paired_FOVs = character()))
  }
  
  stats_dt <- pair_wide[
    ,
    .(
      n_paired_FOV = .N,
      mean_Within = mean(Within, na.rm = TRUE),
      mean_Beyond = mean(Beyond, na.rm = TRUE),
      median_Within = as.numeric(median(as.numeric(Within), na.rm = TRUE)),
      median_Beyond = as.numeric(median(as.numeric(Beyond), na.rm = TRUE)),
      mean_delta_Within_minus_Beyond = mean(Within - Beyond, na.rm = TRUE),
      median_delta_Within_minus_Beyond = as.numeric(median(as.numeric(Within - Beyond), na.rm = TRUE)),
      p_raw = run_paired_test(Within, Beyond, test = paired_test)
    ),
    by = cd8_subtype
  ]
  
  stats_dt[, p_BH5 := p.adjust(p_raw, method = "BH")]
  stats_dt[, star := vapply(p_raw, p_to_star, character(1))]
  stats_dt[, scope := scope_i]
  stats_dt[, denominator := denominator_name]
  stats_dt[, paired_test := paired_test]
  
  list(
    stats = stats_dt,
    paired_FOVs = unique(pair_wide$FOV)
  )
}

# ============================================================
# 16 violin plot function
# ============================================================
make_violin_page <- function(scope_i, denominator = c("CD8", "nontumor")) {
  denominator <- match.arg(denominator)
  d <- copy(analysis_long[scope == scope_i])
  
  if (denominator == "CD8") {
    value_col <- "pct_within_CD8"
    valid_col <- "valid_cd8_denom"
    y_lab <- "% of CD8 cells"
    subtitle_txt <- paste0(
      "Each point = one FOV/patient; subtype / total CD8; ",
      "FOV-distance groups require >= ",
      min_cd8_cells_per_FOV_group,
      " total CD8 cells"
    )
  } else {
    value_col <- "pct_within_nontumor"
    valid_col <- "valid_nontumor_denom"
    y_lab <- "% of all non-tumor cells"
    subtitle_txt <- paste0(
      "Each point = one FOV/patient; subtype / all non-tumor cells; ",
      "FOV-distance groups require >= ",
      min_nontumor_cells_per_FOV_group,
      " non-tumor cells"
    )
  }
  
  d <- d[get(valid_col) == TRUE & is.finite(get(value_col))]
  
  stats_obj <- run_violin_stats(
    d = d,
    value_col = value_col,
    scope_i = scope_i,
    denominator_name = denominator
  )
  
  stats_dt <- stats_obj$stats
  
  if (length(stats_obj$paired_FOVs) > 0) {
    d <- d[FOV %in% stats_obj$paired_FOVs]
  } else {
    d <- d[0]
  }
  
  if (nrow(d) == 0) {
    return(list(
      plot = ggplot() + theme_void() +
        labs(title = paste0(scope_labels[[scope_i]], ": no paired FOVs available")),
      stats = stats_dt
    ))
  }
  
  d[, x_num := fifelse(as.character(distance_group) == "Within", 1, 2)]
  
  d[, subtype_label := factor(
    as.character(cd8_subtype),
    levels = cd8_subtypes,
    labels = unname(cd8_subtype_labels[cd8_subtypes])
  )]
  
  heights <- d[
    ,
    .(
      ymin = min(get(value_col), na.rm = TRUE),
      ymax = max(get(value_col), na.rm = TRUE)
    ),
    by = cd8_subtype
  ]
  
  heights[, yrange := pmax(ymax - ymin, 1e-6)]
  heights[, y_bracket := ymax + 0.10 * yrange]
  heights[, y_text := ymax + 0.13 * yrange]
  heights[, y_tip := 0.03 * yrange]
  
  ann <- merge(stats_dt, heights, by = "cd8_subtype", all.x = TRUE)
  
  ann[, subtype_label := factor(
    as.character(cd8_subtype),
    levels = cd8_subtypes,
    labels = unname(cd8_subtype_labels[cd8_subtypes])
  )]
  
  ann[, p_label := paste0(
    ifelse(paired_test == "wilcoxon", "paired Wilcoxon p=", "paired t p="),
    vapply(p_raw, format_p, character(1)),
    " ",
    star
  )]
  
  p <- ggplot(
    d,
    aes(
      x = x_num,
      y = get(value_col),
      fill = distance_group,
      group = distance_group
    )
  ) +
    geom_violin(
      scale = "width",
      trim = FALSE,
      color = "black",
      linewidth = 0.3,
      alpha = 0.88
    ) +
    geom_boxplot(
      width = 0.18,
      outlier.shape = NA,
      alpha = 0.30,
      color = "black",
      linewidth = 0.35
    )
  
  if (show_paired_lines) {
    p <- p + geom_line(
      aes(group = FOV),
      linewidth = 0.25,
      alpha = 0.20,
      color = "grey45"
    )
  }
  
  p <- p +
    geom_jitter(
      width = 0.10,
      height = 0,
      size = 1.2,
      alpha = 0.58,
      color = "black"
    ) +
    stat_summary(
      fun = median,
      geom = "point",
      shape = 21,
      size = 1.7,
      fill = "white",
      color = "black",
      stroke = 0.4
    ) +
    geom_segment(
      data = ann,
      aes(x = 1, xend = 2, y = y_bracket, yend = y_bracket),
      inherit.aes = FALSE,
      linewidth = 0.40
    ) +
    geom_segment(
      data = ann,
      aes(x = 1, xend = 1, y = y_bracket - y_tip, yend = y_bracket),
      inherit.aes = FALSE,
      linewidth = 0.40
    ) +
    geom_segment(
      data = ann,
      aes(x = 2, xend = 2, y = y_bracket - y_tip, yend = y_bracket),
      inherit.aes = FALSE,
      linewidth = 0.40
    ) +
    geom_text(
      data = ann,
      aes(x = 1.5, y = y_text, label = p_label),
      inherit.aes = FALSE,
      size = 3.0,
      fontface = "bold",
      vjust = 0
    ) +
    scale_fill_manual(
      values = distance_cols,
      breaks = c("Within", "Beyond"),
      labels = c(within_label, beyond_label),
      drop = FALSE,
      guide = "none"
    ) +
    scale_x_continuous(
      breaks = c(1, 2),
      labels = c(within_label, beyond_label),
      limits = c(0.55, 2.45)
    ) +
    facet_wrap(
      ~ subtype_label,
      scales = "free_y",
      ncol = violin_ncol
    ) +
    labs(
      title = paste0(
        "CD8 subtype proportions by tumor distance — ",
        scope_labels[[scope_i]]
      ),
      subtitle = subtitle_txt,
      x = NULL,
      y = y_lab
    ) +
    theme_bw(base_size = 12) +
    theme(
      strip.background = element_rect(fill = "grey92", color = NA),
      strip.text = element_text(face = "bold", size = 10),
      axis.text.x = element_text(size = 10, face = "bold"),
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 9.5),
      plot.margin = margin(8, 16, 8, 8),
      panel.spacing.y = grid::unit(0.9, "lines")
    ) +
    coord_cartesian(clip = "off")
  
  list(plot = p, stats = stats_dt)
}

# ============================================================
# 17 generate pie PDF: 4 pages
# ============================================================
pie_pdf <- file.path(
  fig_dir,
  paste0(
    "03_4_4_CD8_subtype_pies_cutoff_",
    safe_tag(as.character(distance_cutoff_um)),
    "um_",
    pie_summary_mode,
    ".pdf"
  )
)

if (file.exists(pie_pdf)) file.remove(pie_pdf)

grDevices::cairo_pdf(
  filename = pie_pdf,
  width = pie_pdf_width,
  height = pie_pdf_height
)

for (scope_i in scope_order) {
  print(make_pie_page(scope_i))
}

grDevices::dev.off()
message("Saved: ", pie_pdf)

# ============================================================
# 18 generate violin PDF 1: subtype / all CD8
# ============================================================
violin_cd8_pdf <- file.path(
  fig_dir,
  paste0(
    "03_4_4_violin_subtype_over_CD8_cutoff_",
    safe_tag(as.character(distance_cutoff_um)),
    "um.pdf"
  )
)

if (file.exists(violin_cd8_pdf)) file.remove(violin_cd8_pdf)

grDevices::cairo_pdf(
  filename = violin_cd8_pdf,
  width = violin_pdf_width,
  height = violin_pdf_height
)

violin_stats_CD8 <- list()

for (scope_i in scope_order) {
  obj <- make_violin_page(scope_i = scope_i, denominator = "CD8")
  print(obj$plot)
  if (nrow(obj$stats) > 0) violin_stats_CD8[[scope_i]] <- obj$stats
}

grDevices::dev.off()
message("Saved: ", violin_cd8_pdf)

# ============================================================
# 19 generate violin PDF 2: subtype / all non-tumor cells
# ============================================================
violin_nontumor_pdf <- file.path(
  fig_dir,
  paste0(
    "03_4_4_violin_subtype_over_nontumor_cutoff_",
    safe_tag(as.character(distance_cutoff_um)),
    "um.pdf"
  )
)

if (file.exists(violin_nontumor_pdf)) file.remove(violin_nontumor_pdf)

grDevices::cairo_pdf(
  filename = violin_nontumor_pdf,
  width = violin_pdf_width,
  height = violin_pdf_height
)

violin_stats_nontumor <- list()

for (scope_i in scope_order) {
  obj <- make_violin_page(scope_i = scope_i, denominator = "nontumor")
  print(obj$plot)
  if (nrow(obj$stats) > 0) violin_stats_nontumor[[scope_i]] <- obj$stats
}

grDevices::dev.off()
message("Saved: ", violin_nontumor_pdf)

# ============================================================
# 20 save violin statistics
# ============================================================
stats_CD8_dt <- if (length(violin_stats_CD8) > 0) {
  rbindlist(violin_stats_CD8, use.names = TRUE, fill = TRUE)
} else {
  data.table()
}

stats_nontumor_dt <- if (length(violin_stats_nontumor) > 0) {
  rbindlist(violin_stats_nontumor, use.names = TRUE, fill = TRUE)
} else {
  data.table()
}

stats_all <- rbindlist(
  list(stats_CD8_dt, stats_nontumor_dt),
  use.names = TRUE,
  fill = TRUE
)

if (nrow(stats_all) > 0) {
  stats_all[, subtype_label := cd8_subtype_labels[as.character(cd8_subtype)]]
  
  wanted_order <- c(
    "scope",
    "denominator",
    "cd8_subtype",
    "subtype_label",
    "paired_test",
    "n_paired_FOV",
    "mean_Within",
    "mean_Beyond",
    "median_Within",
    "median_Beyond",
    "mean_delta_Within_minus_Beyond",
    "median_delta_Within_minus_Beyond",
    "p_raw",
    "p_BH5",
    "star"
  )
  
  setcolorder(stats_all, intersect(wanted_order, names(stats_all)))
  save_csv(stats_all, "03_4_4_violin_paired_statistics")
}

# ============================================================
# 21 paired-FOV availability QC
# ============================================================
paired_qc_list <- list()

for (scope_i in scope_order) {
  d1 <- unique(
    analysis_long[
      scope == scope_i & valid_cd8_denom == TRUE,
      .(FOV, distance_group)
    ]
  )
  d1[, present := 1L]
  
  wide1 <- dcast(
    d1,
    FOV ~ distance_group,
    value.var = "present",
    fun.aggregate = sum,
    fill = 0
  )
  
  n_pair_cd8 <- if (all(c("Within", "Beyond") %in% names(wide1))) {
    nrow(wide1[Within > 0 & Beyond > 0])
  } else {
    0L
  }
  
  d2 <- unique(
    analysis_long[
      scope == scope_i & valid_nontumor_denom == TRUE,
      .(FOV, distance_group)
    ]
  )
  d2[, present := 1L]
  
  wide2 <- dcast(
    d2,
    FOV ~ distance_group,
    value.var = "present",
    fun.aggregate = sum,
    fill = 0
  )
  
  n_pair_nt <- if (all(c("Within", "Beyond") %in% names(wide2))) {
    nrow(wide2[Within > 0 & Beyond > 0])
  } else {
    0L
  }
  
  paired_qc_list[[scope_i]] <- data.table(
    scope = scope_i,
    paired_FOV_CD8_denominator = n_pair_cd8,
    paired_FOV_nontumor_denominator = n_pair_nt
  )
}

paired_qc <- rbindlist(paired_qc_list)
save_csv(paired_qc, "03_4_4_paired_FOV_QC")

# ============================================================
# 22 parameters
# ============================================================
parameter_dt <- data.table(
  parameter = c(
    "metadata_file",
    "distance_cutoff_um",
    "distance_group_Within",
    "distance_group_Beyond",
    "min_cd8_cells_per_FOV_group",
    "min_nontumor_cells_per_FOV_group",
    "pie_summary_mode",
    "paired_test",
    "biological_replicate",
    "CD8_broad_definition",
    "CD8_subtypes",
    "nontumor_definition",
    "immunotype_order",
    "immunotype_colors",
    "pie_pdf_width",
    "pie_pdf_height",
    "violin_pdf_width",
    "violin_pdf_height"
  ),
  value = c(
    meta_file,
    as.character(distance_cutoff_um),
    paste0("distance <= ", distance_cutoff_um, " um"),
    paste0("distance > ", distance_cutoff_um, " um"),
    as.character(min_cd8_cells_per_FOV_group),
    as.character(min_nontumor_cells_per_FOV_group),
    pie_summary_mode,
    paired_test,
    "one patient = one FOV",
    'BC6k_refined == "T.cell.CD8"',
    paste(cd8_subtypes, collapse = ";"),
    'BC6k_refined != "tumor"',
    paste(immunotype_order, collapse = " -> "),
    "ID=#000000; IE=#0000FF; IN=#FF0000",
    as.character(pie_pdf_width),
    as.character(pie_pdf_height),
    as.character(violin_pdf_width),
    as.character(violin_pdf_height)
  )
)

save_csv(parameter_dt, "03_4_4_parameters")

# ============================================================
# 23 finished
# ============================================================
message("")
message("============================================================")
message("03_4_4 finished.")
message("============================================================")
message("Distance cutoff: ", distance_cutoff_um, " um")
message("Group 1: distance <= ", distance_cutoff_um, " um")
message("Group 2: distance >  ", distance_cutoff_um, " um")
message("")
message("Pie summary mode: ", pie_summary_mode)
message("Violin paired test: ", paired_test)
message("")
message("Outputs:")
message("  1) 4-page CD8 subtype pie PDF: All FOVs / ID / IE / IN")
message("  2) 4-page violin PDF, denominator = all CD8")
message("  3) 4-page violin PDF, denominator = all non-tumor cells")
message("  4) FOV-level counts/proportions CSV + RDS")
message("  5) FOV-balanced and pooled pie-summary CSVs")
message("  6) paired-test statistics CSV")
message("  7) QC and parameter CSVs")
message("")
message("Figures: ", fig_dir)
message("Data:    ", data_dir)
message("============================================================")

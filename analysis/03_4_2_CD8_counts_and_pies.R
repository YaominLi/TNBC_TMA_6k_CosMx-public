rm(list = ls())
gc()
graphics.off()
options(stringsAsFactors = FALSE)

# =========================================================
# 03_4_2
# CD8 subtype summary by immunotype and selected distance windows
#
# Outputs:
# 1) Overall total CD8 cell counts across IN / IE / ID
# 2) Overall CD8 subtype composition pie charts across IN / IE / ID
# 3) For each user-defined distance window:
#    - total CD8 cell counts across IN / IE / ID
#    - CD8 subtype composition pie charts across IN / IE / ID
#
# Composition logic:
#   For each immunotype, first calculate subtype proportion within each FOV
#   among all CD8 cells in the relevant scope, then average across FOVs
#   (FOV-balanced composition).
#
# Count logic:
#   Raw pooled CD8 cell counts.
#
# Current metadata:
#   results/03_2/data/03_2_meta_with_immunotype_BW_top_40.rds
#   disttonearesttumor is already treated as um.
# =========================================================

# =========================================================
# 01 paths
# =========================================================
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

proj <- tryCatch(normalizePath(getwd(), winslash = "/"), error = function(e) NA)
if (is.na(proj) || !file.exists(file.path(proj, "results", "03_2", "data"))) {
  proj <- find_project_root()
}

in_dir  <- file.path(proj, "results", "03_2", "data")
out_dir <- file.path(proj, "results", "03_4_2")
fig_dir <- file.path(out_dir, "figures")
dat_dir <- file.path(out_dir, "data")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dat_dir, recursive = TRUE, showWarnings = FALSE)

# =========================================================
# 02 packages
# =========================================================
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# =========================================================
# 03 user settings
# =========================================================

# Arbitrary continuous distance windows in um.
# Examples:
# distance_windows <- list(c(30, 50))
# distance_windows <- list(c(45, 55))
# distance_windows <- list(c(0, 20), c(20, 40), c(40, 60))
#
# Intervals are treated as [low, high), preventing overlap
# when adjacent windows are used together.
distance_windows <- list(
  c(0, 30)
)

fine_col <- "celltype_new"

cd8_levels <- c(
  "cd8_cytotoxic",
  "cd8_exhausted",
  "cd8_tem",
  "cd8_tcm",
  "cd8_naive"
)

# Immunotype display order: ID -> IE -> IN.
# Colors: ID black, IE blue, IN red.
immunotype_cols <- c(
  "IN" = "#FF0000",
  "IE" = "#0000FF",
  "ID" = "#000000"
)

# CD8 subtype colors:
# Cytotoxic red, naive blue, intermediate states in between.
cd8_cols <- c(
  "cd8_cytotoxic" = "#D55E00",  # vermillion
  "cd8_exhausted" = "#CC79A7",  # reddish purple
  "cd8_tem"       = "#E69F00",  # orange
  "cd8_tcm"       = "#56B4E9",  # sky blue
  "cd8_naive"     = "#0072B2"   # deep blue
)

# Minimum total CD8 cells required in one FOV for inclusion in
# FOV-balanced composition summaries.
min_cd8_cells_per_fov_all <- 3L

# Minimum total CD8 cells required in one FOV-window for inclusion in
# distance-window FOV-balanced composition summaries.
min_cd8_cells_per_fov_window <- 3L

# Hide very small pie-sector labels.
label_min_pct <- 4

# =========================================================
# 04 helpers
# =========================================================
save_dat <- function(x, filename, ext = c("csv", "rds")) {
  ext <- match.arg(ext)
  out <- file.path(dat_dir, paste0(filename, ".", ext))
  if (file.exists(out)) file.remove(out)
  
  if (ext == "csv") {
    write.csv(as.data.frame(x), out, row.names = FALSE)
  } else {
    saveRDS(x, out)
  }
  
  message("Saved: ", out)
}

save_pdf <- function(p, filename, w = 8.0, h = 4.8) {
  out <- file.path(fig_dir, paste0(filename, ".pdf"))
  if (file.exists(out)) file.remove(out)
  
  ggplot2::ggsave(
    out,
    plot = p,
    width = w,
    height = h,
    device = cairo_pdf
  )
  
  message("Saved: ", out)
}

safe_tag <- function(x) gsub("[^A-Za-z0-9_]+", "_", x)

pretty_cd8_label <- function(x) {
  out <- gsub("^cd8_", "", x)
  out <- gsub("_", " ", out)
  tools::toTitleCase(out)
}

# =========================================================
# 05 load data
# =========================================================
meta_file <- file.path(in_dir, "03_2_meta_with_immunotype_BW_top_40.rds")
if (!file.exists(meta_file)) {
  stop("Missing file: ", meta_file)
}

dt <- as.data.table(readRDS(meta_file))

need_cols <- c(
  "cell_id",
  "FOV",
  "Run_Tissue_name",
  "BC6k_refined",
  fine_col,
  "immunotype_BW_top_40",
  "disttonearesttumor"
)

miss_cols <- setdiff(need_cols, names(dt))
if (length(miss_cols) > 0) {
  stop("Missing required columns: ", paste(miss_cols, collapse = ", "))
}

dt <- dt[
  !is.na(FOV) &
    !is.na(Run_Tissue_name) &
    !is.na(get(fine_col)) &
    !is.na(immunotype_BW_top_40) &
    !is.na(disttonearesttumor)
]

dt[, FOV := as.character(FOV)]
dt[, Run_Tissue_name := as.character(Run_Tissue_name)]
dt[, immunotype_BW_top_40 := as.character(immunotype_BW_top_40)]
dt[, (fine_col) := as.character(get(fine_col))]
dt[, dist_to_tumor_um := as.numeric(disttonearesttumor)]

dt <- dt[
  immunotype_BW_top_40 %in% c("ID", "IE", "IN") &
    is.finite(dist_to_tumor_um) &
    dist_to_tumor_um >= 0
]

valid_cd8 <- intersect(cd8_levels, sort(unique(dt[[fine_col]])))
if (length(valid_cd8) == 0) {
  stop("No requested CD8 subtype labels were found in ", fine_col)
}

missing_cd8 <- setdiff(cd8_levels, valid_cd8)
if (length(missing_cd8) > 0) {
  warning("Missing CD8 subtypes and skipped: ", paste(missing_cd8, collapse = ", "))
}

dt_cd8 <- copy(dt[get(fine_col) %in% valid_cd8])
dt_cd8[, CD8_subtype := factor(get(fine_col), levels = cd8_levels)]
dt_cd8[, immunotype_BW_top_40 := factor(immunotype_BW_top_40, levels = c("ID", "IE", "IN"))]

if (nrow(dt_cd8) == 0) {
  stop("No CD8 subtype cells remain after filtering.")
}

# =========================================================
# 06 general summary functions
# =========================================================

# Raw pooled counts
make_count_summary <- function(x, scope_label = "all", low_um = NA_real_, high_um = NA_real_) {
  total_dt <- x[, .(n_cd8 = .N), by = .(immunotype_BW_top_40)]
  fill_dt <- data.table(immunotype_BW_top_40 = factor(c("ID", "IE", "IN"), levels = c("ID", "IE", "IN")))
  total_dt <- merge(fill_dt, total_dt, by = "immunotype_BW_top_40", all.x = TRUE)
  total_dt[is.na(n_cd8), n_cd8 := 0L]
  
  subtype_dt <- x[, .(n_subtype = .N), by = .(immunotype_BW_top_40, CD8_subtype)]
  fill_sub <- CJ(
    immunotype_BW_top_40 = factor(c("ID", "IE", "IN"), levels = c("ID", "IE", "IN")),
    CD8_subtype = factor(cd8_levels, levels = cd8_levels),
    unique = TRUE
  )
  subtype_dt <- merge(fill_sub, subtype_dt, by = c("immunotype_BW_top_40", "CD8_subtype"), all.x = TRUE)
  subtype_dt[is.na(n_subtype), n_subtype := 0L]
  subtype_dt <- merge(subtype_dt, total_dt, by = "immunotype_BW_top_40", all.x = TRUE)
  subtype_dt[, pooled_pct := fifelse(n_cd8 > 0, n_subtype / n_cd8, 0)]
  subtype_dt[, `:=`(
    scope = scope_label,
    distance_low_um = low_um,
    distance_high_um = high_um
  )]
  
  total_dt[, `:=`(
    scope = scope_label,
    distance_low_um = low_um,
    distance_high_um = high_um
  )]
  
  list(total = total_dt, subtype = subtype_dt)
}

# FOV-balanced composition
make_fov_balanced_composition <- function(
    x,
    scope_label = "all",
    low_um = NA_real_,
    high_um = NA_real_,
    min_cells = 3L
) {
  fov_total <- x[, .(n_cd8_fov = .N), by = .(FOV, Run_Tissue_name, immunotype_BW_top_40)]
  fov_total <- fov_total[n_cd8_fov >= min_cells]
  
  if (nrow(fov_total) == 0) {
    return(list(fov = data.table(), summary = data.table()))
  }
  
  x2 <- merge(
    x,
    fov_total,
    by = c("FOV", "Run_Tissue_name", "immunotype_BW_top_40"),
    all = FALSE
  )
  
  subtype_count <- x2[, .(n_subtype = .N), by = .(FOV, Run_Tissue_name, immunotype_BW_top_40, CD8_subtype)]
  
  base <- unique(fov_total[, .(FOV, Run_Tissue_name, immunotype_BW_top_40, n_cd8_fov)])
  subtype_grid <- data.table(CD8_subtype = factor(valid_cd8, levels = cd8_levels))
  
  base[, join_key := 1L]
  subtype_grid[, join_key := 1L]
  full_grid <- merge(base, subtype_grid, by = "join_key", allow.cartesian = TRUE)[, join_key := NULL]
  
  fov_prop <- merge(
    full_grid,
    subtype_count,
    by = c("FOV", "Run_Tissue_name", "immunotype_BW_top_40", "CD8_subtype"),
    all.x = TRUE
  )
  
  fov_prop[is.na(n_subtype), n_subtype := 0L]
  fov_prop[, proportion := n_subtype / n_cd8_fov]
  fov_prop[, `:=`(
    scope = scope_label,
    distance_low_um = low_um,
    distance_high_um = high_um
  )]
  
  summary_dt <- fov_prop[
    ,
    .(
      n_FOV = uniqueN(FOV),
      mean_prop = mean(proportion, na.rm = TRUE),
      sd_prop = sd(proportion, na.rm = TRUE),
      sem_prop = sd(proportion, na.rm = TRUE) / sqrt(uniqueN(FOV)),
      total_subtype_cells = sum(n_subtype, na.rm = TRUE)
    ),
    by = .(scope, distance_low_um, distance_high_um, immunotype_BW_top_40, CD8_subtype)
  ]
  
  summary_dt[is.na(sd_prop), sd_prop := 0]
  summary_dt[is.na(sem_prop), sem_prop := 0]
  
  list(fov = fov_prop, summary = summary_dt)
}

# =========================================================
# 07 overall summaries (all distances)
# =========================================================
overall_count <- make_count_summary(
  dt_cd8,
  scope_label = "all_distance",
  low_um = NA_real_,
  high_um = NA_real_
)

overall_comp <- make_fov_balanced_composition(
  dt_cd8,
  scope_label = "all_distance",
  low_um = NA_real_,
  high_um = NA_real_,
  min_cells = min_cd8_cells_per_fov_all
)

overall_total_dt <- overall_count$total
overall_subtype_pooled_dt <- overall_count$subtype
overall_fov_comp_dt <- overall_comp$fov
overall_pie_dt <- overall_comp$summary

save_dat(overall_total_dt, "03_4_2_overall_CD8_total_counts", "csv")
save_dat(overall_total_dt, "03_4_2_overall_CD8_total_counts", "rds")
save_dat(overall_subtype_pooled_dt, "03_4_2_overall_CD8_subtype_pooled_counts", "csv")
save_dat(overall_subtype_pooled_dt, "03_4_2_overall_CD8_subtype_pooled_counts", "rds")
save_dat(overall_fov_comp_dt, "03_4_2_overall_CD8_subtype_FOVlevel", "csv")
save_dat(overall_fov_comp_dt, "03_4_2_overall_CD8_subtype_FOVlevel", "rds")
save_dat(overall_pie_dt, "03_4_2_overall_CD8_subtype_FOVbalanced_summary", "csv")
save_dat(overall_pie_dt, "03_4_2_overall_CD8_subtype_FOVbalanced_summary", "rds")

# =========================================================
# 08 window summaries
# =========================================================
all_window_total <- list()
all_window_subtype_pooled <- list()
all_window_fov <- list()
all_window_pie <- list()

for (i in seq_along(distance_windows)) {
  rng <- distance_windows[[i]]
  if (length(rng) != 2 || any(!is.finite(rng))) {
    stop("Each distance window must contain exactly two finite numbers.")
  }
  
  low_um  <- min(rng)
  high_um <- max(rng)
  
  if (low_um < 0 || high_um <= low_um) {
    stop("Invalid distance window: ", low_um, "-", high_um)
  }
  
  scope_label <- paste0(low_um, "-", high_um, "um")
  message("Processing distance window: ", scope_label)
  
  sub <- dt_cd8[dist_to_tumor_um >= low_um & dist_to_tumor_um < high_um]
  
  if (nrow(sub) == 0) {
    warning("No CD8 cells found in window ", scope_label, ". Skipped.")
    next
  }
  
  tmp_count <- make_count_summary(
    sub,
    scope_label = scope_label,
    low_um = low_um,
    high_um = high_um
  )
  
  tmp_comp <- make_fov_balanced_composition(
    sub,
    scope_label = scope_label,
    low_um = low_um,
    high_um = high_um,
    min_cells = min_cd8_cells_per_fov_window
  )
  
  all_window_total[[length(all_window_total) + 1]] <- tmp_count$total
  all_window_subtype_pooled[[length(all_window_subtype_pooled) + 1]] <- tmp_count$subtype
  
  if (nrow(tmp_comp$fov) > 0) {
    all_window_fov[[length(all_window_fov) + 1]] <- tmp_comp$fov
  }
  
  if (nrow(tmp_comp$summary) > 0) {
    all_window_pie[[length(all_window_pie) + 1]] <- tmp_comp$summary
  }
}

window_total_dt <- if (length(all_window_total) > 0) rbindlist(all_window_total, fill = TRUE) else data.table()
window_subtype_pooled_dt <- if (length(all_window_subtype_pooled) > 0) rbindlist(all_window_subtype_pooled, fill = TRUE) else data.table()
window_fov_comp_dt <- if (length(all_window_fov) > 0) rbindlist(all_window_fov, fill = TRUE) else data.table()
window_pie_dt <- if (length(all_window_pie) > 0) rbindlist(all_window_pie, fill = TRUE) else data.table()

if (nrow(window_total_dt) > 0) {
  save_dat(window_total_dt, "03_4_2_window_CD8_total_counts", "csv")
  save_dat(window_total_dt, "03_4_2_window_CD8_total_counts", "rds")
}
if (nrow(window_subtype_pooled_dt) > 0) {
  save_dat(window_subtype_pooled_dt, "03_4_2_window_CD8_subtype_pooled_counts", "csv")
  save_dat(window_subtype_pooled_dt, "03_4_2_window_CD8_subtype_pooled_counts", "rds")
}
if (nrow(window_fov_comp_dt) > 0) {
  save_dat(window_fov_comp_dt, "03_4_2_window_CD8_subtype_FOVlevel", "csv")
  save_dat(window_fov_comp_dt, "03_4_2_window_CD8_subtype_FOVlevel", "rds")
}
if (nrow(window_pie_dt) > 0) {
  save_dat(window_pie_dt, "03_4_2_window_CD8_subtype_FOVbalanced_summary", "csv")
  save_dat(window_pie_dt, "03_4_2_window_CD8_subtype_FOVbalanced_summary", "rds")
}

# =========================================================
# 09 plotting helpers
# =========================================================
plot_count_bar <- function(plot_dt, title_text, subtitle_text) {
  plot_dt <- copy(plot_dt)
  plot_dt[, immunotype_BW_top_40 := factor(immunotype_BW_top_40, levels = c("ID", "IE", "IN"))]
  
  ggplot(plot_dt, aes(x = immunotype_BW_top_40, y = n_cd8, fill = immunotype_BW_top_40)) +
    geom_col(width = 0.72, color = "black", linewidth = 0.35) +
    geom_text(aes(label = n_cd8), vjust = -0.35, size = 4.2) +
    scale_fill_manual(values = immunotype_cols, drop = FALSE) +
    labs(
      x = NULL,
      y = "Number of CD8 cells",
      title = title_text,
      subtitle = subtitle_text
    ) +
    expand_limits(y = max(plot_dt$n_cd8, na.rm = TRUE) * 1.12) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(size = 15, face = "bold"),
      plot.subtitle = element_text(size = 10.5),
      axis.title.y = element_text(size = 12),
      axis.text.x = element_text(size = 11, face = "bold"),
      axis.text.y = element_text(size = 10.5),
      legend.position = "none"
    )
}

plot_pie_facets <- function(plot_dt, title_text, subtitle_text) {
  plot_dt <- copy(plot_dt)
  plot_dt[, pct := mean_prop * 100]
  plot_dt[, label_text := ifelse(pct >= label_min_pct, sprintf("%.1f%%", pct), "")]
  plot_dt[, subtype_label := factor(pretty_cd8_label(as.character(CD8_subtype)),
                                    levels = pretty_cd8_label(cd8_levels))]
  
  facet_info <- unique(plot_dt[, .(immunotype_BW_top_40, n_FOV)])
  facet_info[, facet_label := paste0(as.character(immunotype_BW_top_40), "\n", "n=", n_FOV, " FOVs")]
  plot_dt <- merge(plot_dt, facet_info, by = c("immunotype_BW_top_40", "n_FOV"), all.x = TRUE)
  facet_order <- facet_info[order(factor(immunotype_BW_top_40, levels = c("ID", "IE", "IN"))), facet_label]
  plot_dt[, facet_label := factor(facet_label, levels = facet_order)]
  
  ggplot(plot_dt, aes(x = "", y = mean_prop, fill = CD8_subtype)) +
    geom_col(width = 1, color = "white", linewidth = 0.7) +
    coord_polar(theta = "y") +
    geom_text(
      aes(label = label_text),
      position = position_stack(vjust = 0.5),
      size = 3.4,
      fontface = "bold",
      color = "white"
    ) +
    facet_wrap(~ facet_label, nrow = 1) +
    scale_fill_manual(
      values = cd8_cols,
      breaks = cd8_levels,
      drop = FALSE,
      labels = pretty_cd8_label(cd8_levels)
    ) +
    labs(
      title = title_text,
      subtitle = subtitle_text,
      fill = "CD8 subtype"
    ) +
    theme_void(base_size = 11) +
    theme(
      plot.title = element_text(size = 15, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 10, hjust = 0),
      strip.text = element_text(size = 11, face = "bold"),
      legend.title = element_text(size = 10.5, face = "bold"),
      legend.text = element_text(size = 9.5),
      legend.position = "right"
    )
}

# =========================================================
# 10 overall plots
# =========================================================
p_overall_bar <- plot_count_bar(
  overall_total_dt,
  title_text = "Total CD8 cell counts across immunotypes",
  subtitle_text = "Pooled raw cell counts"
)

save_pdf(
  p_overall_bar,
  "03_4_2_overall_CD8_total_counts_bar",
  w = 6.0,
  h = 4.5
)

if (nrow(overall_pie_dt) > 0) {
  p_overall_pie <- plot_pie_facets(
    overall_pie_dt,
    title_text = "Overall CD8 subtype composition across immunotypes",
    subtitle_text = "Mean FOV-level composition among all CD8 cells"
  )
  
  save_pdf(
    p_overall_pie,
    "03_4_2_overall_CD8_subtype_composition_pie",
    w = 9.2,
    h = 4.8
  )
}

# =========================================================
# 11 distance-window plots
# =========================================================
if (nrow(window_total_dt) > 0) {
  for (scope_i in unique(window_total_dt$scope)) {
    bar_dt <- window_total_dt[scope == scope_i]
    
    p_bar <- plot_count_bar(
      bar_dt,
      title_text = paste0("CD8 cell counts at ", scope_i),
      subtitle_text = "Pooled raw cell counts within the selected distance window"
    )
    
    save_pdf(
      p_bar,
      paste0("03_4_2_window_CD8_total_counts_bar_", safe_tag(scope_i)),
      w = 6.0,
      h = 4.5
    )
  }
}

if (nrow(window_pie_dt) > 0) {
  for (scope_i in unique(window_pie_dt$scope)) {
    pie_dt <- window_pie_dt[scope == scope_i]
    
    p_pie <- plot_pie_facets(
      pie_dt,
      title_text = paste0("CD8 subtype composition at ", scope_i),
      subtitle_text = "Mean FOV-level composition among CD8 cells within the selected distance window"
    )
    
    save_pdf(
      p_pie,
      paste0("03_4_2_window_CD8_subtype_composition_pie_", safe_tag(scope_i)),
      w = 9.2,
      h = 4.8
    )
  }
}

# =========================================================
# 12 combined summary tables
# =========================================================
combined_total_dt <- rbindlist(
  list(overall_total_dt, window_total_dt),
  fill = TRUE
)

combined_subtype_pooled_dt <- rbindlist(
  list(overall_subtype_pooled_dt, window_subtype_pooled_dt),
  fill = TRUE
)

combined_pie_dt <- rbindlist(
  list(overall_pie_dt, window_pie_dt),
  fill = TRUE
)

save_dat(combined_total_dt, "03_4_2_combined_CD8_total_counts", "csv")
save_dat(combined_total_dt, "03_4_2_combined_CD8_total_counts", "rds")

save_dat(combined_subtype_pooled_dt, "03_4_2_combined_CD8_subtype_pooled_counts", "csv")
save_dat(combined_subtype_pooled_dt, "03_4_2_combined_CD8_subtype_pooled_counts", "rds")

save_dat(combined_pie_dt, "03_4_2_combined_CD8_subtype_FOVbalanced_summary", "csv")
save_dat(combined_pie_dt, "03_4_2_combined_CD8_subtype_FOVbalanced_summary", "rds")

message("03_4_2 finished successfully.")

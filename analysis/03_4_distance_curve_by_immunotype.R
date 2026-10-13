rm(list = ls())
gc()
graphics.off()

options(stringsAsFactors = FALSE)

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
out_dir <- file.path(proj, "results", "03_4")
fig_dir <- file.path(out_dir, "figures")
dat_dir <- file.path(out_dir, "data")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dat_dir, recursive = TRUE, showWarnings = FALSE)

save_dat_overwrite <- function(x, filename, ext = c("csv", "rds")) {
  ext <- match.arg(ext)
  out <- file.path(dat_dir, paste0(filename, ".", ext))
  if (file.exists(out)) file.remove(out)
  if (ext == "csv") write.csv(x, out, row.names = FALSE)
  if (ext == "rds") saveRDS(x, out)
  message("Saved: ", out)
  invisible(out)
}

save_plot_pdf <- function(plot_obj, filename, w = 7, h = 6) {
  out <- file.path(fig_dir, paste0(filename, ".pdf"))
  if (file.exists(out)) file.remove(out)
  ggplot2::ggsave(out, plot = plot_obj, width = w, height = h, device = cairo_pdf)
  message("Saved: ", out)
  invisible(out)
}

# =========================================================
# 02 packages
# =========================================================
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# =========================================================
# 03 parameters
# =========================================================
bin_width <- 10
x_max_use <- 50
min_cells_per_fov_bin <- 5

# colors consistent with 03_2
immunotype_cols <- c(
  "IN" = "#FF0000",
  "IE" = "#0000FF",
  "ID" = "#000000"
)

# set NULL to plot all non-tumor cell types
cell_types_to_plot <- NULL

# fitted curve
use_loess_curve <- TRUE
loess_span <- 0.75

# =========================================================
# 04 load data
# =========================================================
meta_file <- file.path(in_dir, "03_2_meta_with_immunotype_BW_top_40.rds")
if (!file.exists(meta_file)) {
  stop("Missing file: ", meta_file)
}

dt <- as.data.table(readRDS(meta_file))

# =========================================================
# 05 checks
# =========================================================
need_cols <- c(
  "cell_id", "FOV", "Run_Tissue_name",
  "BC6k_refined", "immunotype_BW_top_40",
  "disttonearesttumor"
)
miss <- setdiff(need_cols, names(dt))
if (length(miss) > 0) {
  stop("Missing required columns: ", paste(miss, collapse = ", "))
}

dt <- dt[
  !is.na(FOV) &
    !is.na(Run_Tissue_name) &
    !is.na(BC6k_refined) &
    !is.na(immunotype_BW_top_40) &
    !is.na(disttonearesttumor)
]

dt[, FOV := as.character(FOV)]
dt[, Run_Tissue_name := as.character(Run_Tissue_name)]
dt[, BC6k_refined := as.character(BC6k_refined)]
dt[, immunotype_BW_top_40 := as.character(immunotype_BW_top_40)]
dt[, disttonearesttumor := as.numeric(disttonearesttumor)]

dt <- dt[immunotype_BW_top_40 %in% c("ID", "IE", "IN")]

# =========================================================
# 06 follow 03_2 logic: exclude tumor cells from the start
# Finalized metadata distances are already in micrometers.
dt <- dt[BC6k_refined != "tumor"]

dt[, dist_to_tumor_um := as.numeric(disttonearesttumor)]

dt <- dt[dist_to_tumor_um >= 0]
dt <- dt[dist_to_tumor_um <= x_max_use]

# =========================================================
# 07 define distance bins
# =========================================================
breaks_use <- seq(0, x_max_use, by = bin_width)
if (tail(breaks_use, 1) < x_max_use) {
  breaks_use <- c(breaks_use, x_max_use)
}

dt[, dist_bin := cut(
  dist_to_tumor_um,
  breaks = breaks_use,
  include.lowest = TRUE,
  right = FALSE
)]

dt <- dt[!is.na(dist_bin)]

bin_map <- data.table(
  dist_bin = levels(dt$dist_bin),
  bin_left = head(breaks_use, -1),
  bin_right = tail(breaks_use, -1)
)
bin_map[, dist_mid := (bin_left + bin_right) / 2]

# =========================================================
# 08 choose cell types to plot
# =========================================================
all_cell_types <- sort(unique(dt$BC6k_refined))

if (is.null(cell_types_to_plot)) {
  cell_types_to_plot <- all_cell_types
} else {
  cell_types_to_plot <- intersect(cell_types_to_plot, all_cell_types)
}

if (length(cell_types_to_plot) == 0) {
  stop("No valid non-tumor cell types found in cell_types_to_plot.")
}

dt <- dt[BC6k_refined %in% cell_types_to_plot]

# =========================================================
# 09 FOV-bin level composition
# denominator = all non-tumor cells in that FOV-bin
# numerator   = each BC6k_refined cell type count in that FOV-bin
# =========================================================
fov_bin_total <- dt[, .(
  n_bin_total = .N
), by = .(FOV, Run_Tissue_name, immunotype_BW_top_40, dist_bin)]

fov_bin_total <- fov_bin_total[n_bin_total >= min_cells_per_fov_bin]

dt_use <- merge(
  dt,
  fov_bin_total,
  by = c("FOV", "Run_Tissue_name", "immunotype_BW_top_40", "dist_bin"),
  all = FALSE
)

if (nrow(dt_use) == 0) {
  stop("No data left after applying min_cells_per_fov_bin filter.")
}

fov_bin_celltype <- dt_use[, .(
  n_celltype = .N
), by = .(FOV, Run_Tissue_name, immunotype_BW_top_40, dist_bin, BC6k_refined)]

fov_bin_pct <- merge(
  fov_bin_celltype,
  fov_bin_total,
  by = c("FOV", "Run_Tissue_name", "immunotype_BW_top_40", "dist_bin"),
  all.x = TRUE
)

fov_bin_pct[, proportion := n_celltype / n_bin_total]

# make missing cell types explicit as 0 within valid FOV-bin
fov_bin_base <- unique(fov_bin_total[, .(FOV, Run_Tissue_name, immunotype_BW_top_40, dist_bin)])
celltype_dt <- data.table(BC6k_refined = cell_types_to_plot)

fov_bin_base[, join_dummy := 1L]
celltype_dt[, join_dummy := 1L]

all_grid <- merge(
  fov_bin_base,
  celltype_dt,
  by = "join_dummy",
  allow.cartesian = TRUE
)[, join_dummy := NULL]

fov_bin_pct_full <- merge(
  all_grid,
  fov_bin_pct,
  by = c("FOV", "Run_Tissue_name", "immunotype_BW_top_40", "dist_bin", "BC6k_refined"),
  all.x = TRUE
)

fov_bin_pct_full <- merge(
  fov_bin_pct_full,
  fov_bin_total,
  by = c("FOV", "Run_Tissue_name", "immunotype_BW_top_40", "dist_bin"),
  all.x = TRUE
)

fov_bin_pct_full[is.na(n_celltype), n_celltype := 0L]
fov_bin_pct_full[is.na(proportion), proportion := 0]

fov_bin_pct_full <- merge(fov_bin_pct_full, bin_map, by = "dist_bin", all.x = TRUE)

# =========================================================
# 10 summarize across FOVs: mean ± SEM
# =========================================================
summary_dt <- fov_bin_pct_full[, .(
  n_FOV = uniqueN(FOV),
  mean_prop = mean(proportion, na.rm = TRUE),
  sd_prop = sd(proportion, na.rm = TRUE),
  sem_prop = sd(proportion, na.rm = TRUE) / sqrt(uniqueN(FOV))
), by = .(immunotype_BW_top_40, BC6k_refined, dist_bin, dist_mid, bin_left, bin_right)]

summary_dt[is.na(sd_prop), sd_prop := 0]
summary_dt[is.na(sem_prop), sem_prop := 0]

setorder(summary_dt, BC6k_refined, immunotype_BW_top_40, dist_mid)

# =========================================================
# 11 diagnostics
# =========================================================
diag_bins <- summary_dt[, .(
  n_dist_bins = uniqueN(dist_mid),
  min_dist = min(dist_mid, na.rm = TRUE),
  max_dist = max(dist_mid, na.rm = TRUE)
), by = .(BC6k_refined, immunotype_BW_top_40)]

save_dat_overwrite(as.data.frame(diag_bins), "03_4_distance_curve_diagnostic_bins", ext = "csv")

# =========================================================
# 12 save data
# =========================================================
save_dat_overwrite(as.data.frame(fov_bin_pct_full), "03_4_FOV_bin_celltype_proportion_long", ext = "csv")
save_dat_overwrite(as.data.frame(fov_bin_pct_full), "03_4_FOV_bin_celltype_proportion_long", ext = "rds")

save_dat_overwrite(as.data.frame(summary_dt), "03_4_distance_curve_summary", ext = "csv")
save_dat_overwrite(as.data.frame(summary_dt), "03_4_distance_curve_summary", ext = "rds")

# =========================================================
# 13 plotting helper
# =========================================================
make_curve_plot <- function(plot_dt, celltype_name = NULL) {
  
  p <- ggplot(
    plot_dt,
    aes(
      x = dist_mid,
      y = mean_prop,
      color = immunotype_BW_top_40,
      fill = immunotype_BW_top_40,
      group = immunotype_BW_top_40
    )
  ) +
    geom_ribbon(
      aes(ymin = pmax(mean_prop - sem_prop, 0), ymax = mean_prop + sem_prop),
      alpha = 0.18,
      color = NA
    ) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.2) +
    scale_color_manual(values = immunotype_cols, breaks = c("IN", "IE", "ID")) +
    scale_fill_manual(values = immunotype_cols, breaks = c("IN", "IE", "ID")) +
    scale_x_continuous(limits = c(0, x_max_use), expand = c(0, 0)) +
    scale_y_continuous(labels = function(x) paste0(round(x * 100), "%")) +
    labs(
      x = "Distance to nearest tumor (um)",
      y = "Cell relative frequency (%)",
      title = if (!is.null(celltype_name)) celltype_name else NULL
    ) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", size = 10),
      legend.title = element_blank(),
      legend.text = element_text(size = 8),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 8)
    )
  
  loess_ok <- plot_dt[, .(n_x = uniqueN(dist_mid)), by = immunotype_BW_top_40]
  if (use_loess_curve && any(loess_ok$n_x >= 4)) {
    smooth_dt <- plot_dt[immunotype_BW_top_40 %in% loess_ok[n_x >= 4, immunotype_BW_top_40]]
    p <- p + geom_smooth(
      data = smooth_dt,
      aes(group = immunotype_BW_top_40),
      method = "loess",
      formula = y ~ x,
      span = loess_span,
      se = FALSE,
      linewidth = 0.8,
      linetype = 2,
      alpha = 0.8
    )
  }
  
  p
}

# =========================================================
# 14 combined facet plot
# =========================================================
facet_dt <- copy(summary_dt)
facet_dt[, BC6k_refined := factor(BC6k_refined, levels = sort(unique(BC6k_refined)))]

p_facet <- ggplot(
  facet_dt,
  aes(
    x = dist_mid,
    y = mean_prop,
    color = immunotype_BW_top_40,
    fill = immunotype_BW_top_40,
    group = immunotype_BW_top_40
  )
) +
  geom_ribbon(
    aes(ymin = pmax(mean_prop - sem_prop, 0), ymax = mean_prop + sem_prop),
    alpha = 0.15,
    color = NA
  ) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 0.9) +
  scale_color_manual(values = immunotype_cols, breaks = c("IN", "IE", "ID")) +
  scale_fill_manual(values = immunotype_cols, breaks = c("IN", "IE", "ID")) +
  scale_x_continuous(limits = c(0, x_max_use), expand = c(0, 0)) +
  scale_y_continuous(labels = function(x) paste0(round(x * 100), "%")) +
  labs(
    x = "Distance to nearest tumor (um)",
    y = "Cell relative frequency (%)",
    title = "Cell-type distribution across continuous tumor distance",
    subtitle = "03_2 logic: tumor cells excluded; FOV-level bin proportions summarized as mean ± SEM"
  ) +
  facet_wrap(~ BC6k_refined, scales = "free_y") +
  theme_bw(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(size = 8, face = "bold"),
    plot.title = element_text(face = "bold", size = 10),
    plot.subtitle = element_text(size = 8),
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8)
  )

loess_ok_facet <- facet_dt[, .(n_x = uniqueN(dist_mid)), by = .(BC6k_refined, immunotype_BW_top_40)]
valid_loess_facet <- loess_ok_facet[n_x >= 4, .(BC6k_refined, immunotype_BW_top_40)]

if (use_loess_curve && nrow(valid_loess_facet) > 0) {
  smooth_facet_dt <- merge(
    facet_dt,
    valid_loess_facet,
    by = c("BC6k_refined", "immunotype_BW_top_40"),
    all = FALSE
  )
  
  p_facet <- p_facet + geom_smooth(
    data = smooth_facet_dt,
    aes(group = immunotype_BW_top_40),
    method = "loess",
    formula = y ~ x,
    span = loess_span,
    se = FALSE,
    linewidth = 0.7,
    linetype = 2,
    alpha = 0.8
  )
}

save_plot_pdf(
  p_facet,
  "03_4_distance_curve_facet_all_celltypes",
  w = 12,
  h = max(8, 2.2 * ceiling(length(unique(facet_dt$BC6k_refined)) / 3))
)

# =========================================================
# 15 one pdf per cell type
# =========================================================
for (ct in unique(summary_dt$BC6k_refined)) {
  plot_dt <- summary_dt[BC6k_refined == ct]
  p_ct <- make_curve_plot(plot_dt, celltype_name = ct)
  save_plot_pdf(
    p_ct,
    paste0("03_4_distance_curve_", gsub("[^A-Za-z0-9_]+", "_", ct)),
    w = 6.2,
    h = 4.6
  )
}

# =========================================================
# 16 extra plots: smooth-curve-only version
# keep only fitted loess curves
# =========================================================

make_smooth_only_plot <- function(plot_dt, celltype_name = NULL) {
  
  p <- ggplot(
    plot_dt,
    aes(
      x = dist_mid,
      y = mean_prop,
      color = immunotype_BW_top_40,
      group = immunotype_BW_top_40
    )
  ) +
    scale_color_manual(values = immunotype_cols, breaks = c("IN", "IE", "ID")) +
    scale_x_continuous(limits = c(0, x_max_use), expand = c(0, 0)) +
    scale_y_continuous(labels = function(x) paste0(round(x * 100), "%")) +
    labs(
      x = "Distance to nearest tumor (um)",
      y = "Cell relative frequency (%)",
      title = if (!is.null(celltype_name)) celltype_name else NULL
    ) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", size = 10),
      legend.title = element_blank(),
      legend.text = element_text(size = 8),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 8)
    )
  
  loess_ok <- plot_dt[, .(n_x = uniqueN(dist_mid)), by = immunotype_BW_top_40]
  if (any(loess_ok$n_x >= 4)) {
    smooth_dt <- plot_dt[immunotype_BW_top_40 %in% loess_ok[n_x >= 4, immunotype_BW_top_40]]
    p <- p + geom_smooth(
      data = smooth_dt,
      aes(group = immunotype_BW_top_40),
      method = "loess",
      formula = y ~ x,
      span = loess_span,
      se = FALSE,
      linewidth = 1.0,
      alpha = 0.95
    )
  }
  
  p
}

# ---------- combined facet plot: smooth only ----------
facet_dt_smooth <- copy(summary_dt)
facet_dt_smooth[, BC6k_refined := factor(BC6k_refined, levels = sort(unique(BC6k_refined)))]

p_facet_smooth <- ggplot(
  facet_dt_smooth,
  aes(
    x = dist_mid,
    y = mean_prop,
    color = immunotype_BW_top_40,
    group = immunotype_BW_top_40
  )
) +
  scale_color_manual(values = immunotype_cols, breaks = c("IN", "IE", "ID")) +
  scale_x_continuous(limits = c(0, x_max_use), expand = c(0, 0)) +
  scale_y_continuous(labels = function(x) paste0(round(x * 100), "%")) +
  labs(
    x = "Distance to nearest tumor (um)",
    y = "Cell relative frequency (%)",
    title = "Cell-type distribution across continuous tumor distance",
    subtitle = "Smooth fitted curves only"
  ) +
  facet_wrap(~ BC6k_refined, scales = "free_y") +
  theme_bw(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(size = 8, face = "bold"),
    plot.title = element_text(face = "bold", size = 10),
    plot.subtitle = element_text(size = 8),
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8)
  )

loess_ok_facet_smooth <- facet_dt_smooth[, .(n_x = uniqueN(dist_mid)), by = .(BC6k_refined, immunotype_BW_top_40)]
valid_loess_facet_smooth <- loess_ok_facet_smooth[n_x >= 4, .(BC6k_refined, immunotype_BW_top_40)]

if (nrow(valid_loess_facet_smooth) > 0) {
  smooth_facet_dt <- merge(
    facet_dt_smooth,
    valid_loess_facet_smooth,
    by = c("BC6k_refined", "immunotype_BW_top_40"),
    all = FALSE
  )
  
  p_facet_smooth <- p_facet_smooth + geom_smooth(
    data = smooth_facet_dt,
    aes(group = immunotype_BW_top_40),
    method = "loess",
    formula = y ~ x,
    span = loess_span,
    se = FALSE,
    linewidth = 0.95,
    alpha = 0.95
  )
}

save_plot_pdf(
  p_facet_smooth,
  "03_4_distance_curve_facet_all_celltypes_smooth_only",
  w = 12,
  h = max(8, 2.2 * ceiling(length(unique(facet_dt_smooth$BC6k_refined)) / 3))
)

# ---------- one pdf per cell type: smooth only ----------
for (ct in unique(summary_dt$BC6k_refined)) {
  plot_dt <- summary_dt[BC6k_refined == ct]
  p_ct_smooth <- make_smooth_only_plot(plot_dt, celltype_name = ct)
  save_plot_pdf(
    p_ct_smooth,
    paste0("03_4_distance_curve_", gsub("[^A-Za-z0-9_]+", "_", ct), "_smooth_only"),
    w = 6.2,
    h = 4.6
  )
}


message("03_4 smooth-only plots finished successfully.")

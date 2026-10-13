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
out_dir <- file.path(proj, "results", "03_4_1")
fig_dir <- file.path(out_dir, "figures")
dat_dir <- file.path(out_dir, "data")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dat_dir, recursive = TRUE, showWarnings = FALSE)

save_dat <- function(x, nm, ext = c("csv", "rds")) {
  ext <- match.arg(ext)
  f <- file.path(dat_dir, paste0(nm, ".", ext))
  if (file.exists(f)) file.remove(f)
  if (ext == "csv") write.csv(x, f, row.names = FALSE)
  if (ext == "rds") saveRDS(x, f)
  message("Saved: ", f)
}

save_pdf <- function(p, nm, w = 6.2, h = 4.6) {
  f <- file.path(fig_dir, paste0(nm, ".pdf"))
  if (file.exists(f)) file.remove(f)
  ggplot2::ggsave(f, p, width = w, height = h, device = cairo_pdf)
  message("Saved: ", f)
}

safe_tag <- function(x) gsub("[^A-Za-z0-9_]+", "_", x)

# =========================================================
# 02 packages / parameters
# =========================================================
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

bin_width <- 10
x_max_use <- 50
min_non_tumor_cells_per_fov_bin <- 5
min_cd8_cells_per_fov_bin <- 3
loess_span <- 0.75

fine_col <- "celltype_new"
cd8_subtypes <- c(
  "cd8_cytotoxic",
  "cd8_exhausted",
  "cd8_naive",
  "cd8_tcm",
  "cd8_tem"
)

immunotype_cols <- c(
  "IN" = "#FF0000",
  "IE" = "#0000FF",
  "ID" = "#000000"
)

# =========================================================
# 03 load / check
# =========================================================
meta_file <- file.path(in_dir, "03_2_meta_with_immunotype_BW_top_40.rds")
if (!file.exists(meta_file)) stop("Missing file: ", meta_file)

dt <- as.data.table(readRDS(meta_file))

need <- c(
  "cell_id", "FOV", "Run_Tissue_name",
  "BC6k_refined", fine_col,
  "immunotype_BW_top_40", "disttonearesttumor"
)
miss <- setdiff(need, names(dt))
if (length(miss) > 0) stop("Missing required columns: ", paste(miss, collapse = ", "))

dt <- dt[
  !is.na(FOV) &
    !is.na(Run_Tissue_name) &
    !is.na(BC6k_refined) &
    !is.na(get(fine_col)) &
    !is.na(immunotype_BW_top_40) &
    !is.na(disttonearesttumor)
]

dt[, FOV := as.character(FOV)]
dt[, Run_Tissue_name := as.character(Run_Tissue_name)]
dt[, BC6k_refined := as.character(BC6k_refined)]
dt[, (fine_col) := as.character(get(fine_col))]
dt[, immunotype_BW_top_40 := as.character(immunotype_BW_top_40)]
dt[, dist_to_tumor_um := as.numeric(disttonearesttumor)]

dt <- dt[
  immunotype_BW_top_40 %in% c("ID", "IE", "IN") &
    BC6k_refined != "tumor" &
    dist_to_tumor_um >= 0 &
    dist_to_tumor_um <= x_max_use
]

if (nrow(dt) == 0) stop("No non-tumor cells remain after filtering.")

# =========================================================
# 04 distance bins
# =========================================================
breaks_use <- seq(0, x_max_use, by = bin_width)
if (tail(breaks_use, 1) < x_max_use) breaks_use <- c(breaks_use, x_max_use)

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

valid_cd8 <- intersect(cd8_subtypes, sort(unique(dt[[fine_col]])))
if (length(valid_cd8) == 0) stop("No requested CD8 subtype found in ", fine_col)

save_dat(
  data.frame(
    CD8_subtype = cd8_subtypes,
    available = cd8_subtypes %in% valid_cd8
  ),
  "03_4_1_CD8_subtype_manifest",
  "csv"
)

# =========================================================
# 05 helper: complete missing subtypes as zero
# =========================================================
complete_grid <- function(valid_bins, count_dt, denom_col) {
  base <- unique(valid_bins[, .(
    FOV, Run_Tissue_name, immunotype_BW_top_40, dist_bin
  )])
  subtype_dt <- data.table(CD8_subtype = valid_cd8)
  
  base[, key := 1L]
  subtype_dt[, key := 1L]
  grid <- merge(base, subtype_dt, by = "key", allow.cartesian = TRUE)[, key := NULL]
  
  out <- merge(
    grid, count_dt,
    by = c("FOV", "Run_Tissue_name", "immunotype_BW_top_40", "dist_bin", "CD8_subtype"),
    all.x = TRUE
  )
  
  denom <- unique(valid_bins[, c(
    "FOV", "Run_Tissue_name", "immunotype_BW_top_40", "dist_bin", denom_col
  ), with = FALSE])
  
  out <- merge(
    out, denom,
    by = c("FOV", "Run_Tissue_name", "immunotype_BW_top_40", "dist_bin"),
    all.x = TRUE
  )
  
  out[is.na(n_subtype), n_subtype := 0L]
  out
}

# =========================================================
# 06 denominator A: all non-tumor cells
# =========================================================
non_tumor_total <- dt[, .(
  n_non_tumor_total = .N
), by = .(FOV, Run_Tissue_name, immunotype_BW_top_40, dist_bin)]

non_tumor_total <- non_tumor_total[
  n_non_tumor_total >= min_non_tumor_cells_per_fov_bin
]

dt_valid_non_tumor <- merge(
  dt, non_tumor_total,
  by = c("FOV", "Run_Tissue_name", "immunotype_BW_top_40", "dist_bin"),
  all = FALSE
)

cd8_counts_A <- dt_valid_non_tumor[
  get(fine_col) %in% valid_cd8,
  .(n_subtype = .N),
  by = .(
    FOV, Run_Tissue_name, immunotype_BW_top_40, dist_bin,
    CD8_subtype = get(fine_col)
  )
]

fov_A <- complete_grid(
  non_tumor_total,
  cd8_counts_A,
  "n_non_tumor_total"
)
fov_A[, proportion := n_subtype / n_non_tumor_total]
fov_A[, denominator_mode := "all_non_tumor"]
fov_A <- merge(fov_A, bin_map, by = "dist_bin", all.x = TRUE)

# =========================================================
# 07 denominator B: all CD8 cells
# =========================================================
dt_cd8 <- dt[get(fine_col) %in% valid_cd8]

cd8_total <- dt_cd8[, .(
  n_cd8_total = .N
), by = .(FOV, Run_Tissue_name, immunotype_BW_top_40, dist_bin)]

cd8_total <- cd8_total[
  n_cd8_total >= min_cd8_cells_per_fov_bin
]

dt_valid_cd8 <- merge(
  dt_cd8, cd8_total,
  by = c("FOV", "Run_Tissue_name", "immunotype_BW_top_40", "dist_bin"),
  all = FALSE
)

cd8_counts_B <- dt_valid_cd8[, .(
  n_subtype = .N
), by = .(
  FOV, Run_Tissue_name, immunotype_BW_top_40, dist_bin,
  CD8_subtype = get(fine_col)
)]

fov_B <- complete_grid(
  cd8_total,
  cd8_counts_B,
  "n_cd8_total"
)
fov_B[, proportion := n_subtype / n_cd8_total]
fov_B[, denominator_mode := "all_CD8"]
fov_B <- merge(fov_B, bin_map, by = "dist_bin", all.x = TRUE)

# =========================================================
# 08 summarize across FOVs
# =========================================================
summarize_one <- function(x) {
  y <- x[, .(
    n_FOV = uniqueN(FOV),
    mean_prop = mean(proportion, na.rm = TRUE),
    sd_prop = sd(proportion, na.rm = TRUE),
    sem_prop = sd(proportion, na.rm = TRUE) / sqrt(uniqueN(FOV))
  ), by = .(
    denominator_mode,
    immunotype_BW_top_40,
    CD8_subtype,
    dist_bin, dist_mid, bin_left, bin_right
  )]
  
  y[is.na(sd_prop), sd_prop := 0]
  y[is.na(sem_prop), sem_prop := 0]
  y
}

summary_A <- summarize_one(fov_A)
summary_B <- summarize_one(fov_B)
summary_dt <- rbindlist(list(summary_A, summary_B), fill = TRUE)

setorder(
  summary_dt,
  denominator_mode,
  CD8_subtype,
  immunotype_BW_top_40,
  dist_mid
)

save_dat(fov_A, "03_4_1_FOV_bin_CD8_subtype_vs_all_non_tumor", "csv")
save_dat(fov_A, "03_4_1_FOV_bin_CD8_subtype_vs_all_non_tumor", "rds")
save_dat(fov_B, "03_4_1_FOV_bin_CD8_subtype_within_CD8", "csv")
save_dat(fov_B, "03_4_1_FOV_bin_CD8_subtype_within_CD8", "rds")
save_dat(summary_dt, "03_4_1_CD8_distance_curve_summary", "csv")
save_dat(summary_dt, "03_4_1_CD8_distance_curve_summary", "rds")

diag_dt <- summary_dt[, .(
  n_dist_bins = uniqueN(dist_mid),
  min_n_FOV = min(n_FOV, na.rm = TRUE),
  max_n_FOV = max(n_FOV, na.rm = TRUE)
), by = .(
  denominator_mode,
  CD8_subtype,
  immunotype_BW_top_40
)]
save_dat(diag_dt, "03_4_1_CD8_distance_curve_diagnostic", "csv")

# =========================================================
# 09 plot helpers
# =========================================================
y_lab_for <- function(mode) {
  if (mode == "all_CD8") {
    "Fraction within all CD8 cells (%)"
  } else {
    "CD8 subtype frequency among non-tumor cells (%)"
  }
}

subtitle_for <- function(mode) {
  if (mode == "all_CD8") {
    "Denominator: all CD8 cells in each FOV-distance bin"
  } else {
    "Denominator: all non-tumor cells in each FOV-distance bin"
  }
}

add_loess <- function(p, plot_dt, smooth_only = FALSE) {
  ok <- plot_dt[, .(n_x = uniqueN(dist_mid)), by = immunotype_BW_top_40]
  valid_groups <- ok[n_x >= 4, immunotype_BW_top_40]
  
  if (length(valid_groups) > 0) {
    sm <- plot_dt[immunotype_BW_top_40 %in% valid_groups]
    p <- p + geom_smooth(
      data = sm,
      aes(group = immunotype_BW_top_40),
      method = "loess",
      formula = y ~ x,
      span = loess_span,
      se = FALSE,
      linewidth = if (smooth_only) 1.0 else 0.8,
      linetype = if (smooth_only) 1 else 2,
      alpha = 0.9
    )
  }
  p
}

make_one_plot <- function(plot_dt, subtype, mode, smooth_only = FALSE) {
  p <- ggplot(
    plot_dt,
    aes(
      x = dist_mid,
      y = mean_prop,
      color = immunotype_BW_top_40,
      fill = immunotype_BW_top_40,
      group = immunotype_BW_top_40
    )
  )
  
  if (!smooth_only) {
    p <- p +
      geom_ribbon(
        aes(
          ymin = pmax(mean_prop - sem_prop, 0),
          ymax = mean_prop + sem_prop
        ),
        alpha = 0.18,
        color = NA
      ) +
      geom_line(linewidth = 0.9) +
      geom_point(size = 1.2)
  }
  
  p <- p +
    scale_color_manual(values = immunotype_cols, breaks = c("IN", "IE", "ID")) +
    scale_fill_manual(values = immunotype_cols, breaks = c("IN", "IE", "ID")) +
    scale_x_continuous(limits = c(0, x_max_use), expand = c(0, 0)) +
    scale_y_continuous(labels = function(x) paste0(round(x * 100), "%")) +
    labs(
      x = "Distance to nearest tumor (um)",
      y = y_lab_for(mode),
      title = subtype,
      subtitle = if (smooth_only) {
        paste0(subtitle_for(mode), "; smooth fitted curves only")
      } else {
        subtitle_for(mode)
      }
    ) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", size = 10),
      plot.subtitle = element_text(size = 8),
      legend.title = element_blank(),
      legend.text = element_text(size = 8),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 8)
    )
  
  add_loess(p, plot_dt, smooth_only)
}

make_facet <- function(plot_dt, mode, smooth_only = FALSE) {
  tmp <- copy(plot_dt)
  tmp[, CD8_subtype := factor(CD8_subtype, levels = valid_cd8)]
  
  p <- ggplot(
    tmp,
    aes(
      x = dist_mid,
      y = mean_prop,
      color = immunotype_BW_top_40,
      fill = immunotype_BW_top_40,
      group = immunotype_BW_top_40
    )
  )
  
  if (!smooth_only) {
    p <- p +
      geom_ribbon(
        aes(
          ymin = pmax(mean_prop - sem_prop, 0),
          ymax = mean_prop + sem_prop
        ),
        alpha = 0.15,
        color = NA
      ) +
      geom_line(linewidth = 0.8) +
      geom_point(size = 0.9)
  }
  
  ok <- tmp[, .(n_x = uniqueN(dist_mid)), by = .(
    CD8_subtype, immunotype_BW_top_40
  )]
  valid <- ok[n_x >= 4, .(CD8_subtype, immunotype_BW_top_40)]
  
  if (nrow(valid) > 0) {
    sm <- merge(
      tmp, valid,
      by = c("CD8_subtype", "immunotype_BW_top_40"),
      all = FALSE
    )
    
    p <- p + geom_smooth(
      data = sm,
      aes(group = immunotype_BW_top_40),
      method = "loess",
      formula = y ~ x,
      span = loess_span,
      se = FALSE,
      linewidth = if (smooth_only) 0.95 else 0.7,
      linetype = if (smooth_only) 1 else 2,
      alpha = 0.9
    )
  }
  
  p +
    scale_color_manual(values = immunotype_cols, breaks = c("IN", "IE", "ID")) +
    scale_fill_manual(values = immunotype_cols, breaks = c("IN", "IE", "ID")) +
    scale_x_continuous(limits = c(0, x_max_use), expand = c(0, 0)) +
    scale_y_continuous(labels = function(x) paste0(round(x * 100), "%")) +
    labs(
      x = "Distance to nearest tumor (um)",
      y = y_lab_for(mode),
      title = "CD8 subtype distribution across continuous tumor distance",
      subtitle = if (smooth_only) {
        paste0(subtitle_for(mode), "; smooth fitted curves only")
      } else {
        paste0(subtitle_for(mode), "; mean +/- SEM across FOVs")
      }
    ) +
    facet_wrap(~ CD8_subtype, scales = "free_y") +
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
}

# =========================================================
# 10 plots
# =========================================================
for (mode in c("all_non_tumor", "all_CD8")) {
  mode_dt <- summary_dt[denominator_mode == mode]
  if (nrow(mode_dt) == 0) next
  
  p1 <- make_facet(mode_dt, mode, smooth_only = FALSE)
  save_pdf(
    p1,
    paste0("03_4_1_CD8_subtypes_facet_", mode),
    11, 7.5
  )
  
  p2 <- make_facet(mode_dt, mode, smooth_only = TRUE)
  save_pdf(
    p2,
    paste0("03_4_1_CD8_subtypes_facet_", mode, "_smooth_only"),
    11, 7.5
  )
  
  for (subtype in valid_cd8) {
    one <- mode_dt[CD8_subtype == subtype]
    if (nrow(one) == 0) next
    
    p3 <- make_one_plot(one, subtype, mode, smooth_only = FALSE)
    save_pdf(
      p3,
      paste0("03_4_1_distance_curve_", safe_tag(subtype), "_", mode)
    )
    
    p4 <- make_one_plot(one, subtype, mode, smooth_only = TRUE)
    save_pdf(
      p4,
      paste0("03_4_1_distance_curve_", safe_tag(subtype), "_", mode, "_smooth_only")
    )
  }
}

# =========================================================
# 11 overview: all CD8 subtypes together within each immunotype
# =========================================================
overview_dt <- copy(summary_B)
overview_dt[, CD8_subtype := factor(CD8_subtype, levels = valid_cd8)]

p_overview <- ggplot(
  overview_dt,
  aes(
    x = dist_mid,
    y = mean_prop,
    color = CD8_subtype,
    group = CD8_subtype
  )
) +
  geom_line(linewidth = 0.85) +
  geom_point(size = 1.0) +
  scale_x_continuous(limits = c(0, x_max_use), expand = c(0, 0)) +
  scale_y_continuous(labels = function(x) paste0(round(x * 100), "%")) +
  labs(
    x = "Distance to nearest tumor (um)",
    y = "Fraction within all CD8 cells (%)",
    color = "CD8 subtype",
    title = "CD8 subtype composition across tumor distance",
    subtitle = "Each panel represents one immunotype"
  ) +
  facet_wrap(~ immunotype_BW_top_40, nrow = 1) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(size = 9, face = "bold"),
    plot.title = element_text(face = "bold", size = 10),
    plot.subtitle = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8)
  )

save_pdf(
  p_overview,
  "03_4_1_CD8_subtype_composition_within_CD8_overview",
  11, 4.6
)

message("03_4_1 CD8 subtype distance analysis finished successfully.")
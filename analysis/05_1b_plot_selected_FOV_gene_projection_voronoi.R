rm(list = ls())
gc()
graphics.off()

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(ggplot2)
  library(deldir)
  library(grid)
  library(cowplot)
})

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
if (is.na(proj) || !file.exists(file.path(proj, "raw"))) {
  proj <- find_project_root()
}

raw_dir <- file.path(proj, "raw")
out_dir <- file.path(proj, "results", "05_1b")
data_dir <- file.path(out_dir, "data")
fig_dir  <- file.path(out_dir, "figures")

dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir,  recursive = TRUE, showWarnings = FALSE)

save_csv_overwrite <- function(x, filename) {
  out <- file.path(data_dir, paste0(filename, ".csv"))
  if (file.exists(out)) file.remove(out)
  fwrite(x, out)
  invisible(out)
}

save_rds_overwrite <- function(x, filename) {
  out <- file.path(data_dir, paste0(filename, ".rds"))
  if (file.exists(out)) file.remove(out)
  saveRDS(x, out)
  invisible(out)
}

# =========================================================
# 02 user input
# =========================================================
target_fovs <- c("s3f1", "s4f79", "s2f2")

gene_set_name <- "Neuroactive ligand-receptor interaction"
score_genes <- c(
  "ADRA1D","GCGR","MTNR1A","PRLHR","NMBR","GRIK1","NTSR2","GRIN2B","GALR2","CRH",
  "C5","TBXA2R","SSTR5","F2","TSHR","SSTR1","CHRNA4","PRL","APELA","POMC",
  "OPRM1","MC4R","GRM8","GALR1","GRIA3","GRIA1","CHRM2","CCK","GRIN2A","CALCR",
  "AVPR1A","SSTR3","PYY","ADRA2B","GRM5","CNR2","GRID1","HTR6","MTNR1B","CHRM5",
  "UCN3","AVP","P2RY1","NPS","DRD1","VIPR2","GLP2R","CYSLTR2","GRIK2","GRID2",
  "GALR3","VIPR1","GH1","HRH3","ADRB3","TACR3","DRD4","ADRB1","HTR1A","MC5R",
  "GRIN1","GALP","MC2R","CORT","PTGER3","CCKBR","NPY2R","GRM3","P2RY14","NPVF",
  "OPRK1","AGTR1","SSTR2","PENK","CTSG","GCG","SSTR4","RXFP4","ADRB2","RLN3",
  "CCKAR","GABRA2","FPR1","AGTR2","QRFPR","PTGIR","CNR1","FSHB","NMU","HTR2C",
  "VIP","PMCH","MAS1","GRM7","AVPR1B","LPAR2","ADRA1A","GRIA4","TACR1","MC3R",
  "P2RX1","P2RX7","OXT","HRH1","GLRA2","CALCB","NPFF","EDNRA","SST","PTGDR",
  "ADM2","SPX","HTR7","GHR","F2RL1","P2RY8","P2RY13","FPR2","ADORA1","GRM1",
  "CHRM3","PNOC","APLN","GABRG3","GABRB2","GLP1R","HTR1B","ADRA2C","TAC1","PTGER2",
  "P2RY11","CHRNA7","NPY5R","PDYN","DRD2","APLNR","GRIA2","LEP","ADORA2A","P2RY10",
  "CALCRL","TRPV1","KNG1","THRA","CRHR2","EDNRB"
)

celltype_col <- "BC6k_refined"
tumor_label <- "tumor"
display_celltype_col <- "BC6k_refined"
cd8_label <- "T.cell.CD8"

# Additional immune-only projection.
# These populations are merged into one "Immune" display class.
# The immune display uses the exact same pink as the existing CD8-only panel.
immune_labels <- c(
  "B.cell",
  "Plasma",
  "Plasmablast",
  "Plasmacytoid.dendritic.cell",
  "Dendritic.cell",
  "Conventional.dendritic.cell",
  "Macrophage",
  "Mast.cell",
  "Monocyte",
  "Neutrophil",
  "NK.cell",
  "T.cell.CD4",
  "T.cell.CD8",
  "T.cell.regulatory"
)

expr_group_method <- "quantile_nonzero_tumor"
high_q <- 0.75

min_tumor_cells <- 10
min_genes_present <- 2

expand_ratio <- 0.01
nn_keep_quantile <- 0.90
nn_keep_mult <- 1.15
polygon_linewidth <- 0.12

plot_bg <- "black"
panel_bg <- "black"
axis_text_col <- "grey80"
axis_line_col <- "grey55"
border_col <- "grey70"

# =========================================================
# 03 palettes
# =========================================================
celltype_cols <- c(
  "B.cell"                      = "#F2A65A",
  "Plasma"                      = "#D98C2B",
  "Plasmablast"                 = "#B87333",
  "Plasmacytoid.dendritic.cell" = "#8C6D31",
  "Conventional.dendritic.cell" = "#8BC34A",
  "Macrophage"                  = "#6DAA2C",
  "Mast.cell"                   = "grey50",
  "Monocyte"                    = "#35B56A",
  "Neutrophil"                  = "#2EB7B0",
  "NK.cell"                     = "#23A6D5",
  "T.cell.CD4"                  = "#4A97C9",
  "T.cell.regulatory"           = "#7B8CC7",
  "T.cell.CD8"                  = "#D864A6",
  "Fibroblast"                  = "#B1AA34",
  "Endothelial"                 = "#C99A2A",
  "Pericyte"                    = "#8E7D22",
  "tumor"                       = "#EE7068"
)

projection_cols <- c(
  "microenvironment" = "grey65",
  "tumor_low"        = "#6FA8DC",
  "tumor_high"       = "#EE7068"
)

cont_low  <- "#FFFFFF"
cont_high <- "#08306B"

cd8_only_cols <- c(
  "other_cells" = "grey65",
  "CD8" = "#D864A6"
)

immune_only_cols <- c(
  "other_cells" = "grey65",
  "Immune" = "#D864A6"
)

# =========================================================
# 04 helper functions
# =========================================================
base_theme <- function() {
  theme_bw(base_size = 9) +
    theme(
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text = element_text(size = 8, color = axis_text_col),
      axis.ticks = element_line(color = axis_line_col, linewidth = 0.2),
      plot.title = element_text(size = 16, face = "plain", color = "white", hjust = 0.5),
      legend.title = element_text(size = 9, color = "white"),
      legend.text = element_text(size = 7.5, color = "white"),
      legend.background = element_rect(fill = plot_bg, color = plot_bg),
      legend.box.background = element_rect(fill = plot_bg, color = plot_bg),
      plot.background = element_rect(fill = plot_bg, color = plot_bg),
      panel.background = element_rect(fill = panel_bg, color = panel_bg),
      panel.border = element_rect(color = border_col, linewidth = 0.4, fill = NA)
    )
}

harmonize_celltype_names <- function(x) {
  x <- as.character(x)
  name_map <- c("Dendritic.cell" = "Conventional.dendritic.cell")
  idx <- x %in% names(name_map)
  x[idx] <- name_map[x[idx]]
  x
}

get_nn1_dist <- function(x, y) {
  mat <- as.matrix(data.frame(x = x, y = y))
  d <- as.matrix(dist(mat))
  diag(d) <- Inf
  apply(d, 1, min, na.rm = TRUE)
}

zscore_safe <- function(x) {
  s <- sd(x, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(0, length(x)))
  (x - m) / s
}

make_voronoi_dt <- function(df, x_col = "x", y_col = "y", id_col = "cell_id", expand_ratio = 0.01) {
  stopifnot(all(c(x_col, y_col, id_col) %in% colnames(df)))
  
  dt0 <- copy(df[, .(
    cell_id = get(id_col),
    x = get(x_col),
    y = get(y_col)
  )])
  
  dt0 <- dt0[!is.na(x) & !is.na(y)]
  
  if (nrow(dt0) < 3) {
    stop("Too few cells in selected FOV for Voronoi tessellation.")
  }
  
  dup_idx <- duplicated(dt0[, .(x, y)]) | duplicated(dt0[, .(x, y)], fromLast = TRUE)
  if (any(dup_idx)) {
    set.seed(1)
    dt0[dup_idx, x := x + runif(.N, -1e-8, 1e-8)]
    dt0[dup_idx, y := y + runif(.N, -1e-8, 1e-8)]
  }
  
  xr <- range(dt0$x, na.rm = TRUE)
  yr <- range(dt0$y, na.rm = TRUE)
  
  xpad <- diff(xr) * expand_ratio
  ypad <- diff(yr) * expand_ratio
  
  rw <- c(xr[1] - xpad, xr[2] + xpad, yr[1] - ypad, yr[2] + ypad)
  
  dd <- deldir::deldir(dt0$x, dt0$y, rw = rw)
  tt <- deldir::tile.list(dd)
  
  poly_list <- vector("list", length(tt))
  for (i in seq_along(tt)) {
    tile_i <- tt[[i]]
    xi <- tile_i$x
    yi <- tile_i$y
    if (length(xi) < 3) next
    poly_list[[i]] <- data.table(
      cell_id = dt0$cell_id[i],
      x = xi,
      y = yi,
      polygon_id = i,
      vertex_order = seq_along(xi)
    )
  }
  
  rbindlist(poly_list, use.names = TRUE, fill = TRUE)
}

build_fov_polygon_dt <- function(dt_fov) {
  dt_fov <- copy(dt_fov)
  dt_fov[, nn1_dist := get_nn1_dist(x, y)]
  nn_keep_cut <- as.numeric(quantile(dt_fov$nn1_dist, probs = nn_keep_quantile, na.rm = TRUE)) * nn_keep_mult
  dt_fov[, keep_polygon := nn1_dist <= nn_keep_cut]
  
  poly_dt <- make_voronoi_dt(
    df = dt_fov[, .(cell_id, x, y)],
    x_col = "x",
    y_col = "y",
    id_col = "cell_id",
    expand_ratio = expand_ratio
  )
  
  annot_dt <- copy(dt_fov[, !c("x", "y"), with = FALSE])
  poly_dt <- merge(poly_dt, annot_dt, by = "cell_id", all.x = TRUE)
  setorder(poly_dt, polygon_id, vertex_order)
  poly_dt_plot <- copy(poly_dt[keep_polygon == TRUE])
  
  list(dt_fov = dt_fov, poly_dt_plot = poly_dt_plot, nn_keep_cut = nn_keep_cut)
}

make_row_label_plot <- function(label_text) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = label_text, angle = 90, color = "white", size = 5) +
    xlim(0, 1) + ylim(0, 1) +
    theme_void() +
    theme(
      plot.background = element_rect(fill = plot_bg, color = plot_bg),
      panel.background = element_rect(fill = plot_bg, color = plot_bg)
    )
}

make_legend_panel <- function(plot_obj) {
  leg <- cowplot::get_legend(plot_obj)
  cowplot::ggdraw() +
    cowplot::draw_grob(leg, x = 0, y = 0, width = 1, height = 1) +
    theme(
      plot.background = element_rect(fill = plot_bg, color = plot_bg),
      panel.background = element_rect(fill = plot_bg, color = plot_bg)
    )
}

make_row_block <- function(plot_list, row_label, legend_plot, legend_width = 0.95) {
  main_plots <- lapply(plot_list, function(p) p + theme(legend.position = "none"))
  left_label <- make_row_label_plot(row_label)
  main_block <- cowplot::plot_grid(plotlist = main_plots, nrow = 1, align = "hv", axis = "tblr")
  right_leg  <- make_legend_panel(legend_plot)
  
  cowplot::plot_grid(
    left_label, main_block, right_leg,
    nrow = 1,
    rel_widths = c(0.35, length(plot_list) * 1.0, legend_width),
    align = "hv",
    axis = "tblr"
  )
}



save_legend_only_pdf <- function(plot_obj, filename, width = 3.2, height = 2.2, bg = "black") {
  leg <- cowplot::get_legend(
    plot_obj + theme(
      legend.position = "right",
      legend.box = "vertical"
    )
  )
  
  p_leg <- cowplot::ggdraw() +
    cowplot::draw_grob(leg, x = 0, y = 0, width = 1, height = 1) +
    theme(
      plot.background = element_rect(fill = bg, color = bg),
      panel.background = element_rect(fill = bg, color = bg)
    )
  
  if (file.exists(filename)) file.remove(filename)
  ggsave(filename, p_leg, width = width, height = height, bg = bg)
}
# =========================================================
# 05 input gene set
# =========================================================
score_name <- as.character(gene_set_name)
score_genes <- unique(as.character(score_genes))
score_genes <- score_genes[!is.na(score_genes) & score_genes != ""]

if (length(score_genes) == 0) stop("No genes provided for scoring.")

save_csv_overwrite(
  data.table(score_name = score_name, input_gene = score_genes),
  paste0("05_1b_", score_name, "_input_gene_set")
)

# =========================================================
# 06 load data
# =========================================================
meta <- as.data.table(readRDS(file.path(raw_dir, "TNBC_6k_metadata_Bogang_selected.RDS")))
xy   <- readRDS(file.path(raw_dir, "TNBC_6k_xy_Bogang.RDS"))
norm <- readRDS(file.path(raw_dir, "TNBC_6k_norm_Bogang.RDS"))

xy_dt <- as.data.table(xy)
xy_dt[, cell_id := rownames(xy)]
setcolorder(xy_dt, c("cell_id", colnames(xy_dt)[1:2]))
setnames(xy_dt, old = c("x_mm", "y_mm"), new = c("x", "y"))

dt <- merge(meta, xy_dt, by = "cell_id", all.x = TRUE)

need_cols <- c("cell_id", "Run_Tissue_name", "FOV", "BC6k_refined", "x", "y")
miss <- setdiff(need_cols, names(dt))
if (length(miss) > 0) stop("Missing required columns: ", paste(miss, collapse = ", "))

dt <- dt[!is.na(FOV) & !is.na(Run_Tissue_name) & !is.na(x) & !is.na(y)]
dt[, FOV := as.character(FOV)]
dt[, Run_Tissue_name := as.character(Run_Tissue_name)]
dt[, BC6k_refined := as.character(BC6k_refined)]

missing_fovs <- setdiff(target_fovs, unique(dt$FOV))
if (length(missing_fovs) > 0) stop("These target_fovs are not present: ", paste(missing_fovs, collapse = ", "))

# =========================================================
# 07 expression orientation
# =========================================================
if (!(inherits(norm, "matrix") || inherits(norm, "Matrix") || inherits(norm, "dgCMatrix"))) {
  stop("Unsupported expression object class.")
}

expr_cells_as_rows <- any(dt$cell_id %in% rownames(norm))
expr_cells_as_cols <- any(dt$cell_id %in% colnames(norm))

if (!expr_cells_as_rows && !expr_cells_as_cols) stop("Expression matrix does not match metadata cell_id.")
if (expr_cells_as_cols && !expr_cells_as_rows) norm <- t(norm)
if (!any(dt$cell_id %in% rownames(norm))) stop("Expression matrix rownames do not contain cell_id after orientation check.")

genes_present <- intersect(score_genes, colnames(norm))
save_csv_overwrite(
  data.table(score_name = score_name, input_gene = score_genes, in_expression = score_genes %in% genes_present),
  paste0("05_1b_", score_name, "_gene_presence")
)

if (length(genes_present) < min_genes_present) {
  stop("Too few genes present in expression matrix. Present genes: ", paste(genes_present, collapse = ", "))
}

# =========================================================
# 08 tumor cell-level pathway score
# =========================================================
dt_tumor <- copy(dt[get(celltype_col) == tumor_label])
if (nrow(dt_tumor) == 0) stop("No tumor cells found in metadata.")

common_tumor_cells <- intersect(dt_tumor$cell_id, rownames(norm))
dt_tumor <- dt_tumor[match(common_tumor_cells, cell_id)]
dt_tumor <- dt_tumor[!is.na(cell_id)]

expr_use <- norm[dt_tumor$cell_id, genes_present, drop = FALSE]
if (!identical(rownames(expr_use), dt_tumor$cell_id)) stop("Tumor expression row order mismatch.")

expr_long <- rbindlist(lapply(genes_present, function(g) {
  data.table(
    cell_id = dt_tumor$cell_id,
    Run_Tissue_name = dt_tumor$Run_Tissue_name,
    gene = g,
    expr_value = as.numeric(expr_use[, g])
  )
}))
expr_long[, gene_z_within_slide := zscore_safe(expr_value), by = .(Run_Tissue_name, gene)]

tumor_score_dt <- expr_long[, .(
  pathway_score = mean(gene_z_within_slide, na.rm = TRUE),
  n_genes_used = uniqueN(gene)
), by = cell_id]

save_csv_overwrite(expr_long, paste0("05_1b_", score_name, "_tumor_gene_long"))
save_csv_overwrite(tumor_score_dt, paste0("05_1b_", score_name, "_tumor_cell_score"))
save_rds_overwrite(tumor_score_dt, paste0("05_1b_", score_name, "_tumor_cell_score"))

# =========================================================
# 09 merge tumor score back
# =========================================================
dt_plot_all <- merge(dt, tumor_score_dt, by = "cell_id", all.x = TRUE)

# =========================================================
# 10 global standards for panel 2 and panel 3
# =========================================================
dt_target <- copy(dt_plot_all[FOV %in% target_fovs])
dt_target[, is_tumor := get(celltype_col) == tumor_label]

global_tumor_scores <- dt_target[is_tumor == TRUE & !is.na(pathway_score), pathway_score]
if (length(global_tumor_scores) == 0) stop("No tumor pathway scores available across selected FOVs.")

if (expr_group_method == "quantile_nonzero_tumor") {
  global_cut_high <- as.numeric(quantile(global_tumor_scores, probs = high_q, na.rm = TRUE))
} else if (expr_group_method == "median_nonzero_tumor") {
  global_cut_high <- median(global_tumor_scores, na.rm = TRUE)
} else {
  stop("Unsupported expr_group_method.")
}

global_score_min <- min(global_tumor_scores, na.rm = TRUE)
global_score_max <- max(global_tumor_scores, na.rm = TRUE)

save_csv_overwrite(
  data.table(
    score_name = score_name,
    expr_group_method = expr_group_method,
    high_q = high_q,
    global_cut_high = global_cut_high,
    global_score_min = global_score_min,
    global_score_max = global_score_max,
    n_target_fovs = length(target_fovs),
    n_target_tumor_cells = length(global_tumor_scores)
  ),
  paste0("05_1b_", score_name, "_global_plotting_standards")
)

# =========================================================
# 11 plotting function
# =========================================================
plot_one_fov <- function(fov_id, dt_all, score_name, tumor_score_dt,
                         global_cut_high, global_score_min, global_score_max) {
  dt_fov <- copy(dt_all[FOV == fov_id])
  dt_fov[, is_tumor := get(celltype_col) == tumor_label]
  dt_fov[, celltype_plot := harmonize_celltype_names(get(display_celltype_col))]
  
  n_tumor_fov <- dt_fov[is_tumor == TRUE, .N]
  if (n_tumor_fov < min_tumor_cells) {
    return(list(summary = data.table(FOV = fov_id, score_name = score_name, status = "failed",
                                     error_message = "Too few tumor cells in selected FOV."),
                p1 = NULL, p2 = NULL, p3 = NULL, p4 = NULL, p5 = NULL))
  }
  
  tumor_scores <- dt_fov[is_tumor == TRUE, pathway_score]
  tumor_score_non_na <- tumor_scores[!is.na(tumor_scores)]
  if (length(tumor_score_non_na) == 0) {
    return(list(summary = data.table(FOV = fov_id, score_name = score_name, status = "failed",
                                     error_message = "No tumor pathway score available in selected FOV."),
                p1 = NULL, p2 = NULL, p3 = NULL, p4 = NULL, p5 = NULL))
  }
  
  dt_fov[, pathway_group := fifelse(
    is_tumor == FALSE, "microenvironment",
    fifelse(is.na(pathway_score) | pathway_score < global_cut_high, "tumor_low", "tumor_high")
  )]
  
  dt_fov[, cd8_panel_group := fifelse(celltype_plot == cd8_label, "CD8", "other_cells")]
  dt_fov[, immune_panel_group := fifelse(celltype_plot %in% immune_labels, "Immune", "other_cells")]
  
  build_obj <- build_fov_polygon_dt(dt_fov)
  dt_fov <- build_obj$dt_fov
  poly_dt_plot <- build_obj$poly_dt_plot
  nn_keep_cut <- build_obj$nn_keep_cut
  
  celltype_levels_present <- unique(dt_fov$celltype_plot)
  celltype_levels_present <- names(celltype_cols)[names(celltype_cols) %in% celltype_levels_present]
  
  poly_dt_plot[, celltype_plot := factor(celltype_plot, levels = celltype_levels_present)]
  poly_dt_plot[, pathway_group := factor(pathway_group, levels = c("microenvironment", "tumor_low", "tumor_high"))]
  poly_dt_plot[, cd8_panel_group := factor(cd8_panel_group, levels = c("other_cells", "CD8"))]
  poly_dt_plot[, immune_panel_group := factor(immune_panel_group, levels = c("other_cells", "Immune"))]
  
  p3_micro <- poly_dt_plot[is_tumor == FALSE]
  p3_tumor <- poly_dt_plot[is_tumor == TRUE]
  
  # single-FOV figure
  p1_single <- ggplot(poly_dt_plot, aes(x = x, y = y, group = polygon_id, fill = celltype_plot)) +
    geom_polygon(color = "black", linewidth = polygon_linewidth) +
    scale_fill_manual(values = celltype_cols, na.value = "grey50", drop = FALSE) +
    coord_fixed() +
    labs(title = paste0("FOV ", fov_id, " | ", display_celltype_col), fill = display_celltype_col) +
    base_theme()
  
  p2_single <- ggplot(poly_dt_plot, aes(x = x, y = y, group = polygon_id, fill = pathway_group)) +
    geom_polygon(color = "black", linewidth = polygon_linewidth) +
    scale_fill_manual(values = projection_cols, drop = FALSE) +
    coord_fixed() +
    labs(title = paste0("FOV ", fov_id, " | tumor ", score_name, " high/low"), fill = "pathway_group") +
    base_theme()
  
  p3_single <- ggplot() +
    geom_polygon(data = p3_micro, aes(x = x, y = y, group = polygon_id),
                 fill = "grey65", color = "black", linewidth = polygon_linewidth) +
    geom_polygon(data = p3_tumor, aes(x = x, y = y, group = polygon_id, fill = pathway_score),
                 color = "black", linewidth = polygon_linewidth) +
    scale_fill_gradient(low = cont_low, high = cont_high,
                        limits = c(global_score_min, global_score_max),
                        oob = scales::squish, na.value = "grey65") +
    coord_fixed() +
    labs(title = paste0("FOV ", fov_id, " | tumor ", score_name, " continuous"), fill = score_name) +
    base_theme()
  
  p4_single <- ggplot(poly_dt_plot, aes(x = x, y = y, group = polygon_id, fill = cd8_panel_group)) +
    geom_polygon(color = "black", linewidth = polygon_linewidth) +
    scale_fill_manual(values = cd8_only_cols, drop = FALSE) +
    coord_fixed() +
    labs(title = paste0("FOV ", fov_id, " | CD8 only"), fill = "CD8_panel") +
    base_theme()
  
  single_fig <- cowplot::plot_grid(
    cowplot::plot_grid(p1_single, p2_single, nrow = 1, align = "hv"),
    cowplot::plot_grid(p3_single, p4_single, nrow = 1, align = "hv"),
    ncol = 1,
    rel_heights = c(1, 1)
  )
  
  fig_file <- file.path(fig_dir, paste0("05_1b_", fov_id, "_", score_name, "_4panel.pdf"))
  if (file.exists(fig_file)) file.remove(fig_file)
  ggsave(fig_file, single_fig, width = 16, height = 12, bg = plot_bg)
  
  
  # Additional immune-only projection.
  # All immune cells are pink; all other cells remain grey65.
  p5_single <- ggplot(
    poly_dt_plot,
    aes(x = x, y = y, group = polygon_id, fill = immune_panel_group)
  ) +
    geom_polygon(color = "black", linewidth = polygon_linewidth) +
    scale_fill_manual(values = immune_only_cols, drop = FALSE) +
    coord_fixed() +
    labs(
      title = paste0("FOV ", fov_id, " | all immune cells"),
      fill = "Immune_panel"
    ) +
    base_theme()
  
  immune_single_file <- file.path(
    fig_dir,
    paste0("05_1b_", fov_id, "_", score_name, "_immune_only.pdf")
  )
  if (file.exists(immune_single_file)) file.remove(immune_single_file)
  ggsave(
    immune_single_file,
    p5_single,
    width = 8,
    height = 6,
    bg = plot_bg
  )
  
  # comparison plots
  p1_compare <- ggplot(poly_dt_plot, aes(x = x, y = y, group = polygon_id, fill = celltype_plot)) +
    geom_polygon(color = "black", linewidth = polygon_linewidth) +
    scale_fill_manual(
      values = celltype_cols, na.value = "grey50", drop = FALSE,
      guide = guide_legend(ncol = 2, byrow = TRUE, override.aes = list(color = NA))
    ) +
    coord_fixed() +
    labs(title = paste0("FOV ", fov_id), fill = display_celltype_col) +
    base_theme() +
    theme(
      legend.text = element_text(size = 6.6, color = "white"),
      legend.title = element_text(size = 8.3, color = "white"),
      legend.key.height = unit(2.8, "mm"),
      legend.key.width = unit(2.8, "mm")
    )
  
  p2_compare <- ggplot(poly_dt_plot, aes(x = x, y = y, group = polygon_id, fill = pathway_group)) +
    geom_polygon(color = "black", linewidth = polygon_linewidth) +
    scale_fill_manual(values = projection_cols, drop = FALSE, guide = guide_legend(override.aes = list(color = NA))) +
    coord_fixed() +
    labs(title = NULL, fill = "pathway_group") +
    base_theme() +
    theme(plot.title = element_blank())
  
  p3_compare <- ggplot() +
    geom_polygon(data = p3_micro, aes(x = x, y = y, group = polygon_id),
                 fill = "grey65", color = "black", linewidth = polygon_linewidth) +
    geom_polygon(data = p3_tumor, aes(x = x, y = y, group = polygon_id, fill = pathway_score),
                 color = "black", linewidth = polygon_linewidth) +
    scale_fill_gradient(
      low = cont_low, high = cont_high,
      limits = c(global_score_min, global_score_max),
      oob = scales::squish, na.value = "grey65",
      guide = guide_colorbar(barheight = unit(22, "mm"), barwidth = unit(4, "mm"),
                             frame.colour = NA, ticks.colour = "white")
    ) +
    coord_fixed() +
    labs(title = NULL, fill = score_name) +
    base_theme() +
    theme(plot.title = element_blank())
  
  p4_compare <- ggplot(poly_dt_plot, aes(x = x, y = y, group = polygon_id, fill = cd8_panel_group)) +
    geom_polygon(color = "black", linewidth = polygon_linewidth) +
    scale_fill_manual(values = cd8_only_cols, drop = FALSE, guide = guide_legend(override.aes = list(color = NA))) +
    coord_fixed() +
    labs(title = NULL, fill = "CD8_panel") +
    base_theme() +
    theme(plot.title = element_blank())
  
  
  p5_compare <- ggplot(
    poly_dt_plot,
    aes(x = x, y = y, group = polygon_id, fill = immune_panel_group)
  ) +
    geom_polygon(color = "black", linewidth = polygon_linewidth) +
    scale_fill_manual(
      values = immune_only_cols,
      drop = FALSE,
      guide = guide_legend(override.aes = list(color = NA))
    ) +
    coord_fixed() +
    labs(title = NULL, fill = "Immune_panel") +
    base_theme() +
    theme(plot.title = element_blank())
  
  save_csv_overwrite(dt_fov, paste0("05_1b_", fov_id, "_", score_name, "_plotting_table"))
  save_csv_overwrite(dt_fov[, .(cell_id, x, y, nn1_dist, keep_polygon)],
                     paste0("05_1b_", fov_id, "_", score_name, "_local_density_filter_table"))
  save_csv_overwrite(poly_dt_plot, paste0("05_1b_", fov_id, "_", score_name, "_voronoi_polygon_table"))
  
  summary_dt <- data.table(
    FOV = fov_id,
    score_name = score_name,
    status = "ok",
    high_q = high_q,
    global_cut_high = global_cut_high,
    global_score_min = global_score_min,
    global_score_max = global_score_max,
    n_cells = nrow(dt_fov),
    n_cells_kept = uniqueN(dt_fov[keep_polygon == TRUE, cell_id]),
    n_tumor = sum(dt_fov$is_tumor, na.rm = TRUE),
    n_microenvironment = sum(!dt_fov$is_tumor, na.rm = TRUE),
    n_tumor_low = sum(dt_fov$pathway_group == "tumor_low", na.rm = TRUE),
    n_tumor_high = sum(dt_fov$pathway_group == "tumor_high", na.rm = TRUE),
    n_cd8 = sum(dt_fov$celltype_plot == cd8_label, na.rm = TRUE),
    n_immune = sum(dt_fov$celltype_plot %in% immune_labels, na.rm = TRUE),
    n_genes_used = unique(tumor_score_dt$n_genes_used),
    nn_keep_cut = nn_keep_cut
  )
  save_csv_overwrite(summary_dt, paste0("05_1b_", fov_id, "_", score_name, "_summary"))
  
  list(
    summary = summary_dt,
    p1 = p1_compare,
    p2 = p2_compare,
    p3 = p3_compare,
    p4 = p4_compare,
    p5 = p5_compare
  )
}

# =========================================================
# 12 run all FOVs
# =========================================================
result_list <- lapply(target_fovs, function(fov_id) {
  message("Processing FOV: ", fov_id)
  tryCatch(
    plot_one_fov(
      fov_id = fov_id,
      dt_all = dt_plot_all,
      score_name = score_name,
      tumor_score_dt = tumor_score_dt,
      global_cut_high = global_cut_high,
      global_score_min = global_score_min,
      global_score_max = global_score_max
    ),
    error = function(e) {
      list(
        summary = data.table(FOV = fov_id, score_name = score_name, status = "failed",
                             error_message = conditionMessage(e)),
        p1 = NULL, p2 = NULL, p3 = NULL, p4 = NULL, p5 = NULL
      )
    }
  )
})

summary_all <- rbindlist(lapply(result_list, function(x) x$summary), fill = TRUE)
save_csv_overwrite(summary_all, paste0("05_1b_allFOV_summary_", score_name))
save_rds_overwrite(summary_all, paste0("05_1b_allFOV_summary_", score_name))

# =========================================================
# 13 combined all-FOV panel-only figure and separate legends
# =========================================================
ok_idx <- vapply(result_list, function(x) !is.null(x$p1), logical(1))
result_ok <- result_list[ok_idx]

if (length(result_ok) > 0) {
  
  # -------------------------------------------------------
  # 13.1 panel-only rows
  # -------------------------------------------------------
  row1_plots <- lapply(result_ok, function(x) x$p1 + theme(legend.position = "none"))
  row2_plots <- lapply(result_ok, function(x) x$p2 + theme(legend.position = "none"))
  row3_plots <- lapply(result_ok, function(x) x$p3 + theme(legend.position = "none"))
  row4_plots <- lapply(result_ok, function(x) x$p4 + theme(legend.position = "none"))
  row5_plots <- lapply(result_ok, function(x) x$p5 + theme(legend.position = "none"))
  
  row1_main <- cowplot::plot_grid(plotlist = row1_plots, nrow = 1, align = "hv", axis = "tblr")
  row2_main <- cowplot::plot_grid(plotlist = row2_plots, nrow = 1, align = "hv", axis = "tblr")
  row3_main <- cowplot::plot_grid(plotlist = row3_plots, nrow = 1, align = "hv", axis = "tblr")
  row4_main <- cowplot::plot_grid(plotlist = row4_plots, nrow = 1, align = "hv", axis = "tblr")
  row5_main <- cowplot::plot_grid(plotlist = row5_plots, nrow = 1, align = "hv", axis = "tblr")
  
  row1 <- cowplot::plot_grid(
    make_row_label_plot("All cell types"),
    row1_main,
    nrow = 1,
    rel_widths = c(0.32, length(row1_plots)),
    align = "hv",
    axis = "tblr"
  )
  
  row2 <- cowplot::plot_grid(
    make_row_label_plot("High / low"),
    row2_main,
    nrow = 1,
    rel_widths = c(0.32, length(row2_plots)),
    align = "hv",
    axis = "tblr"
  )
  
  row3 <- cowplot::plot_grid(
    make_row_label_plot("Continuous"),
    row3_main,
    nrow = 1,
    rel_widths = c(0.32, length(row3_plots)),
    align = "hv",
    axis = "tblr"
  )
  
  row4 <- cowplot::plot_grid(
    make_row_label_plot("CD8 only"),
    row4_main,
    nrow = 1,
    rel_widths = c(0.32, length(row4_plots)),
    align = "hv",
    axis = "tblr"
  )
  
  row5 <- cowplot::plot_grid(
    make_row_label_plot("All immune cells"),
    row5_main,
    nrow = 1,
    rel_widths = c(0.32, length(row5_plots)),
    align = "hv",
    axis = "tblr"
  )
  
  p_panels_only <- cowplot::plot_grid(
    row1, row2, row3, row4,
    ncol = 1,
    rel_heights = c(1, 1, 1, 1),
    align = "v",
    axis = "lr"
  )
  
  panels_only_file <- file.path(
    fig_dir,
    paste0("05_1b_allFOV_", score_name, "_4row_panels_only.pdf")
  )
  if (file.exists(panels_only_file)) file.remove(panels_only_file)
  ggsave(
    filename = panels_only_file,
    plot = p_panels_only,
    width = 4.3 * length(row1_plots) + 1.2,
    height = 17.5,
    bg = plot_bg
  )
  
  
  # -------------------------------------------------------
  # 13.1b additional 5-row panel-only figure
  # -------------------------------------------------------
  # Same layout as the original 4-row combined figure,
  # with "All immune cells" added as the fifth/bottom row.
  p_panels_only_5row <- cowplot::plot_grid(
    row1, row2, row3, row4, row5,
    ncol = 1,
    rel_heights = c(1, 1, 1, 1, 1),
    align = "v",
    axis = "lr"
  )
  
  panels_only_5row_file <- file.path(
    fig_dir,
    paste0("05_1b_allFOV_", score_name, "_5row_panels_only.pdf")
  )
  
  if (file.exists(panels_only_5row_file)) {
    file.remove(panels_only_5row_file)
  }
  
  ggsave(
    filename = panels_only_5row_file,
    plot = p_panels_only_5row,
    width = 4.3 * length(row1_plots) + 1.2,
    height = 21.5,
    bg = plot_bg
  )
  
  # -------------------------------------------------------
  # 13.1c additional 4-row figure WITHOUT the continuous row
  # -------------------------------------------------------
  # Rows:
  #   1) All cell types
  #   2) pathway high / low
  #   3) CD8 only
  #   4) All immune cells
  #
  # The pathway name is explicitly shown in the row label,
  # rather than using the generic label "High / low".
  row_label_width_nocont <- 0.60
  
  row1_nocont <- cowplot::plot_grid(
    make_row_label_plot("All cell types"),
    row1_main,
    nrow = 1,
    rel_widths = c(row_label_width_nocont, length(row1_plots)),
    align = "hv",
    axis = "tblr"
  )
  
  row2_nocont <- cowplot::plot_grid(
    make_row_label_plot(
      paste0(score_name, "\nHigh / low")
    ),
    row2_main,
    nrow = 1,
    rel_widths = c(row_label_width_nocont, length(row2_plots)),
    align = "hv",
    axis = "tblr"
  )
  
  row4_nocont <- cowplot::plot_grid(
    make_row_label_plot("CD8 only"),
    row4_main,
    nrow = 1,
    rel_widths = c(row_label_width_nocont, length(row4_plots)),
    align = "hv",
    axis = "tblr"
  )
  
  row5_nocont <- cowplot::plot_grid(
    make_row_label_plot("All immune cells"),
    row5_main,
    nrow = 1,
    rel_widths = c(row_label_width_nocont, length(row5_plots)),
    align = "hv",
    axis = "tblr"
  )
  
  p_panels_only_4row_no_continuous <- cowplot::plot_grid(
    row1_nocont,
    row2_nocont,
    row4_nocont,
    row5_nocont,
    ncol = 1,
    rel_heights = c(1, 1, 1, 1),
    align = "v",
    axis = "lr"
  )
  
  panels_only_4row_no_continuous_file <- file.path(
    fig_dir,
    paste0(
      "05_1b_allFOV_",
      score_name,
      "_4row_no_continuous_panels_only.pdf"
    )
  )
  
  if (file.exists(panels_only_4row_no_continuous_file)) {
    file.remove(panels_only_4row_no_continuous_file)
  }
  
  ggsave(
    filename = panels_only_4row_no_continuous_file,
    plot = p_panels_only_4row_no_continuous,
    width = 4.3 * length(row1_plots) + 1.5,
    height = 17.5,
    bg = plot_bg
  )
  
  # -------------------------------------------------------
  # 13.2 separate legends
  # -------------------------------------------------------
  p1_legend_plot <- result_ok[[1]]$p1 +
    theme(
      legend.position = "right",
      legend.text = element_text(size = 6.6, color = "white"),
      legend.title = element_text(size = 8.3, color = "white"),
      legend.key.height = unit(2.8, "mm"),
      legend.key.width = unit(2.8, "mm")
    ) +
    guides(fill = guide_legend(ncol = 2, byrow = TRUE))
  
  p2_legend_plot <- result_ok[[1]]$p2 +
    theme(legend.position = "right")
  
  p3_legend_plot <- result_ok[[1]]$p3 +
    theme(legend.position = "right")
  
  p4_legend_plot <- result_ok[[1]]$p4 +
    theme(legend.position = "right")
  
  save_legend_only_pdf(
    p1_legend_plot,
    file.path(fig_dir, paste0("05_1b_allFOV_", score_name, "_legend_all_cell_types.pdf")),
    width = 4.4,
    height = 2.8,
    bg = plot_bg
  )
  
  save_legend_only_pdf(
    p2_legend_plot,
    file.path(fig_dir, paste0("05_1b_allFOV_", score_name, "_legend_high_low.pdf")),
    width = 2.8,
    height = 2.0,
    bg = plot_bg
  )
  
  save_legend_only_pdf(
    p3_legend_plot,
    file.path(fig_dir, paste0("05_1b_allFOV_", score_name, "_legend_continuous.pdf")),
    width = 3.0,
    height = 2.4,
    bg = plot_bg
  )
  
  save_legend_only_pdf(
    p4_legend_plot,
    file.path(fig_dir, paste0("05_1b_allFOV_", score_name, "_legend_cd8_only.pdf")),
    width = 2.6,
    height = 1.8,
    bg = plot_bg
  )
  
  # -------------------------------------------------------
  # 13.3 additional all-immune-only projection
  # -------------------------------------------------------
  # The original 4-row combined figure remains unchanged.
  # This new output is analogous to the existing CD8-only view.
  immune_plots <- lapply(
    result_ok,
    function(x) x$p5 + theme(legend.position = "none")
  )
  
  immune_main <- cowplot::plot_grid(
    plotlist = immune_plots,
    nrow = 1,
    align = "hv",
    axis = "tblr"
  )
  
  immune_row <- cowplot::plot_grid(
    make_row_label_plot("All immune cells"),
    immune_main,
    nrow = 1,
    rel_widths = c(0.32, length(immune_plots)),
    align = "hv",
    axis = "tblr"
  )
  
  immune_panels_file <- file.path(
    fig_dir,
    paste0("05_1b_allFOV_", score_name, "_immune_only_panels_only.pdf")
  )
  if (file.exists(immune_panels_file)) file.remove(immune_panels_file)
  
  ggsave(
    filename = immune_panels_file,
    plot = immune_row,
    width = 4.3 * length(immune_plots) + 1.2,
    height = 4.8,
    bg = plot_bg
  )
  
  p5_legend_plot <- result_ok[[1]]$p5 +
    theme(legend.position = "right")
  
  save_legend_only_pdf(
    p5_legend_plot,
    file.path(fig_dir, paste0("05_1b_allFOV_", score_name, "_legend_immune_only.pdf")),
    width = 2.8,
    height = 1.8,
    bg = plot_bg
  )
  
}

message("05_1b finished successfully, including immune-only, 5-row, and 4-row-no-continuous combined figures.")
message("Output data directory: ", data_dir)
message("Output figure directory: ", fig_dir)
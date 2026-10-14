rm(list = ls())
gc()
graphics.off()

options(stringsAsFactors = FALSE)
set.seed(1)

# ---- 01 paths ----
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
raw <- file.path(proj, "raw")

out_dir <- file.path(proj, "results", "03_2")
fig_dir <- file.path(out_dir, "figures")
dat_dir <- file.path(out_dir, "data")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dat_dir, recursive = TRUE, showWarnings = FALSE)

save_dat_overwrite <- function(x, filename, ext = c("csv", "rds")) {
  ext <- match.arg(ext)
  out <- file.path(dat_dir, paste0(filename, ".", ext))
  if (ext == "csv") write.csv(x, out, row.names = FALSE)
  if (ext == "rds") saveRDS(x, out)
  message("Saved data: ", out)
  invisible(out)
}

save_fig_overwrite <- function(p, filename, w = 8, h = 6, dpi = 300, ext = "pdf") {
  out <- file.path(fig_dir, paste0(filename, ".", ext))
  if (tolower(ext) == "pdf") {
    ggplot2::ggsave(out, plot = p, width = w, height = h, device = cairo_pdf)
  } else {
    ggplot2::ggsave(out, plot = p, width = w, height = h, dpi = dpi)
  }
  message("Saved figure: ", out)
  invisible(out)
}

# ---- 02 packages ----
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(ggplot2)
  library(scales)
  library(grid)
})

# ---- 03 load data ----
# Finalized local metadata integrating the original input with Yi's updated
# 172-FOV annotations; the filename is retained from the original data package.
# celltype_new: updated fine labels; BC6k_refined: mapped broad labels.
# disttonearesttumor is already in micrometres. Do not multiply by 1000 again.
# See README.md, 'Input provenance and annotation versions'.
meta <- readRDS(file.path(raw, "TNBC_6k_metadata_Bogang_selected.RDS"))
meta_dt <- as.data.table(meta)

# Supplied normalized expression; align to finalized metadata by cell_id below.
norm <- readRDS(file.path(raw, "TNBC_6k_norm_Bogang.RDS"))

sig_file <- file.path(raw, "selected_subtype_specific_genes_tumor.csv")
if (!file.exists(sig_file)) {
  stop("Signature file not found: ", sig_file)
}
sig <- fread(sig_file)

# ---- 04 checks ----
need_meta_cols <- c(
  "cell_id", "FOV", "BC6k_refined", "immunotype_BW",
  "disttonearesttumor_3group"
)
stopifnot(all(need_meta_cols %in% names(meta_dt)))
stopifnot(all(c("Gene", "Subtype") %in% names(sig)))

# ---- 05 extract expression matrix ----
get_expr <- function(norm) {
  if (is.matrix(norm) || inherits(norm, "Matrix")) return(norm)
  
  if (inherits(norm, "Seurat")) {
    assay <- if ("RNA" %in% names(norm@assays)) "RNA" else Seurat::DefaultAssay(norm)
    m <- tryCatch(
      Seurat::GetAssayData(norm, assay = assay, slot = "data"),
      error = function(e) Seurat::GetAssayData(norm, assay = assay, slot = "counts")
    )
    return(m)
  }
  
  if (is.list(norm)) {
    for (nm in c("data", "expr", "expression", "mat", "matrix", "norm", "logcounts", "counts")) {
      if (!is.null(norm[[nm]]) && (is.matrix(norm[[nm]]) || inherits(norm[[nm]], "Matrix"))) {
        return(norm[[nm]])
      }
    }
  }
  
  stop("Cannot extract expression matrix from norm.")
}

expr <- get_expr(norm)

if (!any(meta_dt$cell_id %in% colnames(expr)) && any(meta_dt$cell_id %in% rownames(expr))) {
  expr <- t(expr)
}

common_cells <- intersect(meta_dt$cell_id, colnames(expr))
common_cells <- common_cells[match(common_cells, meta_dt$cell_id)]
common_cells <- common_cells[!is.na(common_cells)]

expr <- expr[, common_cells, drop = FALSE]
meta_dt <- meta_dt[match(common_cells, meta_dt$cell_id)]

stopifnot(identical(colnames(expr), meta_dt$cell_id))

# =========================================================
# Part A. redefine immunotype_BW_top_40 using tumor-only pseudobulk
# =========================================================

# ---- 06 keep tumor cells only for scoring ----
tumor_dt <- meta_dt[BC6k_refined == "tumor" & !is.na(FOV)]

min_tumor_cells <- 30
tumor_fov_counts <- tumor_dt[, .N, by = FOV][order(-N)]
save_dat_overwrite(as.data.frame(tumor_fov_counts), "03_2_tumor_cell_counts_per_FOV", ext = "csv")

valid_fov <- tumor_fov_counts[N >= min_tumor_cells, FOV]
tumor_dt <- tumor_dt[FOV %in% valid_fov]

if (nrow(tumor_dt) == 0) {
  stop("No FOV passes min_tumor_cells threshold.")
}

tumor_cells <- tumor_dt$cell_id
expr_tumor <- expr[, tumor_cells, drop = FALSE]

stopifnot(identical(colnames(expr_tumor), tumor_dt$cell_id))

# ---- 07 prepare signature ----
sig[, Gene := as.character(Gene)]
sig[, Subtype := as.character(Subtype)]

subtype_map <- c(
  "tumor_desert"   = "ID",
  "tumor_excluded" = "IE",
  "tumor_inflamed" = "IN"
)

sig <- sig[Subtype %in% names(subtype_map)]
sig[, immunotype_new := subtype_map[Subtype]]

expr_genes <- rownames(expr_tumor)
sig[, gene_in_expr := Gene %in% expr_genes]

sig_present <- sig[gene_in_expr == TRUE]
sig_missing <- sig[gene_in_expr == FALSE]

save_dat_overwrite(as.data.frame(sig_present), "03_2_signature_genes_present", ext = "csv")
save_dat_overwrite(as.data.frame(sig_missing), "03_2_signature_genes_missing", ext = "csv")

gene_usage <- sig_present[, .N, by = .(Subtype, immunotype_new)][order(immunotype_new)]
save_dat_overwrite(as.data.frame(gene_usage), "03_2_signature_gene_usage_by_subtype", ext = "csv")

required_types <- c("ID", "IE", "IN")
usable_types <- gene_usage[N >= 3, immunotype_new]
missing_required_types <- setdiff(required_types, usable_types)

if (length(missing_required_types) > 0) {
  stop(
    "All three immunotype signatures must retain >=3 measured genes. Missing/insufficient: ",
    paste(missing_required_types, collapse = ", ")
  )
}
sig_present <- sig_present[immunotype_new %in% required_types]

# ---- 08 tumor-only FOV pseudobulk ----
fov_list <- split(seq_len(nrow(tumor_dt)), tumor_dt$FOV)

fov_expr_list <- lapply(names(fov_list), function(fv) {
  idx <- fov_list[[fv]]
  cells_use <- tumor_dt$cell_id[idx]
  mat <- expr_tumor[, cells_use, drop = FALSE]
  
  if (ncol(mat) == 1) {
    val <- as.numeric(mat[, 1])
  } else {
    val <- Matrix::rowMeans(mat)
  }
  
  data.table(FOV = fv, gene = rownames(expr_tumor), expr_mean = val)
})

fov_expr_long <- rbindlist(fov_expr_list)
save_dat_overwrite(as.data.frame(fov_expr_long), "03_2_FOV_tumor_pseudobulk_long", ext = "csv")

fov_expr_wide <- dcast(fov_expr_long, gene ~ FOV, value.var = "expr_mean")
fov_mat <- as.matrix(fov_expr_wide[, -1])
rownames(fov_mat) <- fov_expr_wide$gene

# ---- 09 z-score per gene across FOVs ----
zscore_vec <- function(x) {
  s <- sd(x, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(0, length(x)))
  (x - m) / s
}
fov_mat_z <- t(apply(fov_mat, 1, zscore_vec))

# ---- 10 compute signature scores ----
score_one_subtype <- function(type_label) {
  genes_use <- unique(sig_present[immunotype_new == type_label, Gene])
  genes_use <- intersect(genes_use, rownames(fov_mat_z))
  if (length(genes_use) == 0) return(NULL)
  
  score <- colMeans(fov_mat_z[genes_use, , drop = FALSE], na.rm = TRUE)
  
  data.table(
    FOV = names(score),
    subtype_score = as.numeric(score),
    subtype = type_label,
    n_genes_used = length(genes_use)
  )
}

score_list <- lapply(sort(unique(sig_present$immunotype_new)), score_one_subtype)
score_dt_long <- rbindlist(score_list, fill = TRUE)

save_dat_overwrite(as.data.frame(score_dt_long), "03_2_FOV_signature_scores_long", ext = "csv")
save_dat_overwrite(as.data.frame(score_dt_long), "03_2_FOV_signature_scores_long", ext = "rds")

score_dt <- dcast(score_dt_long, FOV ~ subtype, value.var = "subtype_score")
gene_n_dt <- dcast(score_dt_long, FOV ~ subtype, value.var = "n_genes_used")

setnames(
  gene_n_dt,
  old = setdiff(names(gene_n_dt), "FOV"),
  new = paste0(setdiff(names(gene_n_dt), "FOV"), "_n_genes")
)

score_dt <- merge(score_dt, gene_n_dt, by = "FOV", all = TRUE)

# ---- 11 assign new immunotype ----
score_cols <- c("ID", "IE", "IN")
if (!all(score_cols %in% names(score_dt))) {
  stop(
    "Missing required signature-score columns: ",
    paste(setdiff(score_cols, names(score_dt)), collapse = ", ")
  )
}

score_mat <- as.matrix(score_dt[, ..score_cols])
rownames(score_mat) <- score_dt$FOV

top1_idx <- apply(score_mat, 1, which.max)
top1_lab <- colnames(score_mat)[top1_idx]

top2_val <- apply(score_mat, 1, function(x) sort(x, decreasing = TRUE)[2])
top1_val <- apply(score_mat, 1, max)

score_dt[, immunotype_BW_top_40 := top1_lab]
score_dt[, score_top1 := top1_val]
score_dt[, score_top2 := top2_val]
score_dt[, score_margin := score_top1 - score_top2]

save_dat_overwrite(as.data.frame(score_dt), "03_2_FOV_scores_and_new_immunotype", ext = "csv")
save_dat_overwrite(as.data.frame(score_dt), "03_2_FOV_scores_and_new_immunotype", ext = "rds")

# ---- 12 merge new FOV label back to meta ----
meta_new <- merge(
  meta_dt,
  score_dt[, .(FOV, immunotype_BW_top_40, score_top1, score_top2, score_margin)],
  by = "FOV",
  all.x = TRUE
)

save_dat_overwrite(as.data.frame(meta_new), "03_2_meta_with_immunotype_BW_top_40", ext = "csv")
save_dat_overwrite(as.data.frame(meta_new), "03_2_meta_with_immunotype_BW_top_40", ext = "rds")

# ---- 13 optional QC: compare old vs new immunotype ----
# The legacy immunotype_BW field is not required for the canonical Top40 classification.
if ("immunotype_BW" %in% names(meta_dt)) {
  old_fov <- unique(meta_dt[, .(FOV, immunotype_BW)])
  compare_dt <- merge(
    old_fov,
    score_dt[, .(FOV, immunotype_BW_top_40)],
    by = "FOV",
    all.y = TRUE
  )
  
  save_dat_overwrite(as.data.frame(compare_dt), "03_2_old_vs_new_immunotype_FOV", ext = "csv")
  
  compare_tab <- compare_dt[, .N, by = .(immunotype_BW, immunotype_BW_top_40)]
  save_dat_overwrite(as.data.frame(compare_tab), "03_2_old_vs_new_immunotype_table", ext = "csv")
} else {
  message("Legacy immunotype_BW column not present; skipping old-vs-new QC comparison.")
}

# ---- 14 visualization of scores ----
score_plot_dt <- melt(
  score_dt[, c("FOV", score_cols), with = FALSE],
  id.vars = "FOV",
  variable.name = "Subtype",
  value.name = "Score"
)

score_plot_dt[, FOV := factor(FOV, levels = score_dt[order(immunotype_BW_top_40, -score_margin)]$FOV)]

p_score <- ggplot(score_plot_dt, aes(x = Subtype, y = FOV, fill = Score)) +
  geom_tile() +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
  labs(title = "Tumor-only FOV signature scores", x = NULL, y = "FOV") +
  theme_bw(base_size = 12) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_text(face = "bold")
  )

save_fig_overwrite(p_score, "03_2_FOV_signature_score_heatmap", w = 5, h = 10, ext = "pdf")

new_counts <- score_dt[, .N, by = immunotype_BW_top_40][order(immunotype_BW_top_40)]
save_dat_overwrite(as.data.frame(new_counts), "03_2_new_immunotype_counts", ext = "csv")

p_count <- ggplot(new_counts, aes(x = immunotype_BW_top_40, y = N, fill = immunotype_BW_top_40)) +
  geom_col(width = 0.75) +
  scale_fill_manual(values = c("IN" = "#FF0000", "IE" = "#0000FF", "ID" = "#000000")) +
  labs(title = "FOV counts by immunotype_BW_top_40", x = NULL, y = "Number of FOVs") +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold")
  )

save_fig_overwrite(p_count, "03_2_new_immunotype_counts", w = 5, h = 4, ext = "pdf")

# =========================================================
# Part B. downstream analyses EXCLUDING tumor cells from the start
# =========================================================

plot_dt2 <- as.data.table(meta_new)[
  !is.na(immunotype_BW_top_40) &
    !is.na(BC6k_refined) &
    !is.na(disttonearesttumor_3group) &
    !is.na(FOV) &
    BC6k_refined != "tumor"
]

immu_levels_new <- c("IN", "IE", "ID")
immu_levels_new <- immu_levels_new[immu_levels_new %in% unique(plot_dt2$immunotype_BW_top_40)]
if (length(immu_levels_new) == 0) immu_levels_new <- sort(unique(plot_dt2$immunotype_BW_top_40))

dist_levels <- c("<=50um", "50-100um", ">100um")
dist_levels <- dist_levels[dist_levels %in% unique(plot_dt2$disttonearesttumor_3group)]
if (length(dist_levels) == 0) dist_levels <- sort(unique(plot_dt2$disttonearesttumor_3group))

bc_levels <- c(
  "B.cell",
  "Plasma",
  "Plasmablast",
  "Plasmacytoid.dendritic.cell",
  "Conventional.dendritic.cell",
  "Macrophage",
  "Monocyte",
  "Neutrophil",
  "NK.cell",
  "T.cell.CD4",
  "T.cell.CD8",
  "T.cell.regulatory",
  "Fibroblast",
  "Endothelial",
  "Pericyte"
)
bc_levels <- bc_levels[bc_levels %in% unique(plot_dt2$BC6k_refined)]
if (length(bc_levels) == 0) bc_levels <- sort(unique(plot_dt2$BC6k_refined))

plot_dt2[, immunotype_BW_top_40 := factor(immunotype_BW_top_40, levels = immu_levels_new)]
plot_dt2[, disttonearesttumor_3group := factor(disttonearesttumor_3group, levels = dist_levels)]
plot_dt2[, BC6k_refined := factor(BC6k_refined, levels = bc_levels)]

bc_cols <- c(
  "B.cell"                      = "#F2A65A",
  "Plasma"                      = "#D98C2B",
  "Plasmablast"                 = "#B87333",
  "Plasmacytoid.dendritic.cell" = "#8C6D31",
  "Conventional.dendritic.cell" = "#8BC34A",
  "Macrophage"                  = "#6DAA2C",
  "Monocyte"                    = "#35B56A",
  "Neutrophil"                  = "#2EB7B0",
  "NK.cell"                     = "#23A6D5",
  "T.cell.CD4"                  = "#4A97C9",
  "T.cell.regulatory"           = "#7B8CC7",
  "T.cell.CD8"                  = "#D864A6",
  "Fibroblast"                  = "#B1AA34",
  "Endothelial"                 = "#C99A2A",
  "Pericyte"                    = "#8E7D22"
)

missing_cols <- setdiff(levels(plot_dt2$BC6k_refined), names(bc_cols))
if (length(missing_cols) > 0) {
  extra_cols <- rep("#BDBDBD", length(missing_cols))
  names(extra_cols) <- missing_cols
  bc_cols <- c(bc_cols, extra_cols)
}
bc_cols <- bc_cols[levels(plot_dt2$BC6k_refined)]

# ---- 15 re-run 03 analysis using non-tumor cells only ----
# pooled
pooled_counts_new <- plot_dt2[
  ,
  .N,
  by = .(immunotype_BW_top_40, disttonearesttumor_3group, BC6k_refined)
]
pooled_counts_new[, Percentage := N / sum(N), by = .(immunotype_BW_top_40, disttonearesttumor_3group)]

save_dat_overwrite(as.data.frame(pooled_counts_new), "03_2_pooled_counts_nonTumor_immunotypeTop40_distance_BC6k_refined", ext = "csv")
save_dat_overwrite(as.data.frame(pooled_counts_new), "03_2_pooled_counts_nonTumor_immunotypeTop40_distance_BC6k_refined", ext = "rds")

p_pooled_new <- ggplot(
  pooled_counts_new,
  aes(x = disttonearesttumor_3group, y = Percentage, fill = BC6k_refined)
) +
  geom_col(width = 0.85, color = NA) +
  facet_grid(. ~ immunotype_BW_top_40, scales = "fixed", space = "fixed") +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = c(0, 0)) +
  scale_fill_manual(values = bc_cols, drop = FALSE) +
  labs(
    x = NULL,
    y = "Percentage",
    fill = "Cell_Type",
    title = "Non-tumor BC6k_refined composition across tumor-distance groups",
    subtitle = "Pooled-cell percentages within each immunotype_BW_top_40"
  ) +
  theme_bw(base_size = 18) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text = element_text(face = "bold", size = 18),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.key.height = unit(0.55, "cm"),
    legend.key.width = unit(0.45, "cm"),
    plot.title = element_text(face = "bold", size = 20),
    plot.subtitle = element_text(size = 14)
  )

print(p_pooled_new)
save_fig_overwrite(p_pooled_new, "03_2_nonTumor_BC6kRefined_distance_by_immunotypeTop40_pooled", w = 13, h = 8, ext = "pdf")

# FOV-balanced
fov_counts_new <- plot_dt2[
  ,
  .N,
  by = .(immunotype_BW_top_40, FOV, disttonearesttumor_3group, BC6k_refined)
]
fov_counts_new[, FOV_total := sum(N), by = .(immunotype_BW_top_40, FOV, disttonearesttumor_3group)]
fov_counts_new[, prop_in_FOV := fifelse(FOV_total > 0, N / FOV_total, NA_real_)]

fov_summary_new <- fov_counts_new[
  ,
  .(
    mean_prop = mean(prop_in_FOV, na.rm = TRUE),
    sd_prop   = sd(prop_in_FOV, na.rm = TRUE),
    n_FOV     = uniqueN(FOV)
  ),
  by = .(immunotype_BW_top_40, disttonearesttumor_3group, BC6k_refined)
]

fov_summary_new[, Percentage := mean_prop / sum(mean_prop), by = .(immunotype_BW_top_40, disttonearesttumor_3group)]

save_dat_overwrite(as.data.frame(fov_counts_new), "03_2_nonTumor_FOVlevel_counts_immunotypeTop40_distance_BC6k_refined", ext = "csv")
save_dat_overwrite(as.data.frame(fov_summary_new), "03_2_nonTumor_FOVbalanced_summary_immunotypeTop40_distance_BC6k_refined", ext = "csv")
save_dat_overwrite(as.data.frame(fov_summary_new), "03_2_nonTumor_FOVbalanced_summary_immunotypeTop40_distance_BC6k_refined", ext = "rds")

p_fov_new <- ggplot(
  fov_summary_new,
  aes(x = disttonearesttumor_3group, y = Percentage, fill = BC6k_refined)
) +
  geom_col(width = 0.85, color = NA) +
  facet_grid(. ~ immunotype_BW_top_40, scales = "fixed", space = "fixed") +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = c(0, 0)) +
  scale_fill_manual(values = bc_cols, drop = FALSE) +
  labs(
    x = NULL,
    y = "Percentage",
    fill = "Cell_Type",
    title = "Non-tumor BC6k_refined composition across tumor-distance groups",
    subtitle = "FOV-balanced mean percentages within each immunotype_BW_top_40"
  ) +
  theme_bw(base_size = 18) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text = element_text(face = "bold", size = 18),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.key.height = unit(0.55, "cm"),
    legend.key.width = unit(0.45, "cm"),
    plot.title = element_text(face = "bold", size = 20),
    plot.subtitle = element_text(size = 14)
  )

print(p_fov_new)
save_fig_overwrite(p_fov_new, "03_2_nonTumor_BC6kRefined_distance_by_immunotypeTop40_FOVbalanced", w = 13, h = 8, ext = "pdf")

# wide tables
pooled_wide_new <- dcast(
  pooled_counts_new,
  immunotype_BW_top_40 + disttonearesttumor_3group ~ BC6k_refined,
  value.var = "Percentage",
  fill = 0
)

fov_wide_new <- dcast(
  fov_summary_new,
  immunotype_BW_top_40 + disttonearesttumor_3group ~ BC6k_refined,
  value.var = "Percentage",
  fill = 0
)

save_dat_overwrite(as.data.frame(pooled_wide_new), "03_2_nonTumor_pooled_percentage_wide_immunotypeTop40_distance_BC6k_refined", ext = "csv")
save_dat_overwrite(as.data.frame(fov_wide_new), "03_2_nonTumor_FOVbalanced_percentage_wide_immunotypeTop40_distance_BC6k_refined", ext = "csv")

# ---- 16 CD8-specific analysis using non-tumor denominator only ----
immu_cols <- c(
  "IN" = "#FF0000",
  "IE" = "#0000FF",
  "ID" = "#000000"
)
immu_cols <- immu_cols[names(immu_cols) %in% unique(as.character(plot_dt2$immunotype_BW_top_40))]

cd8_label <- "T.cell.CD8"
stopifnot(cd8_label %in% unique(as.character(plot_dt2$BC6k_refined)))

# 16.1 CD8 percentage in each FOV x distance
fov_total_cd8 <- plot_dt2[, .N, by = .(immunotype_BW_top_40, FOV, disttonearesttumor_3group)]
setnames(fov_total_cd8, "N", "total_cells")

fov_cd8 <- plot_dt2[BC6k_refined == cd8_label,
                    .N,
                    by = .(immunotype_BW_top_40, FOV, disttonearesttumor_3group)
]
setnames(fov_cd8, "N", "cd8_cells")

cd8_fov_pct2 <- merge(
  fov_total_cd8, fov_cd8,
  by = c("immunotype_BW_top_40", "FOV", "disttonearesttumor_3group"),
  all.x = TRUE
)
cd8_fov_pct2[is.na(cd8_cells), cd8_cells := 0]
cd8_fov_pct2[, cd8_percentage := cd8_cells / total_cells]

save_dat_overwrite(as.data.frame(cd8_fov_pct2), "03_2_nonTumor_CD8_percentage_by_FOV_distance", ext = "csv")
save_dat_overwrite(as.data.frame(cd8_fov_pct2), "03_2_nonTumor_CD8_percentage_by_FOV_distance", ext = "rds")

# 16.2 summary plot
cd8_summary2 <- cd8_fov_pct2[
  ,
  .(
    mean_pct = mean(cd8_percentage, na.rm = TRUE),
    sd_pct   = sd(cd8_percentage, na.rm = TRUE),
    sem_pct  = sd(cd8_percentage, na.rm = TRUE) / sqrt(.N),
    n_FOV    = .N
  ),
  by = .(immunotype_BW_top_40, disttonearesttumor_3group)
]

save_dat_overwrite(as.data.frame(cd8_summary2), "03_2_nonTumor_CD8_percentage_summary", ext = "csv")

p_cd8_2 <- ggplot(
  cd8_summary2,
  aes(x = disttonearesttumor_3group, y = mean_pct, color = immunotype_BW_top_40, group = immunotype_BW_top_40)
) +
  geom_point(size = 3) +
  geom_line(linewidth = 1) +
  geom_errorbar(
    aes(ymin = pmax(mean_pct - sem_pct, 0), ymax = mean_pct + sem_pct),
    width = 0.15, linewidth = 0.7
  ) +
  scale_color_manual(values = immu_cols) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    x = NULL,
    y = "CD8 percentage",
    title = "Non-tumor-normalized CD8 percentage across tumor-distance groups",
    subtitle = "Mean ± SEM across FOVs, grouped by immunotype_BW_top_40"
  ) +
  theme_bw(base_size = 18) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

print(p_cd8_2)
save_fig_overwrite(p_cd8_2, "03_2_nonTumor_CD8_percentage_across_distance_by_immunotypeTop40", w = 9, h = 6, ext = "pdf")

# 16.3 immunotype comparison within each distance
cd8_stats_between_immunotype2 <- rbindlist(
  lapply(levels(plot_dt2$disttonearesttumor_3group), function(dg) {
    sub <- cd8_fov_pct2[disttonearesttumor_3group == dg]
    if (length(unique(sub$immunotype_BW_top_40)) < 2) return(NULL)
    
    kw <- kruskal.test(cd8_percentage ~ immunotype_BW_top_40, data = sub)
    
    pw <- pairwise.wilcox.test(
      x = sub$cd8_percentage,
      g = sub$immunotype_BW_top_40,
      p.adjust.method = "BH",
      exact = FALSE
    )
    
    pw_mat <- pw$p.value
    if (is.null(pw_mat)) return(NULL)
    
    pw_df <- as.data.frame(as.table(pw_mat), stringsAsFactors = FALSE)
    colnames(pw_df) <- c("group1", "group2", "p_adj")
    
    pw_dt <- as.data.table(pw_df)
    pw_dt[, distance_group := dg]
    pw_dt[, kruskal_p := kw$p.value]
    pw_dt
  }),
  fill = TRUE
)

save_dat_overwrite(as.data.frame(cd8_stats_between_immunotype2), "03_2_nonTumor_CD8_stats_between_immunotype_each_distance", ext = "csv")

# 16.4 distance comparison within each immunotype
cd8_stats_between_distance2 <- rbindlist(
  lapply(levels(plot_dt2$immunotype_BW_top_40), function(ig) {
    sub <- cd8_fov_pct2[immunotype_BW_top_40 == ig]
    if (length(unique(sub$disttonearesttumor_3group)) < 2) return(NULL)
    
    kw <- kruskal.test(cd8_percentage ~ disttonearesttumor_3group, data = sub)
    
    pw <- pairwise.wilcox.test(
      x = sub$cd8_percentage,
      g = sub$disttonearesttumor_3group,
      p.adjust.method = "BH",
      exact = FALSE
    )
    
    pw_mat <- pw$p.value
    if (is.null(pw_mat)) return(NULL)
    
    pw_df <- as.data.frame(as.table(pw_mat), stringsAsFactors = FALSE)
    colnames(pw_df) <- c("group1", "group2", "p_adj")
    
    pw_dt <- as.data.table(pw_df)
    pw_dt[, immunotype_BW_top_40 := ig]
    pw_dt[, kruskal_p := kw$p.value]
    pw_dt
  }),
  fill = TRUE
)

save_dat_overwrite(as.data.frame(cd8_stats_between_distance2), "03_2_nonTumor_CD8_stats_between_distance_each_immunotype", ext = "csv")

# 16.5 <=50 / >50 raw and log2
cd8_ratio_base2 <- copy(cd8_fov_pct2)
cd8_ratio_base2[disttonearesttumor_3group == "<=50um", dist_group2 := "<=50um"]
cd8_ratio_base2[disttonearesttumor_3group %in% c("50-100um", ">100um"), dist_group2 := ">50um"]
cd8_ratio_base2 <- cd8_ratio_base2[!is.na(dist_group2)]

cd8_ratio_base2 <- cd8_ratio_base2[
  ,
  .(
    cd8_cells = sum(cd8_cells, na.rm = TRUE),
    total_cells = sum(total_cells, na.rm = TRUE)
  ),
  by = .(immunotype_BW_top_40, FOV, dist_group2)
]
cd8_ratio_base2[, cd8_percentage := fifelse(total_cells > 0, cd8_cells / total_cells, NA_real_)]

cd8_ratio_dt2 <- dcast(
  cd8_ratio_base2,
  immunotype_BW_top_40 + FOV ~ dist_group2,
  value.var = "cd8_percentage"
)

for (nm in c("<=50um", ">50um")) {
  if (!nm %in% names(cd8_ratio_dt2)) cd8_ratio_dt2[, (nm) := NA_real_]
}

nonzero_vals2 <- cd8_ratio_base2$cd8_percentage[cd8_ratio_base2$cd8_percentage > 0]
eps2 <- if (length(nonzero_vals2) > 0) min(nonzero_vals2, na.rm = TRUE) / 2 else 1e-4

cd8_ratio_dt2[, ratio_le50_over_gt50 := (`<=50um` + eps2) / (`>50um` + eps2)]
cd8_ratio_dt2[, log2_ratio := log2(ratio_le50_over_gt50)]

save_dat_overwrite(as.data.frame(cd8_ratio_dt2), "03_2_nonTumor_CD8_ratio_le50_over_gt50_by_FOV", ext = "csv")
save_dat_overwrite(as.data.frame(cd8_ratio_dt2), "03_2_nonTumor_CD8_ratio_le50_over_gt50_by_FOV", ext = "rds")

p_ratio2 <- ggplot(cd8_ratio_dt2, aes(x = immunotype_BW_top_40, y = ratio_le50_over_gt50, fill = immunotype_BW_top_40)) +
  geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.85) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.7) +
  geom_hline(yintercept = 1, linetype = 2, color = "grey40") +
  scale_fill_manual(values = immu_cols) +
  labs(
    x = NULL,
    y = "CD8% <=50um / CD8% >50um",
    title = "Non-tumor-normalized CD8 tumor-adjacent / distal ratio by immunotype_BW_top_40",
    subtitle = "Values < 1 suggest relative exclusion from tumor-adjacent region"
  ) +
  theme_bw(base_size = 18) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.title = element_text(face = "bold")
  )

print(p_ratio2)
save_fig_overwrite(p_ratio2, "03_2_nonTumor_CD8_ratio_le50_over_gt50_by_immunotypeTop40", w = 8, h = 6, ext = "pdf")

cd8_log2_dt <- cd8_ratio_dt2[is.finite(log2_ratio) & !is.na(log2_ratio)]

p_log2ratio <- ggplot(cd8_log2_dt, aes(x = immunotype_BW_top_40, y = log2_ratio, fill = immunotype_BW_top_40)) +
  geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.85) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey40") +
  scale_fill_manual(values = immu_cols) +
  labs(
    x = NULL,
    y = "log2( CD8% <=50um / CD8% >50um )",
    title = "Non-tumor-normalized CD8 tumor-adjacent / distal log2-ratio by immunotype_BW_top_40",
    subtitle = "Values < 0 suggest relative exclusion from tumor-adjacent region"
  ) +
  theme_bw(base_size = 18) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.title = element_text(face = "bold")
  )

print(p_log2ratio)
save_fig_overwrite(p_log2ratio, "03_2_nonTumor_CD8_log2ratio_le50_over_gt50_by_immunotypeTop40", w = 8, h = 6, ext = "pdf")

ratio_kw2 <- kruskal.test(ratio_le50_over_gt50 ~ immunotype_BW_top_40, data = cd8_ratio_dt2)
ratio_pw2 <- pairwise.wilcox.test(
  x = cd8_ratio_dt2$ratio_le50_over_gt50,
  g = cd8_ratio_dt2$immunotype_BW_top_40,
  p.adjust.method = "BH",
  exact = FALSE
)
ratio_pw_mat2 <- ratio_pw2$p.value
if (is.null(ratio_pw_mat2)) {
  ratio_pw_dt2 <- data.table(group1 = character(), group2 = character(), p_adj = numeric(), kruskal_p = numeric())
} else {
  ratio_pw_df2 <- as.data.frame(as.table(ratio_pw_mat2), stringsAsFactors = FALSE)
  colnames(ratio_pw_df2) <- c("group1", "group2", "p_adj")
  ratio_pw_dt2 <- as.data.table(ratio_pw_df2)
  ratio_pw_dt2[, kruskal_p := ratio_kw2$p.value]
}
save_dat_overwrite(as.data.frame(ratio_pw_dt2), "03_2_nonTumor_CD8_ratio_le50_over_gt50_stats_between_immunotype", ext = "csv")

log2ratio_kw <- kruskal.test(log2_ratio ~ immunotype_BW_top_40, data = cd8_log2_dt)
log2ratio_pw <- pairwise.wilcox.test(
  x = cd8_log2_dt$log2_ratio,
  g = cd8_log2_dt$immunotype_BW_top_40,
  p.adjust.method = "BH",
  exact = FALSE
)
log2ratio_pw_mat <- log2ratio_pw$p.value
if (is.null(log2ratio_pw_mat)) {
  log2ratio_pw_dt <- data.table(group1 = character(), group2 = character(), p_adj = numeric(), kruskal_p = numeric())
} else {
  log2ratio_pw_df <- as.data.frame(as.table(log2ratio_pw_mat), stringsAsFactors = FALSE)
  colnames(log2ratio_pw_df) <- c("group1", "group2", "p_adj")
  log2ratio_pw_dt <- as.data.table(log2ratio_pw_df)
  log2ratio_pw_dt[, kruskal_p := log2ratio_kw$p.value]
}
save_dat_overwrite(as.data.frame(log2ratio_pw_dt), "03_2_nonTumor_CD8_log2ratio_le50_over_gt50_stats_between_immunotype", ext = "csv")

# 16.6 <=50 / 50-100 raw + log2
cd8_ratio_base_50_100 <- copy(cd8_fov_pct2)[
  disttonearesttumor_3group %in% c("<=50um", "50-100um")
]

cd8_ratio_dt_50_100 <- dcast(
  cd8_ratio_base_50_100,
  immunotype_BW_top_40 + FOV ~ disttonearesttumor_3group,
  value.var = "cd8_percentage"
)

for (nm in c("<=50um", "50-100um")) {
  if (!nm %in% names(cd8_ratio_dt_50_100)) cd8_ratio_dt_50_100[, (nm) := NA_real_]
}

nonzero_vals_50_100 <- cd8_ratio_base_50_100$cd8_percentage[cd8_ratio_base_50_100$cd8_percentage > 0]
eps_50_100 <- if (length(nonzero_vals_50_100) > 0) min(nonzero_vals_50_100, na.rm = TRUE) / 2 else 1e-4

cd8_ratio_dt_50_100[, ratio_le50_over_50_100 := (`<=50um` + eps_50_100) / (`50-100um` + eps_50_100)]
cd8_ratio_dt_50_100[, log2_ratio_le50_over_50_100 := log2(ratio_le50_over_50_100)]

save_dat_overwrite(as.data.frame(cd8_ratio_dt_50_100), "03_2_nonTumor_CD8_ratio_le50_over_50_100_by_FOV", ext = "csv")
save_dat_overwrite(as.data.frame(cd8_ratio_dt_50_100), "03_2_nonTumor_CD8_ratio_le50_over_50_100_by_FOV", ext = "rds")

p_ratio_50_100 <- ggplot(
  cd8_ratio_dt_50_100,
  aes(x = immunotype_BW_top_40, y = ratio_le50_over_50_100, fill = immunotype_BW_top_40)
) +
  geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.85) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.7) +
  geom_hline(yintercept = 1, linetype = 2, color = "grey40") +
  scale_fill_manual(values = immu_cols) +
  labs(
    x = NULL,
    y = "CD8% <=50um / CD8% 50-100um",
    title = "Non-tumor-normalized CD8 tumor-adjacent / peri-tumoral ratio by immunotype_BW_top_40",
    subtitle = "Values < 1 suggest enrichment in the 50-100um shell rather than the <=50um shell"
  ) +
  theme_bw(base_size = 18) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.title = element_text(face = "bold")
  )

print(p_ratio_50_100)
save_fig_overwrite(p_ratio_50_100, "03_2_nonTumor_CD8_ratio_le50_over_50_100_by_immunotypeTop40", w = 8, h = 6, ext = "pdf")

cd8_log2_50_100_plot <- cd8_ratio_dt_50_100[is.finite(log2_ratio_le50_over_50_100) & !is.na(log2_ratio_le50_over_50_100)]

p_log2ratio_50_100 <- ggplot(
  cd8_log2_50_100_plot,
  aes(x = immunotype_BW_top_40, y = log2_ratio_le50_over_50_100, fill = immunotype_BW_top_40)
) +
  geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.85) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey40") +
  scale_fill_manual(values = immu_cols) +
  labs(
    x = NULL,
    y = "log2( CD8% <=50um / CD8% 50-100um )",
    title = "Non-tumor-normalized CD8 tumor-adjacent / peri-tumoral log2-ratio by immunotype_BW_top_40",
    subtitle = "Values < 0 suggest enrichment in the 50-100um shell rather than the <=50um shell"
  ) +
  theme_bw(base_size = 18) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.title = element_text(face = "bold")
  )

print(p_log2ratio_50_100)
save_fig_overwrite(p_log2ratio_50_100, "03_2_nonTumor_CD8_log2ratio_le50_over_50_100_by_immunotypeTop40", w = 8, h = 6, ext = "pdf")

ratio_kw_50_100 <- kruskal.test(ratio_le50_over_50_100 ~ immunotype_BW_top_40, data = cd8_ratio_dt_50_100)
ratio_pw_50_100 <- pairwise.wilcox.test(
  x = cd8_ratio_dt_50_100$ratio_le50_over_50_100,
  g = cd8_ratio_dt_50_100$immunotype_BW_top_40,
  p.adjust.method = "BH",
  exact = FALSE
)
ratio_pw_mat_50_100 <- ratio_pw_50_100$p.value
if (is.null(ratio_pw_mat_50_100)) {
  ratio_pw_dt_50_100 <- data.table(group1 = character(), group2 = character(), p_adj = numeric(), kruskal_p = numeric())
} else {
  ratio_pw_df_50_100 <- as.data.frame(as.table(ratio_pw_mat_50_100), stringsAsFactors = FALSE)
  colnames(ratio_pw_df_50_100) <- c("group1", "group2", "p_adj")
  ratio_pw_dt_50_100 <- as.data.table(ratio_pw_df_50_100)
  ratio_pw_dt_50_100[, kruskal_p := ratio_kw_50_100$p.value]
}
save_dat_overwrite(as.data.frame(ratio_pw_dt_50_100), "03_2_nonTumor_CD8_ratio_le50_over_50_100_stats_between_immunotype", ext = "csv")

log2ratio_kw_50_100 <- kruskal.test(log2_ratio_le50_over_50_100 ~ immunotype_BW_top_40, data = cd8_log2_50_100_plot)
log2ratio_pw_50_100 <- pairwise.wilcox.test(
  x = cd8_log2_50_100_plot$log2_ratio_le50_over_50_100,
  g = cd8_log2_50_100_plot$immunotype_BW_top_40,
  p.adjust.method = "BH",
  exact = FALSE
)
log2ratio_pw_mat_50_100 <- log2ratio_pw_50_100$p.value
if (is.null(log2ratio_pw_mat_50_100)) {
  log2ratio_pw_dt_50_100 <- data.table(group1 = character(), group2 = character(), p_adj = numeric(), kruskal_p = numeric())
} else {
  log2ratio_pw_df_50_100 <- as.data.frame(as.table(log2ratio_pw_mat_50_100), stringsAsFactors = FALSE)
  colnames(log2ratio_pw_df_50_100) <- c("group1", "group2", "p_adj")
  log2ratio_pw_dt_50_100 <- as.data.table(log2ratio_pw_df_50_100)
  log2ratio_pw_dt_50_100[, kruskal_p := log2ratio_kw_50_100$p.value]
}
save_dat_overwrite(as.data.frame(log2ratio_pw_dt_50_100), "03_2_nonTumor_CD8_log2ratio_le50_over_50_100_stats_between_immunotype", ext = "csv")

# ---- 17 clean mean ± SEM plots for CD8 percentage in <=50um and 50-100um (no stats on plot) ----
get_cell_pct_by_layer <- function(layer_name, cell_label) {
  fov_total <- plot_dt2[
    disttonearesttumor_3group == layer_name,
    .N,
    by = .(immunotype_BW_top_40, FOV)
  ]
  setnames(fov_total, "N", "total_cells")
  
  fov_cell <- plot_dt2[
    disttonearesttumor_3group == layer_name & BC6k_refined == cell_label,
    .N,
    by = .(immunotype_BW_top_40, FOV)
  ]
  setnames(fov_cell, "N", "cell_count")
  
  out <- merge(
    fov_total,
    fov_cell,
    by = c("immunotype_BW_top_40", "FOV"),
    all.x = TRUE
  )
  out[is.na(cell_count), cell_count := 0]
  out[, cell_percentage := cell_count / total_cells]
  out[, Cell_Type := cell_label]
  out
}

get_summary <- function(df) {
  df[
    ,
    .(
      mean_pct = mean(cell_percentage, na.rm = TRUE),
      sd_pct   = sd(cell_percentage, na.rm = TRUE),
      sem_pct  = sd(cell_percentage, na.rm = TRUE) / sqrt(.N),
      n_FOV    = .N
    ),
    by = .(immunotype_BW_top_40)
  ]
}

plot_one_layer_cd8 <- function(df_layer, layer_label, file_stub) {
  sum_dt <- get_summary(df_layer)
  
  save_dat_overwrite(as.data.frame(df_layer), paste0(file_stub, "_data"), ext = "csv")
  save_dat_overwrite(as.data.frame(sum_dt),   paste0(file_stub, "_summary"), ext = "csv")
  
  ymax_data <- max(sum_dt$mean_pct + sum_dt$sem_pct, na.rm = TRUE)
  ymax_use <- ymax_data * 1.10
  if (!is.finite(ymax_use) || ymax_use <= 0) ymax_use <- 0.05
  
  ybreaks <- pretty(c(0, ymax_use), n = 5)
  ybreaks <- ybreaks[ybreaks >= 0]
  
  p <- ggplot(
    sum_dt,
    aes(x = immunotype_BW_top_40, y = mean_pct, color = immunotype_BW_top_40, group = 1)
  ) +
    geom_line(linewidth = 0.8, color = "grey40") +
    geom_point(size = 3.2) +
    geom_errorbar(
      aes(ymin = pmax(mean_pct - sem_pct, 0), ymax = mean_pct + sem_pct),
      width = 0.10,
      linewidth = 0.75
    ) +
    scale_color_manual(values = immu_cols) +
    scale_y_continuous(
      breaks = ybreaks,
      labels = scales::label_percent(accuracy = 0.1),
      limits = c(0, ymax_use),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(
      x = NULL,
      y = "CD8 percentage",
      title = paste0("CD8 percentage in ", layer_label, " region"),
      subtitle = "Mean ± SEM across FOVs"
    ) +
    theme_bw(base_size = 15) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 17),
      plot.subtitle = element_text(size = 12),
      axis.title.y = element_text(size = 16),
      axis.text.x = element_text(size = 15),
      axis.text.y = element_text(size = 14)
    )
  
  print(p)
  save_fig_overwrite(p, file_stub, w = 5.8, h = 6, ext = "pdf")
}

cd8_le50 <- get_cell_pct_by_layer("<=50um", cd8_label)
plot_one_layer_cd8(
  df_layer   = cd8_le50,
  layer_label = "<=50um",
  file_stub   = "03_2_nonTumor_CD8_meanSEM_le50_by_immunotypeTop40"
)

cd8_50_100 <- get_cell_pct_by_layer("50-100um", cd8_label)
plot_one_layer_cd8(
  df_layer   = cd8_50_100,
  layer_label = "50-100um",
  file_stub   = "03_2_nonTumor_CD8_meanSEM_50_100_by_immunotypeTop40"
)

# ---- 18 all non-tumor cell types: one plot per cell type + one merged stats table ----
dt_all <- copy(plot_dt2)

all_cell_types <- sort(unique(as.character(dt_all$BC6k_refined)))

sanitize_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

fov_total_all <- dt_all[
  ,
  .N,
  by = .(immunotype_BW_top_40, FOV, disttonearesttumor_3group)
]
setnames(fov_total_all, c("N", "disttonearesttumor_3group"), c("total_cells", "dist_layer"))

cell_count_all <- dt_all[
  ,
  .N,
  by = .(immunotype_BW_top_40, FOV, disttonearesttumor_3group, BC6k_refined)
]
setnames(cell_count_all, c("N", "disttonearesttumor_3group"), c("cell_count", "dist_layer"))

template_dt <- fov_total_all[
  ,
  .(BC6k_refined = all_cell_types),
  by = .(immunotype_BW_top_40, FOV, dist_layer, total_cells)
]

all_celltype_long <- merge(
  template_dt,
  cell_count_all,
  by = c("immunotype_BW_top_40", "FOV", "dist_layer", "BC6k_refined"),
  all.x = TRUE
)

all_celltype_long[is.na(cell_count), cell_count := 0]
all_celltype_long[, cell_percentage := cell_count / total_cells]

save_dat_overwrite(as.data.frame(all_celltype_long), "03_2_nonTumor_all_celltypes_FOVlevel_percentage_long", ext = "csv")

all_celltype_summary <- all_celltype_long[
  ,
  .(
    mean_pct = mean(cell_percentage, na.rm = TRUE),
    sd_pct   = sd(cell_percentage, na.rm = TRUE),
    sem_pct  = sd(cell_percentage, na.rm = TRUE) / sqrt(.N),
    n_FOV    = .N
  ),
  by = .(BC6k_refined, dist_layer, immunotype_BW_top_40)
]

# p-values: immunotype within each distance
immunotype_stats_list <- list()
for (ct in all_cell_types) {
  for (dl in dist_levels) {
    sub <- all_celltype_long[BC6k_refined == ct & dist_layer == dl]
    if (nrow(sub) == 0) next
    if (length(unique(sub$immunotype_BW_top_40)) < 2) next
    
    kw <- tryCatch(
      kruskal.test(cell_percentage ~ immunotype_BW_top_40, data = sub),
      error = function(e) NULL
    )
    if (is.null(kw)) next
    
    pw <- tryCatch(
      pairwise.wilcox.test(
        x = sub$cell_percentage,
        g = sub$immunotype_BW_top_40,
        p.adjust.method = "BH",
        exact = FALSE
      ),
      error = function(e) NULL
    )
    if (is.null(pw) || is.null(pw$p.value)) next
    
    pw_df <- as.data.frame(as.table(pw$p.value), stringsAsFactors = FALSE)
    colnames(pw_df) <- c("group1", "group2", "p_adj")
    pw_dt <- as.data.table(pw_df)
    pw_dt[, BC6k_refined := ct]
    pw_dt[, dist_layer := as.character(dl)]
    pw_dt[, kruskal_p_immunotype := kw$p.value]
    immunotype_stats_list[[paste(ct, dl, sep = "__")]] <- pw_dt
  }
}
immunotype_stats_dt <- rbindlist(immunotype_stats_list, fill = TRUE)

# p-values: distance within each immunotype
distance_stats_list <- list()
for (ct in all_cell_types) {
  for (ig in c("IN", "IE", "ID")) {
    sub <- all_celltype_long[BC6k_refined == ct & immunotype_BW_top_40 == ig]
    if (nrow(sub) == 0) next
    if (length(unique(sub$dist_layer)) < 2) next
    
    kw <- tryCatch(
      kruskal.test(cell_percentage ~ dist_layer, data = sub),
      error = function(e) NULL
    )
    if (is.null(kw)) next
    
    pw <- tryCatch(
      pairwise.wilcox.test(
        x = sub$cell_percentage,
        g = sub$dist_layer,
        p.adjust.method = "BH",
        exact = FALSE
      ),
      error = function(e) NULL
    )
    if (is.null(pw) || is.null(pw$p.value)) next
    
    pw_df <- as.data.frame(as.table(pw$p.value), stringsAsFactors = FALSE)
    colnames(pw_df) <- c("group1", "group2", "p_adj")
    pw_dt <- as.data.table(pw_df)
    pw_dt[, BC6k_refined := ct]
    pw_dt[, immunotype_BW_top_40 := ig]
    pw_dt[, kruskal_p_distance := kw$p.value]
    distance_stats_list[[paste(ct, ig, sep = "__")]] <- pw_dt
  }
}
distance_stats_dt <- rbindlist(distance_stats_list, fill = TRUE)

# merged tables only
mean_wide <- dcast(
  all_celltype_summary,
  BC6k_refined + dist_layer ~ immunotype_BW_top_40,
  value.var = "mean_pct",
  fill = NA
)
setnames(mean_wide,
         old = intersect(c("IN", "IE", "ID"), names(mean_wide)),
         new = paste0("mean_", intersect(c("IN", "IE", "ID"), names(mean_wide))))

sem_wide <- dcast(
  all_celltype_summary,
  BC6k_refined + dist_layer ~ immunotype_BW_top_40,
  value.var = "sem_pct",
  fill = NA
)
setnames(sem_wide,
         old = intersect(c("IN", "IE", "ID"), names(sem_wide)),
         new = paste0("sem_", intersect(c("IN", "IE", "ID"), names(sem_wide))))

n_wide <- dcast(
  all_celltype_summary,
  BC6k_refined + dist_layer ~ immunotype_BW_top_40,
  value.var = "n_FOV",
  fill = NA
)
setnames(n_wide,
         old = intersect(c("IN", "IE", "ID"), names(n_wide)),
         new = paste0("nFOV_", intersect(c("IN", "IE", "ID"), names(n_wide))))

if (nrow(immunotype_stats_dt) > 0) {
  immunotype_stats_dt[, comparison := paste(group1, "vs", group2, sep = "_")]
  immunotype_p_wide <- dcast(
    immunotype_stats_dt,
    BC6k_refined + dist_layer + kruskal_p_immunotype ~ comparison,
    value.var = "p_adj",
    fill = NA
  )
  
  old_cols <- intersect(c("IN_vs_IE", "IN_vs_ID", "IE_vs_ID"), names(immunotype_p_wide))
  if (length(old_cols) > 0) {
    setnames(
      immunotype_p_wide,
      old = old_cols,
      new = paste0("p_adj_", old_cols)
    )
  }
} else {
  immunotype_p_wide <- unique(all_celltype_summary[, .(BC6k_refined, dist_layer)])
  immunotype_p_wide[, kruskal_p_immunotype := NA_real_]
}

merged_by_distance <- Reduce(
  function(x, y) merge(x, y, by = c("BC6k_refined", "dist_layer"), all = TRUE),
  list(mean_wide, sem_wide, n_wide, immunotype_p_wide)
)

mean_dist_wide <- dcast(
  all_celltype_summary,
  BC6k_refined + immunotype_BW_top_40 ~ dist_layer,
  value.var = "mean_pct",
  fill = NA
)
setnames(
  mean_dist_wide,
  old = intersect(dist_levels, names(mean_dist_wide)),
  new = paste0("mean_", make.names(intersect(dist_levels, names(mean_dist_wide))))
)

sem_dist_wide <- dcast(
  all_celltype_summary,
  BC6k_refined + immunotype_BW_top_40 ~ dist_layer,
  value.var = "sem_pct",
  fill = NA
)
setnames(
  sem_dist_wide,
  old = intersect(dist_levels, names(sem_dist_wide)),
  new = paste0("sem_", make.names(intersect(dist_levels, names(sem_dist_wide))))
)

n_dist_wide <- dcast(
  all_celltype_summary,
  BC6k_refined + immunotype_BW_top_40 ~ dist_layer,
  value.var = "n_FOV",
  fill = NA
)
setnames(
  n_dist_wide,
  old = intersect(dist_levels, names(n_dist_wide)),
  new = paste0("nFOV_", make.names(intersect(dist_levels, names(n_dist_wide))))
)

if (nrow(distance_stats_dt) > 0) {
  distance_stats_dt[, comparison := paste(group1, "vs", group2, sep = "_")]
  distance_p_wide <- dcast(
    distance_stats_dt,
    BC6k_refined + immunotype_BW_top_40 + kruskal_p_distance ~ comparison,
    value.var = "p_adj",
    fill = NA
  )
  
  other_cols <- setdiff(names(distance_p_wide), c("BC6k_refined", "immunotype_BW_top_40", "kruskal_p_distance"))
  if (length(other_cols) > 0) {
    setnames(
      distance_p_wide,
      old = other_cols,
      new = paste0("p_adj_", make.names(other_cols))
    )
  }
} else {
  distance_p_wide <- unique(all_celltype_summary[, .(BC6k_refined, immunotype_BW_top_40)])
  distance_p_wide[, kruskal_p_distance := NA_real_]
}

merged_by_immunotype <- Reduce(
  function(x, y) merge(x, y, by = c("BC6k_refined", "immunotype_BW_top_40"), all = TRUE),
  list(mean_dist_wide, sem_dist_wide, n_dist_wide, distance_p_wide)
)

save_dat_overwrite(as.data.frame(merged_by_distance), "03_2_nonTumor_all_celltypes_compare_immunotype_merged", ext = "csv")
save_dat_overwrite(as.data.frame(merged_by_immunotype), "03_2_nonTumor_all_celltypes_compare_distance_merged", ext = "csv")

# plots
plot_one_celltype <- function(cell_label) {
  sum_dt <- all_celltype_summary[BC6k_refined == cell_label]
  if (nrow(sum_dt) == 0) return(invisible(NULL))
  
  ymax_data <- max(sum_dt$mean_pct + sum_dt$sem_pct, na.rm = TRUE)
  ymax_use <- ymax_data * 1.10
  if (!is.finite(ymax_use) || ymax_use <= 0) ymax_use <- 0.05
  
  ybreaks <- pretty(c(0, ymax_use), n = 5)
  ybreaks <- ybreaks[ybreaks >= 0]
  
  p <- ggplot(
    sum_dt,
    aes(x = dist_layer, y = mean_pct, color = immunotype_BW_top_40, group = immunotype_BW_top_40)
  ) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 3.0) +
    geom_errorbar(
      aes(ymin = pmax(mean_pct - sem_pct, 0), ymax = mean_pct + sem_pct),
      width = 0.10,
      linewidth = 0.75
    ) +
    scale_color_manual(values = immu_cols) +
    scale_y_continuous(
      breaks = ybreaks,
      labels = scales::label_percent(accuracy = 0.1),
      limits = c(0, ymax_use),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(
      x = NULL,
      y = paste0(cell_label, " percentage"),
      title = paste0("Non-tumor ", cell_label, " percentage across tumor-distance groups"),
      subtitle = "Mean ± SEM across FOVs, grouped by immunotype_BW_top_40"
    ) +
    theme_bw(base_size = 14) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 11),
      axis.title.y = element_text(size = 14),
      axis.text.x = element_text(size = 13),
      axis.text.y = element_text(size = 12),
      legend.title = element_blank(),
      legend.text = element_text(size = 12)
    )
  
  file_stub <- paste0("03_2_nonTumor_", sanitize_name(cell_label), "_meanSEM_by_distance_and_immunotypeTop40")
  save_fig_overwrite(p, file_stub, w = 6.2, h = 5.6, ext = "pdf")
}

invisible(lapply(all_cell_types, plot_one_celltype))

message("03_2 full pipeline finished successfully.")


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

raw_dir <- file.path(proj, "raw")
in_dir  <- file.path(proj, "results", "03_2", "data")
out_dir <- file.path(proj, "results", "03_5")
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

save_plot_pdf <- function(plot_obj, filename, w = 6, h = 5.6) {
  out <- file.path(fig_dir, paste0(filename, ".pdf"))
  if (file.exists(out)) file.remove(out)
  ggplot2::ggsave(out, plot = plot_obj, width = w, height = h, device = cairo_pdf)
  message("Saved: ", out)
  invisible(out)
}

safe_tag <- function(x) {
  gsub("[^A-Za-z0-9_]+", "_", x)
}

# =========================================================
# 02 packages
# =========================================================
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(ggplot2)
  library(lme4)
  library(lmerTest)
  library(emmeans)
})

# =========================================================
# 03 user input
# =========================================================
# Keep the original broad cell-type analyses and additionally analyze
# Yi's updated CD8 subtypes stored in celltype_new.
#
# Broad analyses use BC6k_refined.
# CD8 subtype analyses use celltype_new.

broad_celltypes <- c(
  "Fibroblast",
  "Macrophage",
  "T.cell.CD4",
  "T.cell.CD8",
  "T.cell.regulatory",
  "tumor"
)

cd8_subtypes <- c(
  "cd8_cytotoxic",
  "cd8_exhausted",
  "cd8_naive",
  "cd8_tcm",
  "cd8_tem"
)

# Existing immune-regulatory/checkpoint genes plus a compact set of
# CD8 functional-state genes useful for interpreting cytotoxicity,
# exhaustion and memory state. Genes absent from the CosMx matrix are
# automatically skipped later.
target_genes <- c(
  # checkpoint / immune-regulatory axis
  "PDCD1", "CD274", "CTLA4", "LAG3", "HAVCR2", "TIGIT",
  "VSIR", "IDO1", "PVR", "NECTIN2", "LGALS9", "CD47",
  # CD8 cytotoxic / effector state
  "GZMB", "PRF1", "IFNG", "NKG7",
  # CD8 exhaustion / dysfunction state
  "TOX", "ENTPD1", "CXCL13",
  # CD8 naive / memory state
  "TCF7", "CCR7", "IL7R"
)

# FOV-level pseudobulk: minimum number of selected cells in one FOV
min_cells_per_fov <- 10

# cell-level plotting: max jitter points shown per group
max_points_per_group_for_plot <- 500

# x-axis order
immunotype_order <- c("ID", "IE", "IN")

# project colors
immunotype_cols <- c(
  "IN" = "#FF0000",
  "IE" = "#0000FF",
  "ID" = "#000000"
)

# =========================================================
# 04 load data
# =========================================================
meta_file <- file.path(in_dir, "03_2_meta_with_immunotype_BW_top_40.rds")
norm_file <- file.path(raw_dir, "TNBC_6k_norm_Bogang.RDS")

if (!file.exists(meta_file)) stop("Missing file: ", meta_file)
if (!file.exists(norm_file)) stop("Missing file: ", norm_file)

dt <- as.data.table(readRDS(meta_file))
norm <- readRDS(norm_file)

# =========================================================
# 05 checks
# =========================================================
need_cols <- c(
  "cell_id", "FOV", "Run_Tissue_name",
  "immunotype_BW_top_40", "BC6k_refined", "celltype_new"
)
miss <- setdiff(need_cols, names(dt))
if (length(miss) > 0) {
  stop("Missing required columns in metadata: ", paste(miss, collapse = ", "))
}

dt <- dt[
  !is.na(cell_id) &
    !is.na(FOV) &
    !is.na(Run_Tissue_name) &
    !is.na(immunotype_BW_top_40)
]

dt[, cell_id := as.character(cell_id)]
dt[, FOV := as.character(FOV)]
dt[, Run_Tissue_name := as.character(Run_Tissue_name)]
dt[, immunotype_BW_top_40 := as.character(immunotype_BW_top_40)]
dt[, BC6k_refined := as.character(BC6k_refined)]
dt[, celltype_new := as.character(celltype_new)]

dt <- dt[immunotype_BW_top_40 %in% immunotype_order]

# =========================================================
# 06 export available cell-type levels and define analysis sets
# =========================================================
avail_broad <- data.table(
  celltype_col = "BC6k_refined",
  celltype = sort(unique(dt$BC6k_refined))
)

avail_fine <- data.table(
  celltype_col = "celltype_new",
  celltype = sort(unique(dt$celltype_new))
)

save_dat_overwrite(
  as.data.frame(rbindlist(list(avail_broad, avail_fine), fill = TRUE)),
  "03_5_available_celltype_levels",
  ext = "csv"
)

valid_broad <- intersect(unique(broad_celltypes), avail_broad$celltype)
invalid_broad <- setdiff(unique(broad_celltypes), valid_broad)
if (length(invalid_broad) > 0) {
  warning("These broad cell types were not found and will be skipped: ",
          paste(invalid_broad, collapse = ", "))
}

valid_cd8 <- intersect(unique(cd8_subtypes), avail_fine$celltype)
invalid_cd8 <- setdiff(unique(cd8_subtypes), valid_cd8)
if (length(invalid_cd8) > 0) {
  warning("These CD8 subtypes were not found and will be skipped: ",
          paste(invalid_cd8, collapse = ", "))
}

analysis_specs <- rbindlist(list(
  data.table(
    analysis_celltype = valid_broad,
    source_col = "BC6k_refined",
    source_value = valid_broad,
    analysis_level = "broad"
  ),
  data.table(
    analysis_celltype = valid_cd8,
    source_col = "celltype_new",
    source_value = valid_cd8,
    analysis_level = "CD8_subtype"
  )
), fill = TRUE)

if (nrow(analysis_specs) == 0) {
  stop("No valid cell types found for analysis.")
}

save_dat_overwrite(
  as.data.frame(analysis_specs),
  "03_5_analysis_celltype_specs",
  ext = "csv"
)

target_genes <- unique(target_genes)

# =========================================================
# 07 extract expression matrix
# =========================================================
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

if (!any(dt$cell_id %in% colnames(expr)) && any(dt$cell_id %in% rownames(expr))) {
  expr <- t(expr)
}

common_cells <- intersect(dt$cell_id, colnames(expr))
common_cells <- common_cells[match(common_cells, dt$cell_id)]
common_cells <- common_cells[!is.na(common_cells)]

dt <- dt[match(common_cells, cell_id)]
expr <- expr[, common_cells, drop = FALSE]

if (!identical(colnames(expr), dt$cell_id)) {
  stop("Expression matrix cell order does not match metadata.")
}

valid_genes <- intersect(target_genes, rownames(expr))
invalid_genes <- setdiff(target_genes, valid_genes)

if (length(invalid_genes) > 0) {
  warning("These target_genes were not found and will be skipped: ",
          paste(invalid_genes, collapse = ", "))
}
if (length(valid_genes) == 0) {
  stop("No valid target_genes found in expression matrix.")
}

save_dat_overwrite(
  data.frame(gene = target_genes, present_in_expression = target_genes %in% valid_genes),
  "03_5_target_gene_availability",
  ext = "csv"
)

# =========================================================
# 08 helpers
# =========================================================
format_p_value <- function(p) {
  if (is.na(p)) return("p=NA")
  if (p < 1e-4) return(paste0("p=", formatC(p, format = "e", digits = 2)))
  paste0("p=", formatC(p, format = "f", digits = 4))
}

canonical_pair <- function(x) {
  x <- gsub("\\s+", "", x)
  x <- gsub("−", "-", x, fixed = TRUE)
  pair_map <- c(
    "ID-IE" = "ID_vs_IE", "IE-ID" = "ID_vs_IE",
    "ID-IN" = "ID_vs_IN", "IN-ID" = "ID_vs_IN",
    "IE-IN" = "IE_vs_IN", "IN-IE" = "IE_vs_IN"
  )
  unname(pair_map[x])
}

pairwise_p_table <- function(pair_dt) {
  out <- data.table(
    comparison = c("ID_vs_IE", "ID_vs_IN", "IE_vs_IN"),
    p_value = NA_real_
  )
  
  if (is.null(pair_dt) || nrow(pair_dt) == 0 || !"contrast" %in% names(pair_dt)) {
    return(out)
  }
  
  p_col <- intersect(c("p.value", "p_value", "p"), names(pair_dt))
  if (length(p_col) == 0) return(out)
  p_col <- p_col[1]
  
  tmp <- as.data.table(copy(pair_dt))
  tmp[, comparison := vapply(contrast, canonical_pair, character(1))]
  tmp <- tmp[!is.na(comparison)]
  
  if (nrow(tmp) > 0) {
    for (cc in out$comparison) {
      hit <- tmp[comparison == cc]
      if (nrow(hit) > 0) out[comparison == cc, p_value := as.numeric(hit[[p_col]][1])]
    }
  }
  out
}

add_pairwise_p_annotation <- function(p, plot_dt, y_col, pair_dt) {
  ptab <- pairwise_p_table(pair_dt)
  
  y <- as.numeric(plot_dt[[y_col]])
  y <- y[is.finite(y)]
  if (length(y) == 0) return(p)
  
  yr <- range(y, na.rm = TRUE)
  span <- diff(yr)
  if (!is.finite(span) || span == 0) span <- max(abs(yr), 1) * 0.15
  if (!is.finite(span) || span == 0) span <- 1
  
  y0 <- max(y, na.rm = TRUE) + 0.10 * span
  step <- 0.14 * span
  
  ann <- data.table(
    comparison = c("ID_vs_IE", "ID_vs_IN", "IE_vs_IN"),
    x1 = c(1, 1, 2),
    x2 = c(2, 3, 3),
    y = c(y0, y0 + step, y0 + 2 * step)
  )
  ann <- merge(ann, ptab, by = "comparison", all.x = TRUE, sort = FALSE)
  ann[, label := vapply(p_value, format_p_value, character(1))]
  ann[, xmid := (x1 + x2) / 2]
  tick <- 0.025 * span
  
  p +
    geom_segment(
      data = ann,
      aes(x = x1, xend = x2, y = y, yend = y),
      inherit.aes = FALSE,
      linewidth = 0.35
    ) +
    geom_segment(
      data = ann,
      aes(x = x1, xend = x1, y = y, yend = y - tick),
      inherit.aes = FALSE,
      linewidth = 0.35
    ) +
    geom_segment(
      data = ann,
      aes(x = x2, xend = x2, y = y, yend = y - tick),
      inherit.aes = FALSE,
      linewidth = 0.35
    ) +
    geom_text(
      data = ann,
      aes(x = xmid, y = y + 0.025 * span, label = label),
      inherit.aes = FALSE,
      size = 3.1,
      vjust = 0
    ) +
    coord_cartesian(
      ylim = c(yr[1], y0 + 2 * step + 0.12 * span),
      clip = "off"
    ) +
    theme(plot.margin = margin(8, 8, 8, 8))
}

make_violin_plot <- function(plot_dt, y_col, y_lab, title_txt, jitter_dt = NULL, pair_dt = NULL) {
  p <- ggplot(plot_dt, aes(x = immunotype_BW_top_40, y = .data[[y_col]], fill = immunotype_BW_top_40)) +
    geom_violin(scale = "width", trim = FALSE, color = NA, alpha = 0.75) +
    geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.35, color = "black") +
    scale_fill_manual(values = immunotype_cols, breaks = c("IN", "IE", "ID")) +
    scale_x_discrete(drop = FALSE) +
    labs(
      x = NULL,
      y = y_lab,
      title = title_txt,
      subtitle = "Pairwise P values: emmeans, BH-adjusted"
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", size = 11),
      plot.subtitle = element_text(size = 9),
      axis.text.x = element_text(size = 10, face = "bold"),
      axis.text.y = element_text(size = 9),
      axis.title = element_text(size = 10),
      legend.position = "none"
    )
  
  if (!is.null(jitter_dt) && nrow(jitter_dt) > 0) {
    p <- p + geom_jitter(
      data = jitter_dt,
      aes(x = immunotype_BW_top_40, y = .data[[y_col]]),
      width = 0.15,
      size = 0.8,
      alpha = 0.45,
      inherit.aes = FALSE
    )
  }
  
  if (!is.null(pair_dt)) {
    p <- add_pairwise_p_annotation(p, plot_dt, y_col, pair_dt)
  }
  
  p
}

# Robust summary helper used later.
make_pairwise_wide <- function(pair_dt) {
  if (is.null(pair_dt) || nrow(pair_dt) == 0) return(data.table())
  
  dt0 <- copy(as.data.table(pair_dt))
  if (!"contrast" %in% names(dt0)) return(data.table())
  
  p_col <- intersect(c("p.value", "p_value", "p"), names(dt0))
  est_col <- intersect(c("estimate", "emmean", "estimate1"), names(dt0))
  if (length(p_col) == 0 || length(est_col) == 0) return(data.table())
  p_col <- p_col[1]
  est_col <- est_col[1]
  
  dt0[, comparison := vapply(contrast, canonical_pair, character(1))]
  dt0 <- dt0[!is.na(comparison)]
  if (nrow(dt0) == 0) return(data.table())
  
  # Orient estimates to the canonical direction named in comparison.
  dt0[, contrast_clean := gsub("\\s+", "", contrast)]
  dt0[, estimate_use := as.numeric(get(est_col))]
  dt0[comparison == "ID_vs_IE" & contrast_clean == "IE-ID", estimate_use := -estimate_use]
  dt0[comparison == "ID_vs_IN" & contrast_clean == "IN-ID", estimate_use := -estimate_use]
  dt0[comparison == "IE_vs_IN" & contrast_clean == "IN-IE", estimate_use := -estimate_use]
  
  out <- data.table(
    comparison = dt0$comparison,
    estimate = dt0$estimate_use,
    p_value = as.numeric(dt0[[p_col]])
  )
  
  out[, direction := fifelse(
    comparison == "ID_vs_IE" & estimate > 0, "ID > IE",
    fifelse(
      comparison == "ID_vs_IE" & estimate < 0, "IE > ID",
      fifelse(
        comparison == "ID_vs_IN" & estimate > 0, "ID > IN",
        fifelse(
          comparison == "ID_vs_IN" & estimate < 0, "IN > ID",
          fifelse(
            comparison == "IE_vs_IN" & estimate > 0, "IE > IN",
            fifelse(comparison == "IE_vs_IN" & estimate < 0, "IN > IE", "tie")
          )
        )
      )
    )
  )]
  
  keep_contrasts <- c("ID_vs_IE", "ID_vs_IN", "IE_vs_IN")
  wide_est <- dcast(out, . ~ comparison, value.var = "estimate")
  wide_p   <- dcast(out, . ~ comparison, value.var = "p_value")
  wide_dir <- dcast(out, . ~ comparison, value.var = "direction")
  
  for (cc in keep_contrasts) {
    if (!cc %in% names(wide_est)) wide_est[, (cc) := NA_real_]
    if (!cc %in% names(wide_p))   wide_p[, (cc) := NA_real_]
    if (!cc %in% names(wide_dir)) wide_dir[, (cc) := NA_character_]
  }
  
  setcolorder(wide_est, c(".", keep_contrasts))
  setcolorder(wide_p,   c(".", keep_contrasts))
  setcolorder(wide_dir, c(".", keep_contrasts))
  
  setnames(wide_est, keep_contrasts, paste0(keep_contrasts, "_estimate"))
  setnames(wide_p,   keep_contrasts, paste0(keep_contrasts, "_p"))
  setnames(wide_dir, keep_contrasts, paste0(keep_contrasts, "_direction"))
  
  cbind(wide_est, wide_p[, -1], wide_dir[, -1])
}

extract_overall_p <- function(anova_dt) {
  if (is.null(anova_dt) || nrow(anova_dt) == 0) return(NA_real_)
  dt0 <- as.data.table(anova_dt)
  p_col <- intersect(c("Pr(>F)", "Pr..F..", "p.value", "p_value"), names(dt0))
  if (length(p_col) == 0) return(NA_real_)
  p_col <- p_col[1]
  if ("term" %in% names(dt0)) {
    hit <- dt0[grepl("immunotype_BW_top_40", term, fixed = TRUE)]
    if (nrow(hit) > 0) return(as.numeric(hit[[p_col]][1]))
  }
  as.numeric(dt0[[p_col]][1])
}

make_rank_string <- function(mean_id, mean_ie, mean_in) {
  vals <- c(ID = mean_id, IE = mean_ie, IN = mean_in)
  vals <- vals[!is.na(vals)]
  if (length(vals) == 0) return(NA_character_)
  paste(names(sort(vals, decreasing = TRUE)), collapse = " > ")
}

# =========================================================
# 09 batch run: broad cell types + CD8 subtypes x target genes
# =========================================================
run_log <- list()

for (ii in seq_len(nrow(analysis_specs))) {
  selected_celltype <- analysis_specs$analysis_celltype[ii]
  source_col <- analysis_specs$source_col[ii]
  source_value <- analysis_specs$source_value[ii]
  analysis_level <- analysis_specs$analysis_level[ii]
  
  dt_sel0 <- copy(dt)[as.character(get(source_col)) == source_value]
  
  if (nrow(dt_sel0) == 0) {
    run_log[[length(run_log) + 1]] <- data.table(
      celltype = selected_celltype,
      analysis_level = analysis_level,
      gene = NA_character_,
      analysis = NA_character_,
      status = "skip_no_cells"
    )
    next
  }
  
  dt_sel0[, immunotype_BW_top_40 := factor(immunotype_BW_top_40, levels = immunotype_order)]
  
  for (target_gene in valid_genes) {
    tag <- paste0(safe_tag(selected_celltype), "__", safe_tag(target_gene))
    
    dt_sel <- copy(dt_sel0)
    dt_sel[, gene_expr := as.numeric(expr[target_gene, dt_sel$cell_id])]
    
    save_dat_overwrite(
      as.data.frame(dt_sel),
      paste0("03_5_", tag, "_cell_level_data"),
      ext = "csv"
    )
    
    # -----------------------------------------------------
    # 09A FOV-level pseudobulk -- primary analysis
    # -----------------------------------------------------
    fov_dt <- dt_sel[, .(
      n_cells_selected = .N,
      gene_expr_mean = mean(gene_expr, na.rm = TRUE),
      gene_expr_median = median(gene_expr, na.rm = TRUE)
    ), by = .(FOV, Run_Tissue_name, immunotype_BW_top_40)]
    
    fov_dt <- fov_dt[n_cells_selected >= min_cells_per_fov]
    
    if (nrow(fov_dt) > 0 && uniqueN(fov_dt$immunotype_BW_top_40) >= 2) {
      fov_dt[, immunotype_BW_top_40 := factor(immunotype_BW_top_40, levels = immunotype_order)]
      fov_dt[, Run_Tissue_name := factor(Run_Tissue_name)]
      
      save_dat_overwrite(
        as.data.frame(fov_dt),
        paste0("03_5_", tag, "_FOV_pseudobulk"),
        ext = "csv"
      )
      
      n_slide_fov <- uniqueN(fov_dt$Run_Tissue_name)
      fit_fov <- tryCatch({
        if (n_slide_fov >= 2) {
          lmerTest::lmer(
            gene_expr_mean ~ immunotype_BW_top_40 + (1 | Run_Tissue_name),
            data = fov_dt,
            REML = FALSE
          )
        } else {
          lm(gene_expr_mean ~ immunotype_BW_top_40, data = fov_dt)
        }
      }, error = function(e) NULL)
      
      fov_anova <- data.frame()
      fov_pair <- data.table()
      if (!is.null(fit_fov)) {
        fov_anova <- tryCatch(as.data.frame(anova(fit_fov)), error = function(e) data.frame())
        if (nrow(fov_anova) > 0) fov_anova$term <- rownames(fov_anova)
        save_dat_overwrite(
          fov_anova,
          paste0("03_5_", tag, "_FOV_model_anova"),
          ext = "csv"
        )
        
        fov_pair <- tryCatch({
          as.data.table(summary(
            pairs(emmeans::emmeans(fit_fov, specs = ~ immunotype_BW_top_40), adjust = "BH")
          ))
        }, error = function(e) data.table())
        save_dat_overwrite(
          as.data.frame(fov_pair),
          paste0("03_5_", tag, "_FOV_pairwise"),
          ext = "csv"
        )
      }
      
      fov_group_summary <- fov_dt[, .(
        n_FOV = .N,
        mean_expr = mean(gene_expr_mean, na.rm = TRUE),
        median_expr = median(gene_expr_mean, na.rm = TRUE),
        sd_expr = sd(gene_expr_mean, na.rm = TRUE)
      ), by = immunotype_BW_top_40]
      save_dat_overwrite(
        as.data.frame(fov_group_summary),
        paste0("03_5_", tag, "_FOV_group_summary"),
        ext = "csv"
      )
      
      p_fov <- make_violin_plot(
        plot_dt = fov_dt,
        y_col = "gene_expr_mean",
        y_lab = paste0(target_gene, " expression"),
        title_txt = paste0(target_gene, " in ", selected_celltype, " (FOV pseudobulk)"),
        jitter_dt = fov_dt,
        pair_dt = fov_pair
      )
      
      save_plot_pdf(
        p_fov,
        paste0("03_5_", tag, "_violin_FOV_pseudobulk"),
        w = 5.8, h = 5.8
      )
      
      run_log[[length(run_log) + 1]] <- data.table(
        celltype = selected_celltype,
        analysis_level = analysis_level,
        gene = target_gene,
        analysis = "FOV_pseudobulk",
        status = "ok"
      )
      
    } else {
      run_log[[length(run_log) + 1]] <- data.table(
        celltype = selected_celltype,
        analysis_level = analysis_level,
        gene = target_gene,
        analysis = "FOV_pseudobulk",
        status = "skip_insufficient_FOV_data"
      )
    }
    
    # -----------------------------------------------------
    # 09B cell-level mixed model -- secondary analysis
    # -----------------------------------------------------
    cell_dt <- copy(dt_sel)
    cell_dt[, Run_Tissue_name := factor(Run_Tissue_name)]
    cell_dt[, FOV := factor(FOV)]
    cell_dt[, immunotype_BW_top_40 := factor(immunotype_BW_top_40, levels = immunotype_order)]
    
    if (nrow(cell_dt) > 0 && uniqueN(cell_dt$immunotype_BW_top_40) >= 2) {
      cell_group_summary <- cell_dt[, .(
        n_cells = .N,
        n_FOV = uniqueN(FOV),
        mean_expr = mean(gene_expr, na.rm = TRUE),
        median_expr = median(gene_expr, na.rm = TRUE),
        sd_expr = sd(gene_expr, na.rm = TRUE)
      ), by = immunotype_BW_top_40]
      save_dat_overwrite(
        as.data.frame(cell_group_summary),
        paste0("03_5_", tag, "_cell_group_summary"),
        ext = "csv"
      )
      
      n_slide_cell <- uniqueN(cell_dt$Run_Tissue_name)
      n_fov_cell <- uniqueN(cell_dt$FOV)
      
      fit_cell <- tryCatch({
        if (n_slide_cell >= 2 && n_fov_cell >= 2) {
          lmerTest::lmer(
            gene_expr ~ immunotype_BW_top_40 + (1 | Run_Tissue_name) + (1 | FOV),
            data = cell_dt,
            REML = FALSE
          )
        } else if (n_fov_cell >= 2) {
          lmerTest::lmer(
            gene_expr ~ immunotype_BW_top_40 + (1 | FOV),
            data = cell_dt,
            REML = FALSE
          )
        } else {
          lm(gene_expr ~ immunotype_BW_top_40, data = cell_dt)
        }
      }, error = function(e) NULL)
      
      cell_anova <- data.frame()
      cell_pair <- data.table()
      if (!is.null(fit_cell)) {
        cell_anova <- tryCatch(as.data.frame(anova(fit_cell)), error = function(e) data.frame())
        if (nrow(cell_anova) > 0) cell_anova$term <- rownames(cell_anova)
        save_dat_overwrite(
          cell_anova,
          paste0("03_5_", tag, "_cell_model_anova"),
          ext = "csv"
        )
        
        cell_pair <- tryCatch({
          as.data.table(summary(
            pairs(emmeans::emmeans(fit_cell, specs = ~ immunotype_BW_top_40), adjust = "BH")
          ))
        }, error = function(e) data.table())
        save_dat_overwrite(
          as.data.frame(cell_pair),
          paste0("03_5_", tag, "_cell_pairwise"),
          ext = "csv"
        )
      }
      
      set.seed(1)
      cell_jitter_dt <- cell_dt[, {
        if (.N <= max_points_per_group_for_plot) .SD else .SD[sample(.N, max_points_per_group_for_plot)]
      }, by = immunotype_BW_top_40]
      
      p_cell <- make_violin_plot(
        plot_dt = cell_dt,
        y_col = "gene_expr",
        y_lab = paste0(target_gene, " expression"),
        title_txt = paste0(target_gene, " in ", selected_celltype, " (cell-level mixed model)"),
        jitter_dt = cell_jitter_dt,
        pair_dt = cell_pair
      )
      
      save_plot_pdf(
        p_cell,
        paste0("03_5_", tag, "_violin_cell_level"),
        w = 5.8, h = 5.8
      )
      
      run_log[[length(run_log) + 1]] <- data.table(
        celltype = selected_celltype,
        analysis_level = analysis_level,
        gene = target_gene,
        analysis = "cell_level_mixed_model",
        status = "ok"
      )
      
    } else {
      run_log[[length(run_log) + 1]] <- data.table(
        celltype = selected_celltype,
        analysis_level = analysis_level,
        gene = target_gene,
        analysis = "cell_level_mixed_model",
        status = "skip_insufficient_cell_data"
      )
    }
  }
}

# =========================================================
# 10 batch run log
# =========================================================
if (length(run_log) > 0) {
  run_log_dt <- rbindlist(run_log, fill = TRUE)
  save_dat_overwrite(as.data.frame(run_log_dt), "03_5_batch_run_log", ext = "csv")
}

message("03_5 broad + CD8-subtype batch finished successfully.")

# =========================================================
# 11 extra module: whole-FOV bulk expression
# ignore cell type, use all cells in each FOV
# =========================================================
for (target_gene in valid_genes) {
  tag_bulk <- paste0("FOV_bulk__", safe_tag(target_gene))
  
  fov_bulk_dt <- dt[, .(
    n_cells_all = .N,
    gene_expr_mean = mean(as.numeric(expr[target_gene, cell_id]), na.rm = TRUE),
    gene_expr_median = median(as.numeric(expr[target_gene, cell_id]), na.rm = TRUE)
  ), by = .(FOV, Run_Tissue_name, immunotype_BW_top_40)]
  
  fov_bulk_dt <- fov_bulk_dt[!is.na(immunotype_BW_top_40)]
  fov_bulk_dt <- fov_bulk_dt[immunotype_BW_top_40 %in% immunotype_order]
  
  if (nrow(fov_bulk_dt) == 0 || uniqueN(fov_bulk_dt$immunotype_BW_top_40) < 2) next
  
  fov_bulk_dt[, immunotype_BW_top_40 := factor(immunotype_BW_top_40, levels = immunotype_order)]
  fov_bulk_dt[, Run_Tissue_name := factor(Run_Tissue_name)]
  
  save_dat_overwrite(
    as.data.frame(fov_bulk_dt),
    paste0("03_5_", tag_bulk, "_data"),
    ext = "csv"
  )
  
  n_slide_bulk <- uniqueN(fov_bulk_dt$Run_Tissue_name)
  fit_bulk <- tryCatch({
    if (n_slide_bulk >= 2) {
      lmerTest::lmer(
        gene_expr_mean ~ immunotype_BW_top_40 + (1 | Run_Tissue_name),
        data = fov_bulk_dt,
        REML = FALSE
      )
    } else {
      lm(gene_expr_mean ~ immunotype_BW_top_40, data = fov_bulk_dt)
    }
  }, error = function(e) NULL)
  
  bulk_anova <- data.frame()
  bulk_pair <- data.table()
  if (!is.null(fit_bulk)) {
    bulk_anova <- tryCatch(as.data.frame(anova(fit_bulk)), error = function(e) data.frame())
    if (nrow(bulk_anova) > 0) bulk_anova$term <- rownames(bulk_anova)
    save_dat_overwrite(
      bulk_anova,
      paste0("03_5_", tag_bulk, "_anova"),
      ext = "csv"
    )
    
    bulk_pair <- tryCatch({
      as.data.table(summary(
        pairs(emmeans::emmeans(fit_bulk, specs = ~ immunotype_BW_top_40), adjust = "BH")
      ))
    }, error = function(e) data.table())
    save_dat_overwrite(
      as.data.frame(bulk_pair),
      paste0("03_5_", tag_bulk, "_pairwise"),
      ext = "csv"
    )
  }
  
  bulk_group_summary <- fov_bulk_dt[, .(
    n_FOV = .N,
    mean_expr = mean(gene_expr_mean, na.rm = TRUE),
    median_expr = median(gene_expr_mean, na.rm = TRUE),
    sd_expr = sd(gene_expr_mean, na.rm = TRUE)
  ), by = immunotype_BW_top_40]
  
  save_dat_overwrite(
    as.data.frame(bulk_group_summary),
    paste0("03_5_", tag_bulk, "_group_summary"),
    ext = "csv"
  )
  
  p_bulk <- make_violin_plot(
    plot_dt = fov_bulk_dt,
    y_col = "gene_expr_mean",
    y_lab = paste0(target_gene, " expression"),
    title_txt = paste0(target_gene, " across immunotype_BW_top_40 (whole-FOV bulk)"),
    jitter_dt = fov_bulk_dt,
    pair_dt = bulk_pair
  )
  
  save_plot_pdf(
    p_bulk,
    paste0("03_5_", tag_bulk, "_violin"),
    w = 5.8, h = 5.8
  )
}

# =========================================================
# 12 compact summary tables across all analyzed cell types and genes
# =========================================================
summary_rows <- list()

for (ii in seq_len(nrow(analysis_specs))) {
  selected_celltype <- analysis_specs$analysis_celltype[ii]
  analysis_level <- analysis_specs$analysis_level[ii]
  
  for (target_gene in valid_genes) {
    tag <- paste0(safe_tag(selected_celltype), "__", safe_tag(target_gene))
    
    # -----------------------------------------------------
    # FOV pseudobulk summary
    # -----------------------------------------------------
    fov_group_file <- file.path(dat_dir, paste0("03_5_", tag, "_FOV_group_summary.csv"))
    fov_anova_file <- file.path(dat_dir, paste0("03_5_", tag, "_FOV_model_anova.csv"))
    fov_pair_file  <- file.path(dat_dir, paste0("03_5_", tag, "_FOV_pairwise.csv"))
    
    if (file.exists(fov_group_file)) {
      grp <- as.data.table(read.csv(fov_group_file, check.names = FALSE))
      anv <- if (file.exists(fov_anova_file)) as.data.table(read.csv(fov_anova_file, check.names = FALSE)) else data.table()
      prw <- if (file.exists(fov_pair_file))  as.data.table(read.csv(fov_pair_file, check.names = FALSE))  else data.table()
      
      mean_id <- grp[immunotype_BW_top_40 == "ID", mean_expr]
      mean_ie <- grp[immunotype_BW_top_40 == "IE", mean_expr]
      mean_in <- grp[immunotype_BW_top_40 == "IN", mean_expr]
      med_id <- grp[immunotype_BW_top_40 == "ID", median_expr]
      med_ie <- grp[immunotype_BW_top_40 == "IE", median_expr]
      med_in <- grp[immunotype_BW_top_40 == "IN", median_expr]
      
      pw_wide <- make_pairwise_wide(prw)
      if (nrow(pw_wide) == 0) pw_wide <- data.table(
        ID_vs_IE_estimate = NA_real_, ID_vs_IN_estimate = NA_real_, IE_vs_IN_estimate = NA_real_,
        ID_vs_IE_p = NA_real_, ID_vs_IN_p = NA_real_, IE_vs_IN_p = NA_real_,
        ID_vs_IE_direction = NA_character_, ID_vs_IN_direction = NA_character_, IE_vs_IN_direction = NA_character_
      )
      
      summary_rows[[length(summary_rows) + 1]] <- cbind(
        data.table(
          celltype = selected_celltype,
          analysis_level = analysis_level,
          gene = target_gene,
          analysis = "FOV_pseudobulk",
          mean_ID = ifelse(length(mean_id) > 0, mean_id, NA_real_),
          mean_IE = ifelse(length(mean_ie) > 0, mean_ie, NA_real_),
          mean_IN = ifelse(length(mean_in) > 0, mean_in, NA_real_),
          median_ID = ifelse(length(med_id) > 0, med_id, NA_real_),
          median_IE = ifelse(length(med_ie) > 0, med_ie, NA_real_),
          median_IN = ifelse(length(med_in) > 0, med_in, NA_real_),
          rank_high_to_low = make_rank_string(
            ifelse(length(mean_id) > 0, mean_id, NA_real_),
            ifelse(length(mean_ie) > 0, mean_ie, NA_real_),
            ifelse(length(mean_in) > 0, mean_in, NA_real_)
          ),
          overall_p = extract_overall_p(anv)
        ),
        pw_wide
      )
    }
    
    # -----------------------------------------------------
    # cell-level mixed model summary
    # -----------------------------------------------------
    cell_group_file <- file.path(dat_dir, paste0("03_5_", tag, "_cell_group_summary.csv"))
    cell_anova_file <- file.path(dat_dir, paste0("03_5_", tag, "_cell_model_anova.csv"))
    cell_pair_file  <- file.path(dat_dir, paste0("03_5_", tag, "_cell_pairwise.csv"))
    
    if (file.exists(cell_group_file)) {
      grp <- as.data.table(read.csv(cell_group_file, check.names = FALSE))
      anv <- if (file.exists(cell_anova_file)) as.data.table(read.csv(cell_anova_file, check.names = FALSE)) else data.table()
      prw <- if (file.exists(cell_pair_file))  as.data.table(read.csv(cell_pair_file, check.names = FALSE))  else data.table()
      
      mean_id <- grp[immunotype_BW_top_40 == "ID", mean_expr]
      mean_ie <- grp[immunotype_BW_top_40 == "IE", mean_expr]
      mean_in <- grp[immunotype_BW_top_40 == "IN", mean_expr]
      med_id <- grp[immunotype_BW_top_40 == "ID", median_expr]
      med_ie <- grp[immunotype_BW_top_40 == "IE", median_expr]
      med_in <- grp[immunotype_BW_top_40 == "IN", median_expr]
      
      pw_wide <- make_pairwise_wide(prw)
      if (nrow(pw_wide) == 0) pw_wide <- data.table(
        ID_vs_IE_estimate = NA_real_, ID_vs_IN_estimate = NA_real_, IE_vs_IN_estimate = NA_real_,
        ID_vs_IE_p = NA_real_, ID_vs_IN_p = NA_real_, IE_vs_IN_p = NA_real_,
        ID_vs_IE_direction = NA_character_, ID_vs_IN_direction = NA_character_, IE_vs_IN_direction = NA_character_
      )
      
      summary_rows[[length(summary_rows) + 1]] <- cbind(
        data.table(
          celltype = selected_celltype,
          analysis_level = analysis_level,
          gene = target_gene,
          analysis = "cell_level_mixed_model",
          mean_ID = ifelse(length(mean_id) > 0, mean_id, NA_real_),
          mean_IE = ifelse(length(mean_ie) > 0, mean_ie, NA_real_),
          mean_IN = ifelse(length(mean_in) > 0, mean_in, NA_real_),
          median_ID = ifelse(length(med_id) > 0, med_id, NA_real_),
          median_IE = ifelse(length(med_ie) > 0, med_ie, NA_real_),
          median_IN = ifelse(length(med_in) > 0, med_in, NA_real_),
          rank_high_to_low = make_rank_string(
            ifelse(length(mean_id) > 0, mean_id, NA_real_),
            ifelse(length(mean_ie) > 0, mean_ie, NA_real_),
            ifelse(length(mean_in) > 0, mean_in, NA_real_)
          ),
          overall_p = extract_overall_p(anv)
        ),
        pw_wide
      )
    }
  }
}

if (length(summary_rows) > 0) {
  summary_dt_all <- rbindlist(summary_rows, fill = TRUE)
  
  setcolorder(
    summary_dt_all,
    c(
      "celltype", "analysis_level", "gene", "analysis",
      "mean_ID", "mean_IE", "mean_IN",
      "median_ID", "median_IE", "median_IN",
      "rank_high_to_low", "overall_p",
      "ID_vs_IE_estimate", "ID_vs_IE_direction", "ID_vs_IE_p",
      "ID_vs_IN_estimate", "ID_vs_IN_direction", "ID_vs_IN_p",
      "IE_vs_IN_estimate", "IE_vs_IN_direction", "IE_vs_IN_p"
    )
  )
  
  save_dat_overwrite(
    as.data.frame(summary_dt_all),
    "03_5_summary_statistics_all_inputs",
    ext = "csv"
  )
  save_dat_overwrite(
    as.data.frame(summary_dt_all),
    "03_5_summary_statistics_all_inputs",
    ext = "rds"
  )
  
  # Separate CD8-subtype-only summary for easier downstream screening.
  cd8_summary <- summary_dt_all[analysis_level == "CD8_subtype"]
  if (nrow(cd8_summary) > 0) {
    save_dat_overwrite(
      as.data.frame(cd8_summary),
      "03_5_summary_statistics_CD8_subtypes",
      ext = "csv"
    )
  }
}

message("03_5 summary tables finished successfully.")
message("03_5 whole-FOV bulk module finished successfully.")
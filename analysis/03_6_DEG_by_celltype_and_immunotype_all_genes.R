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
out_dir <- file.path(proj, "results", "03_6")
dat_dir <- file.path(out_dir, "data")

dir.create(dat_dir, recursive = TRUE, showWarnings = FALSE)

save_csv_overwrite <- function(x, filename) {
  out <- file.path(dat_dir, paste0(filename, ".csv"))
  if (file.exists(out)) file.remove(out)
  write.csv(x, out, row.names = FALSE)
  message("Saved: ", out)
  invisible(out)
}

save_rds_overwrite <- function(x, filename) {
  out <- file.path(dat_dir, paste0(filename, ".rds"))
  if (file.exists(out)) file.remove(out)
  saveRDS(x, out)
  message("Saved: ", out)
  invisible(out)
}

safe_tag <- function(x) gsub("[^A-Za-z0-9_]+", "_", x)

# =========================================================
# 02 packages
# =========================================================
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(limma)
})

# =========================================================
# 03 user settings
# =========================================================
broad_celltype_col <- "BC6k_refined"
broad_celltypes <- NULL

fine_celltype_col <- "celltype_new"
cd8_subtypes <- c(
  "cd8_cytotoxic",
  "cd8_exhausted",
  "cd8_naive",
  "cd8_tcm",
  "cd8_tem"
)

immunotype_col <- "immunotype_BW_top_40"
immunotype_order <- c("ID", "IE", "IN")

min_cells_per_fov <- 10
min_fov_per_group <- 3

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
need_cols <- c("cell_id", "FOV", "Run_Tissue_name", broad_celltype_col, immunotype_col)
miss <- setdiff(need_cols, names(dt))
if (length(miss) > 0) stop("Missing required metadata columns: ", paste(miss, collapse = ", "))

dt <- dt[
  !is.na(cell_id) &
    !is.na(FOV) &
    !is.na(Run_Tissue_name) &
    !is.na(get(immunotype_col))
]

dt[, cell_id := as.character(cell_id)]
dt[, FOV := as.character(FOV)]
dt[, Run_Tissue_name := as.character(Run_Tissue_name)]
dt[, (immunotype_col) := as.character(get(immunotype_col))]
dt <- dt[get(immunotype_col) %in% immunotype_order]

# =========================================================
# 06 extract expression matrix
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

if (is.null(rownames(expr)) || is.null(colnames(expr))) {
  stop("Expression matrix must have both rownames and colnames.")
}

n_cell_match_col <- sum(dt$cell_id %in% colnames(expr))
n_cell_match_row <- sum(dt$cell_id %in% rownames(expr))

if (n_cell_match_col == 0 && n_cell_match_row > 0) {
  expr <- t(expr)
} else if (n_cell_match_col == 0 && n_cell_match_row == 0) {
  stop("No metadata cell_id matches expression matrix rownames or colnames.")
}

common_cells <- unique(dt$cell_id[dt$cell_id %in% colnames(expr)])
dt <- dt[match(common_cells, cell_id)]
expr <- expr[, common_cells, drop = FALSE]

if (!identical(colnames(expr), dt$cell_id)) {
  stop("Expression matrix cell order does not match metadata after alignment.")
}

message("Cells used: ", nrow(dt))
message("Genes tested: ", nrow(expr))
message("FOVs used: ", uniqueN(dt$FOV))

# =========================================================
# 07 define analysis cell-type sets
# =========================================================
if (is.null(broad_celltypes)) {
  broad_celltypes <- sort(unique(as.character(dt[[broad_celltype_col]])))
  broad_celltypes <- broad_celltypes[!is.na(broad_celltypes) & nzchar(broad_celltypes)]
}

analysis_sets <- list()

for (ct in broad_celltypes) {
  analysis_sets[[paste0("broad__", ct)]] <- list(
    annotation_level = "broad",
    celltype_col = broad_celltype_col,
    celltype = ct
  )
}

if (fine_celltype_col %in% names(dt)) {
  available_fine <- sort(unique(as.character(dt[[fine_celltype_col]])))
  valid_cd8_subtypes <- intersect(cd8_subtypes, available_fine)
  
  for (ct in valid_cd8_subtypes) {
    analysis_sets[[paste0("fine__", ct)]] <- list(
      annotation_level = "fine_CD8",
      celltype_col = fine_celltype_col,
      celltype = ct
    )
  }
  
  missing_cd8_subtypes <- setdiff(cd8_subtypes, valid_cd8_subtypes)
  if (length(missing_cd8_subtypes) > 0) {
    warning("These requested CD8 subtypes were not found in ", fine_celltype_col, ": ",
            paste(missing_cd8_subtypes, collapse = ", "))
  }
} else {
  warning("Fine cell-type column ", fine_celltype_col,
          " was not found. Broad cell-type analyses will still run.")
}

analysis_manifest <- rbindlist(lapply(analysis_sets, function(x) {
  data.table(
    annotation_level = x$annotation_level,
    celltype_col = x$celltype_col,
    celltype = x$celltype
  )
}), fill = TRUE)

save_csv_overwrite(as.data.frame(analysis_manifest), "03_6_analysis_celltypes")

# =========================================================
# 08 helper functions
# =========================================================
row_means_safe <- function(m) {
  if (ncol(m) == 1) as.numeric(m[, 1]) else Matrix::rowMeans(m)
}

build_fov_pseudobulk <- function(meta_sub, expr_mat, min_cells) {
  fov_counts <- meta_sub[, .N, by = .(
    FOV,
    Run_Tissue_name,
    immunotype = get(immunotype_col)
  )]
  
  valid_fov <- fov_counts[N >= min_cells, FOV]
  meta_use <- meta_sub[FOV %in% valid_fov]
  fov_counts <- fov_counts[FOV %in% valid_fov]
  
  if (nrow(meta_use) == 0) return(NULL)
  
  fov_check <- meta_use[, .(
    n_slide = uniqueN(Run_Tissue_name),
    n_immunotype = uniqueN(get(immunotype_col))
  ), by = FOV]
  
  if (any(fov_check$n_slide != 1) || any(fov_check$n_immunotype != 1)) {
    stop("A FOV maps to more than one slide or immunotype.")
  }
  
  fov_ids <- sort(unique(meta_use$FOV))
  pb_list <- vector("list", length(fov_ids))
  names(pb_list) <- fov_ids
  
  for (i in seq_along(fov_ids)) {
    fv <- fov_ids[i]
    cells_use <- meta_use[FOV == fv, cell_id]
    cells_use <- cells_use[cells_use %in% colnames(expr_mat)]
    if (length(cells_use) == 0) next
    pb_list[[i]] <- row_means_safe(expr_mat[, cells_use, drop = FALSE])
  }
  
  keep <- !vapply(pb_list, is.null, logical(1))
  pb_list <- pb_list[keep]
  if (length(pb_list) == 0) return(NULL)
  
  pb_mat <- do.call(cbind, pb_list)
  rownames(pb_mat) <- rownames(expr_mat)
  colnames(pb_mat) <- names(pb_list)
  
  fov_meta <- unique(meta_use[
    FOV %in% colnames(pb_mat),
    .(FOV, Run_Tissue_name, immunotype = get(immunotype_col))
  ])
  fov_meta <- fov_meta[match(colnames(pb_mat), FOV)]
  
  if (!identical(fov_meta$FOV, colnames(pb_mat))) {
    stop("FOV pseudobulk matrix does not align with FOV metadata.")
  }
  
  list(expr = pb_mat, meta = fov_meta, counts = fov_counts)
}

run_pairwise_limma <- function(pb_mat, fov_meta, group_a, group_b) {
  use <- fov_meta$immunotype %in% c(group_a, group_b)
  meta2 <- copy(fov_meta[use])
  mat2 <- pb_mat[, use, drop = FALSE]
  meta2 <- meta2[match(colnames(mat2), FOV)]
  
  n_a <- sum(meta2$immunotype == group_a)
  n_b <- sum(meta2$immunotype == group_b)
  
  if (n_a < min_fov_per_group || n_b < min_fov_per_group) {
    return(list(result = NULL, n_a = n_a, n_b = n_b, status = "skip_insufficient_FOVs"))
  }
  
  meta2[, immunotype_pair := factor(immunotype, levels = c(group_a, group_b))]
  
  if (uniqueN(meta2$Run_Tissue_name) >= 2) {
    meta2[, Run_Tissue_name := factor(Run_Tissue_name)]
    design <- model.matrix(~ 0 + immunotype_pair + Run_Tissue_name, data = meta2)
  } else {
    design <- model.matrix(~ 0 + immunotype_pair, data = meta2)
  }
  
  col_a <- paste0("immunotype_pair", group_a)
  col_b <- paste0("immunotype_pair", group_b)
  
  qr_design <- qr(design)
  if (qr_design$rank < ncol(design)) {
    keep_cols <- qr_design$pivot[seq_len(qr_design$rank)]
    design <- design[, keep_cols, drop = FALSE]
  }
  
  if (!all(c(col_a, col_b) %in% colnames(design))) {
    return(list(result = NULL, n_a = n_a, n_b = n_b,
                status = "skip_nonestimable_group_contrast"))
  }
  
  contrast_vec <- rep(0, ncol(design))
  names(contrast_vec) <- colnames(design)
  contrast_vec[col_a] <- 1
  contrast_vec[col_b] <- -1
  
  fit <- limma::lmFit(mat2, design)
  fit2 <- limma::contrasts.fit(
    fit,
    contrasts = matrix(
      contrast_vec,
      ncol = 1,
      dimnames = list(names(contrast_vec), paste0(group_a, "_vs_", group_b))
    )
  )
  fit2 <- limma::eBayes(fit2)
  
  tt <- limma::topTable(
    fit2,
    coef = 1,
    number = Inf,
    sort.by = "none",
    adjust.method = "BH"
  )
  
  res <- as.data.table(tt, keep.rownames = "gene")
  
  idx_a <- meta2$immunotype == group_a
  idx_b <- meta2$immunotype == group_b
  
  res[, mean_group1 := Matrix::rowMeans(mat2[, idx_a, drop = FALSE])]
  res[, mean_group2 := Matrix::rowMeans(mat2[, idx_b, drop = FALSE])]
  
  setnames(res, "logFC", "difference_mean_expression")
  
  res[, `:=`(
    group1 = group_a,
    group2 = group_b,
    contrast = paste0(group_a, "_vs_", group_b),
    n_FOV_group1 = n_a,
    n_FOV_group2 = n_b
  )]
  
  res[, direction := fifelse(
    difference_mean_expression > 0,
    paste0(group_a, " > ", group_b),
    fifelse(difference_mean_expression < 0, paste0(group_b, " > ", group_a), "equal")
  )]
  
  setcolorder(
    res,
    c(
      "gene","contrast","group1","group2",
      "n_FOV_group1","n_FOV_group2",
      "mean_group1","mean_group2",
      "difference_mean_expression","direction",
      "AveExpr","t","P.Value","adj.P.Val","B"
    )
  )
  
  list(result = res, n_a = n_a, n_b = n_b, status = "ok")
}

# =========================================================
# 09 run DEG for every cell type
# =========================================================
contrast_list <- list(
  c("ID", "IE"),
  c("ID", "IN"),
  c("IE", "IN")
)

run_log <- list()
all_results <- list()
all_fov_counts <- list()

for (set_name in names(analysis_sets)) {
  info <- analysis_sets[[set_name]]
  ct_col <- info$celltype_col
  ct <- info$celltype
  annotation_level <- info$annotation_level
  tag <- paste0(annotation_level, "__", safe_tag(ct))
  
  message("Running 03_6: ", annotation_level, " / ", ct)
  
  meta_sub <- copy(dt)[as.character(get(ct_col)) == ct]
  
  if (nrow(meta_sub) == 0) {
    run_log[[length(run_log) + 1]] <- data.table(
      annotation_level = annotation_level,
      celltype_col = ct_col,
      celltype = ct,
      contrast = NA_character_,
      n_cells = 0L,
      n_FOV = 0L,
      status = "skip_no_cells"
    )
    next
  }
  
  pb <- build_fov_pseudobulk(meta_sub, expr, min_cells_per_fov)
  
  if (is.null(pb)) {
    run_log[[length(run_log) + 1]] <- data.table(
      annotation_level = annotation_level,
      celltype_col = ct_col,
      celltype = ct,
      contrast = NA_character_,
      n_cells = nrow(meta_sub),
      n_FOV = 0L,
      status = "skip_no_FOV_after_min_cells"
    )
    next
  }
  
  pb_mat <- pb$expr
  fov_meta <- pb$meta
  
  fov_count_out <- copy(pb$counts)
  fov_count_out[, `:=`(
    annotation_level = annotation_level,
    celltype = ct
  )]
  all_fov_counts[[length(all_fov_counts) + 1]] <- fov_count_out
  
  save_csv_overwrite(as.data.frame(fov_count_out),
                     paste0("03_6_", tag, "_FOV_cell_counts"))
  save_csv_overwrite(as.data.frame(fov_meta),
                     paste0("03_6_", tag, "_FOVs_used"))
  
  ct_result_list <- list()
  
  for (contrast_pair in contrast_list) {
    group_a <- contrast_pair[1]
    group_b <- contrast_pair[2]
    contrast_name <- paste0(group_a, "_vs_", group_b)
    
    fit_out <- run_pairwise_limma(pb_mat, fov_meta, group_a, group_b)
    
    run_log[[length(run_log) + 1]] <- data.table(
      annotation_level = annotation_level,
      celltype_col = ct_col,
      celltype = ct,
      contrast = contrast_name,
      n_cells = nrow(meta_sub),
      n_FOV = nrow(fov_meta),
      n_FOV_group1 = fit_out$n_a,
      n_FOV_group2 = fit_out$n_b,
      status = fit_out$status
    )
    
    if (is.null(fit_out$result)) next
    
    res <- copy(fit_out$result)
    res[, `:=`(
      annotation_level = annotation_level,
      celltype_col = ct_col,
      celltype = ct
    )]
    
    setcolorder(
      res,
      c(
        "annotation_level","celltype_col","celltype",
        "gene","contrast","group1","group2",
        "n_FOV_group1","n_FOV_group2",
        "mean_group1","mean_group2",
        "difference_mean_expression","direction",
        "AveExpr","t","P.Value","adj.P.Val","B"
      )
    )
    
    save_csv_overwrite(
      as.data.frame(res),
      paste0("03_6_", tag, "__", contrast_name, "__all_genes")
    )
    
    ct_result_list[[length(ct_result_list) + 1]] <- res
    all_results[[length(all_results) + 1]] <- res
  }
  
  if (length(ct_result_list) > 0) {
    this_ct_dt <- rbindlist(ct_result_list, fill = TRUE)
    save_csv_overwrite(
      as.data.frame(this_ct_dt),
      paste0("03_6_", tag, "__all_contrasts_all_genes")
    )
  }
}

# =========================================================
# 10 combined master tables
# =========================================================
if (length(all_results) > 0) {
  all_deg_dt <- rbindlist(all_results, fill = TRUE)
  
  save_csv_overwrite(as.data.frame(all_deg_dt),
                     "03_6_ALL_celltypes_ALL_contrasts_ALL_genes")
  save_rds_overwrite(all_deg_dt,
                     "03_6_ALL_celltypes_ALL_contrasts_ALL_genes")
  
  broad_deg_dt <- all_deg_dt[annotation_level == "broad"]
  if (nrow(broad_deg_dt) > 0) {
    save_csv_overwrite(as.data.frame(broad_deg_dt),
                       "03_6_BROAD_celltypes_ALL_contrasts_ALL_genes")
  }
  
  cd8_deg_dt <- all_deg_dt[annotation_level == "fine_CD8"]
  if (nrow(cd8_deg_dt) > 0) {
    save_csv_overwrite(as.data.frame(cd8_deg_dt),
                       "03_6_CD8_subtypes_ALL_contrasts_ALL_genes")
  }
  
  for (contrast_name in c("ID_vs_IE", "ID_vs_IN", "IE_vs_IN")) {
    tmp <- all_deg_dt[contrast == contrast_name]
    if (nrow(tmp) > 0) {
      save_csv_overwrite(
        as.data.frame(tmp),
        paste0("03_6_ALL_celltypes__", contrast_name, "__ALL_genes")
      )
    }
  }
}

# =========================================================
# 11 run log and FOV-count summary
# =========================================================
if (length(run_log) > 0) {
  run_log_dt <- rbindlist(run_log, fill = TRUE)
  save_csv_overwrite(as.data.frame(run_log_dt), "03_6_run_log")
}

if (length(all_fov_counts) > 0) {
  fov_count_master <- rbindlist(all_fov_counts, fill = TRUE)
  save_csv_overwrite(as.data.frame(fov_count_master),
                     "03_6_ALL_celltypes_FOV_cell_counts")
}

message("03_6 differential-expression analysis finished successfully.")


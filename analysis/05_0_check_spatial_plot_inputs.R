rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# =========================
# paths
# =========================
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
raw_dir     <- file.path(project_dir, "raw")
res_dir     <- file.path(project_dir, "results", "05_0")
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

# =========================
# input files
# =========================
raw_files <- c(
  metadata = file.path(raw_dir, "TNBC_6k_metadata_Bogang_selected.RDS"),
  norm     = file.path(raw_dir, "TNBC_6k_norm_Bogang.RDS"),
  um       = file.path(raw_dir, "TNBC_6k_um_Bogang.RDS"),
  xy       = file.path(raw_dir, "TNBC_6k_xy_Bogang.RDS")
)

# =========================
# helper functions
# =========================
safe_dim <- function(x) {
  d <- dim(x)
  if (is.null(d)) return(NA_character_)
  paste(d, collapse = " x ")
}

safe_nrow <- function(x) {
  n <- tryCatch(nrow(x), error = function(e) NA_integer_)
  if (is.null(n)) NA_integer_ else n
}

safe_ncol <- function(x) {
  n <- tryCatch(ncol(x), error = function(e) NA_integer_)
  if (is.null(n)) NA_integer_ else n
}

safe_colnames <- function(x) {
  cn <- tryCatch(colnames(x), error = function(e) NULL)
  if (is.null(cn)) character(0) else cn
}

safe_rownames <- function(x) {
  rn <- tryCatch(rownames(x), error = function(e) NULL)
  if (is.null(rn)) character(0) else rn
}

capture_str_head <- function(x, max.level = 2, vec.len = 6) {
  out <- capture.output(str(x, max.level = max.level, vec.len = vec.len, list.len = 12))
  paste(out, collapse = "\n")
}

find_keywords_in_names <- function(nms, keywords) {
  if (length(nms) == 0) return(character(0))
  hit <- unique(unlist(lapply(keywords, function(k) {
    nms[grepl(k, nms, ignore.case = TRUE)]
  })))
  hit
}

inspect_object <- function(obj, obj_name = "object", parent = NA_character_, level = 0, max_depth = 2) {
  res <- list()
  
  this_names <- tryCatch(names(obj), error = function(e) NULL)
  if (is.null(this_names)) this_names <- character(0)
  
  cn <- safe_colnames(obj)
  rn <- safe_rownames(obj)
  
  all_names <- unique(c(this_names, cn, rn))
  
  key_hits <- list(
    fov       = find_keywords_in_names(all_names, c("^FOV$", "fov", "field")),
    xy        = find_keywords_in_names(all_names, c("^x$", "^y$", "centroid", "coord", "spatial")),
    cell      = find_keywords_in_names(all_names, c("cell", "cell_id", "barcode", "ID$")),
    boundary  = find_keywords_in_names(all_names, c("bound", "polygon", "vertex", "seg", "outline", "mask")),
    molecule  = find_keywords_in_names(all_names, c("molecule", "transcript", "gene", "feature", "rna"))
  )
  
  res[[1]] <- data.table(
    parent_name = parent,
    object_name = obj_name,
    class = paste(class(obj), collapse = "; "),
    typeof = typeof(obj),
    dim = safe_dim(obj),
    nrow = safe_nrow(obj),
    ncol = safe_ncol(obj),
    length = tryCatch(length(obj), error = function(e) NA_integer_),
    names_n = length(this_names),
    colnames_n = length(cn),
    rownames_n = length(rn),
    hit_fov = paste(key_hits$fov, collapse = "; "),
    hit_xy = paste(key_hits$xy, collapse = "; "),
    hit_cell = paste(key_hits$cell, collapse = "; "),
    hit_boundary = paste(key_hits$boundary, collapse = "; "),
    hit_molecule = paste(key_hits$molecule, collapse = "; ")
  )
  
  if (level < max_depth && is.list(obj) && length(obj) > 0) {
    sub_n <- min(length(obj), 30L)
    for (i in seq_len(sub_n)) {
      nm <- this_names[i]
      if (is.null(nm) || is.na(nm) || nm == "") nm <- paste0("[[", i, "]]")
      child <- tryCatch(obj[[i]], error = function(e) NULL)
      if (is.null(child)) next
      res[[length(res) + 1]] <- inspect_object(
        obj = child,
        obj_name = nm,
        parent = obj_name,
        level = level + 1,
        max_depth = max_depth
      )
    }
  }
  
  rbindlist(res, fill = TRUE)
}

# =========================
# load objects
# =========================
obj_list <- lapply(raw_files, readRDS)

# =========================
# inspect top-level
# =========================
report_lines <- c()
inspect_tables <- list()

for (nm in names(obj_list)) {
  obj <- obj_list[[nm]]
  
  report_lines <- c(
    report_lines,
    paste0("========== ", nm, " =========="),
    paste0("file: ", raw_files[[nm]]),
    paste0("class: ", paste(class(obj), collapse = "; ")),
    paste0("typeof: ", typeof(obj)),
    paste0("dim: ", safe_dim(obj)),
    "---- str() head ----",
    capture_str_head(obj, max.level = 2, vec.len = 6),
    ""
  )
  
  inspect_tables[[nm]] <- inspect_object(obj, obj_name = nm, parent = NA_character_, level = 0, max_depth = 2)
}

inspect_dt <- rbindlist(inspect_tables, fill = TRUE, idcol = "source_object")

fwrite(inspect_dt, file.path(res_dir, "05_0_object_structure_summary.csv"))
writeLines(report_lines, con = file.path(res_dir, "05_0_object_structure_report.txt"))

# =========================
# metadata + xy merge check
# =========================
meta <- as.data.table(obj_list$metadata)

xy_obj <- obj_list$xy
um_obj <- obj_list$um

xy <- as.data.table(xy_obj)
um <- as.data.table(um_obj)

# In the project RDS files, cell IDs are stored as matrix row names.
if (!"cell_id" %in% names(xy) && !is.null(rownames(xy_obj))) {
  xy[, cell_id := rownames(xy_obj)]
}
if (!"cell_id" %in% names(um) && !is.null(rownames(um_obj))) {
  um[, cell_id := rownames(um_obj)]
}

# Standardize spatial coordinate names used below.
if (all(c("x_mm", "y_mm") %in% names(xy))) {
  setnames(xy, c("x_mm", "y_mm"), c("x", "y"))
}

fwrite(data.table(
  object = c("metadata", "xy", "um"),
  nrow = c(nrow(meta), nrow(xy), nrow(um)),
  ncol = c(ncol(meta), ncol(xy), ncol(um))
), file.path(res_dir, "05_0_basic_dimensions.csv"))

writeLines(colnames(meta), con = file.path(res_dir, "05_0_metadata_colnames.txt"))
writeLines(colnames(xy),   con = file.path(res_dir, "05_0_xy_colnames.txt"))
writeLines(colnames(um),   con = file.path(res_dir, "05_0_um_colnames.txt"))

# 尝试按 cell_id 合并
merge_ok <- FALSE
dt <- NULL

if ("cell_id" %in% colnames(meta) && "cell_id" %in% colnames(xy)) {
  dt <- merge(meta, xy, by = "cell_id", all.x = TRUE)
  merge_ok <- TRUE
} else if ("cell_id" %in% colnames(meta) && "cell_id" %in% colnames(um)) {
  dt <- merge(meta, um, by = "cell_id", all.x = TRUE)
  merge_ok <- TRUE
}

if (merge_ok) {
  fwrite(data.table(
    nrow = nrow(dt),
    ncol = ncol(dt),
    miss_x = sum(!("x" %in% colnames(dt)) | is.na(dt$x)),
    miss_y = sum(!("y" %in% colnames(dt)) | is.na(dt$y)),
    n_FOV = if ("FOV" %in% colnames(dt)) uniqueN(dt$FOV) else NA_integer_,
    n_Run_Tissue_name = if ("Run_Tissue_name" %in% colnames(dt)) uniqueN(dt$Run_Tissue_name) else NA_integer_
  ), file.path(res_dir, "05_0_merged_dt_summary.csv"))
  
  fwrite(head(dt, 50), file.path(res_dir, "05_0_merged_dt_head50.csv"))
}

# =========================
# check if centroid plot is possible
# =========================
plot_flag <- FALSE
plot_msg <- c()

if (merge_ok) {
  needed_cols <- c("FOV", "x", "y")
  miss_cols <- setdiff(needed_cols, colnames(dt))
  
  if (length(miss_cols) == 0) {
    dt_plot <- copy(dt)[!is.na(FOV) & !is.na(x) & !is.na(y)]
    if (nrow(dt_plot) > 0) {
      plot_flag <- TRUE
      
      # 选一个细胞数最多的 FOV 作为示例
      fov_stat <- dt_plot[, .N, by = FOV][order(-N)]
      fwrite(fov_stat, file.path(res_dir, "05_0_fov_cell_counts.csv"))
      
      demo_fov <- fov_stat$FOV[1]
      
      # 颜色列优先用 BC6k_refined
      color_col <- NULL
      if ("BC6k_refined" %in% colnames(dt_plot)) {
        color_col <- "BC6k_refined"
      } else if ("Cell_Type_TBcoarse_MPbase" %in% colnames(dt_plot)) {
        color_col <- "Cell_Type_TBcoarse_MPbase"
      } else if ("Cell_Type_TBcoarse_MPbranch" %in% colnames(dt_plot)) {
        color_col <- "Cell_Type_TBcoarse_MPbranch"
      }
      
      dt_demo <- dt_plot[FOV == demo_fov]
      
      if (!is.null(color_col)) {
        p <- ggplot(dt_demo, aes(x = x, y = y, color = .data[[color_col]])) +
          geom_point(size = 0.5, alpha = 0.9) +
          coord_fixed() +
          theme_bw(base_size = 10) +
          labs(
            title = paste0("Centroid spatial plot demo | FOV: ", demo_fov),
            color = color_col
          ) +
          theme(
            plot.title = element_text(size = 11),
            axis.title = element_blank()
          )
        
        ggsave(
          filename = file.path(res_dir, "05_0_demo_centroid_plot.png"),
          plot = p, width = 8, height = 7, dpi = 300
        )
      } else {
        p <- ggplot(dt_demo, aes(x = x, y = y)) +
          geom_point(size = 0.5, alpha = 0.9, color = "black") +
          coord_fixed() +
          theme_bw(base_size = 10) +
          labs(
            title = paste0("Centroid spatial plot demo | FOV: ", demo_fov)
          ) +
          theme(
            plot.title = element_text(size = 11),
            axis.title = element_blank()
          )
        
        ggsave(
          filename = file.path(res_dir, "05_0_demo_centroid_plot.png"),
          plot = p, width = 8, height = 7, dpi = 300
        )
      }
      
      plot_msg <- c(plot_msg, paste0("Centroid plot demo generated for FOV: ", demo_fov))
    } else {
      plot_msg <- c(plot_msg, "Merged dt exists, but no non-missing FOV/x/y rows found.")
    }
  } else {
    plot_msg <- c(plot_msg, paste0("Missing required columns for centroid plot: ", paste(miss_cols, collapse = ", ")))
  }
} else {
  plot_msg <- c(plot_msg, "metadata and xy/um could not be merged by cell_id.")
}

writeLines(plot_msg, con = file.path(res_dir, "05_0_plot_check.txt"))

# =========================
# expression matrix quick check
# =========================
norm_obj <- obj_list$norm

expr_report <- c(
  paste0("class: ", paste(class(norm_obj), collapse = "; ")),
  paste0("typeof: ", typeof(norm_obj)),
  paste0("dim: ", safe_dim(norm_obj))
)

expr_colnames <- safe_colnames(norm_obj)
expr_rownames <- safe_rownames(norm_obj)

expr_report <- c(
  expr_report,
  paste0("n_colnames: ", length(expr_colnames)),
  paste0("n_rownames: ", length(expr_rownames)),
  paste0("has cell_id in colnames: ", "cell_id" %in% expr_colnames),
  paste0("example rownames: ", paste(head(expr_rownames, 10), collapse = ", ")),
  paste0("example colnames: ", paste(head(expr_colnames, 10), collapse = ", "))
)

writeLines(expr_report, con = file.path(res_dir, "05_0_expression_matrix_check.txt"))

cat("05_0 check completed.\n")
cat("Results written to:\n", res_dir, "\n")
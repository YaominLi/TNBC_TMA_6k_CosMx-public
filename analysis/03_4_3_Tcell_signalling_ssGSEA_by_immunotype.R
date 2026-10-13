rm(list = ls())
gc()
graphics.off()
options(stringsAsFactors = FALSE)
set.seed(1)

# =========================================================
# 03_4_3 v4
# T-cell signalling ssGSEA across ID / IE / IN
#
# Three score systems:
#   1) normalized_ssGSEA
#   2) pathway_robust_z
#   3) FOV_centered
#
# Three statistical frameworks:
#   A) LM_Tukey
#      score ~ immunotype + Run_Tissue_name
#      overall partial F-test for immunotype
#      pairwise estimated marginal means with Tukey adjustment
#
#   B) KW_WilcoxonBH
#      Kruskal-Wallis overall test
#      pairwise Wilcoxon rank-sum tests with BH correction across 3 pairs
#
#   C) Pairwise_Welch_rawP
#      each comparison is tested independently as a two-group comparison:
#      ID vs IE, ID vs IN, IE vs IN
#      two-sided Welch t-test; NO three-group omnibus test is used
#      raw pairwise p values are displayed on the figure
#      BH-adjusted values are also saved in CSV for reference
#
# Biological replicate:
#   one patient = one FOV
#
# Populations:
#   Total_T = T.cell.CD4 + T.cell.CD8 + T.cell.regulatory
#   CD8_T   = T.cell.CD8
#
# Figures:
#   9 violin figure sets = 3 scores x 3 statistical frameworks.
#   Each PDF contains Total_T page(s) followed by CD8_T page(s).
#   LM/Tukey and KW/Wilcoxon figures show overall p + all 3 pairwise p values.
#   Pairwise_Welch_rawP figures show only the 3 independently calculated
#   two-group p values (no overall three-group p).
#
# Mean heatmaps:
#   restored for all 3 score systems and both populations.
#   Columns = ID / IE / IN; rows = pathways; values are group means,
#   displayed after row-wise z-scaling as in the previous 03_4_3 workflow.
#
# Immunotype display/colors are fixed:
#   ID = pure black #000000
#   IE = pure blue  #0000FF
#   IN = pure red   #FF0000
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
if (
  is.na(proj) ||
  !file.exists(file.path(proj, "raw")) ||
  !file.exists(file.path(proj, "results", "03_2", "data"))
) {
  proj <- find_project_root()
}

raw_dir <- file.path(proj, "raw")
in_dir  <- file.path(proj, "results", "03_2", "data")
out_dir <- file.path(proj, "results", "03_4_3")
fig_dir <- file.path(out_dir, "figures")
dat_dir <- file.path(out_dir, "data")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dat_dir, recursive = TRUE, showWarnings = FALSE)

meta_file <- file.path(in_dir, "03_2_meta_with_immunotype_BW_top_40.rds")
norm_file <- file.path(raw_dir, "TNBC_6k_norm_Bogang.RDS")
xlsx_path <- file.path(raw_dir, "T_cell_signalling.xlsx")

# =========================================================
# 02 packages
# =========================================================
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(readxl)
  library(GSVA)
  library(BiocParallel)
  library(SummarizedExperiment)
  library(matrixStats)
  library(ggplot2)
  library(emmeans)
})

# Windows: keep GSVA serial to avoid BiocParallel port/worker instability.
BiocParallel::register(BiocParallel::SerialParam(), default = TRUE)

# =========================================================
# 03 user parameters
# =========================================================

# Fixed order: ID left, IE middle, IN right.
immunotype_order <- c("ID", "IE", "IN")
immunotype_cols <- c(
  "ID" = "#000000",
  "IE" = "#0000FF",
  "IN" = "#FF0000"
)

# Current broad T-cell annotations.
total_T_labels <- c(
  "T.cell.CD4",
  "T.cell.CD8",
  "T.cell.regulatory"
)
CD8_label <- "T.cell.CD8"

# Minimum selected cells per FOV before FOV pseudobulk.
min_cells_per_FOV_total_T <- 10L
min_cells_per_FOV_CD8 <- 5L

# Minimum FOVs per immunotype required for statistical tests.
min_FOV_per_group_for_stats <- 3L

# ssGSEA settings.
min_genes_per_pathway <- 3L
kcdf_mode <- "Gaussian"

# Plot layout: intentionally taller than previous versions.
violin_per_page <- 18L
violin_ncol <- 3L
violin_pdf_width <- 11
violin_pdf_height <- 15.5
jitter_width <- 0.14
jitter_size <- 0.95

# Mean-heatmap dimensions.
mean_heatmap_width <- 5.8
mean_heatmap_min_height <- 4.8
mean_heatmap_height_per_pathway <- 0.30

# Whether to label nonsignificant pairwise comparisons as "ns".
show_ns <- TRUE

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
  invisible(out)
}

save_pdf <- function(p, filename, w = 8, h = 6) {
  out <- file.path(fig_dir, paste0(filename, ".pdf"))
  if (file.exists(out)) file.remove(out)
  
  ggplot2::ggsave(
    filename = out,
    plot = p,
    width = w,
    height = h,
    device = grDevices::cairo_pdf
  )
  
  message("Saved: ", out)
  invisible(out)
}

safe_tag <- function(x) {
  gsub("[^A-Za-z0-9_]+", "_", x)
}

format_p <- function(p) {
  if (length(p) == 0 || is.na(p) || !is.finite(p)) return("NA")
  if (p < 1e-4) return(formatC(p, format = "e", digits = 1))
  if (p < 0.001) return(formatC(p, format = "f", digits = 4))
  formatC(p, format = "f", digits = 3)
}

p_to_star <- function(p) {
  if (length(p) == 0 || is.na(p) || !is.finite(p)) return("")
  if (p < 1e-4) return("****")
  if (p < 1e-3) return("***")
  if (p < 1e-2) return("**")
  if (p < 5e-2) return("*")
  if (show_ns) return("ns")
  ""
}

# =========================================================
# 05 input checks
# =========================================================
if (!file.exists(meta_file)) stop("Missing metadata file: ", meta_file)
if (!file.exists(norm_file)) stop("Missing expression file: ", norm_file)
if (!file.exists(xlsx_path)) stop("Missing pathway file: ", xlsx_path)

# =========================================================
# 06 read T-cell signalling pathways
# Same wide-table format as the previous scRNA workflow:
# columns = pathways; cells below each column = HGNC symbols.
# =========================================================
read_pathways_wide_xlsx_sym <- function(xlsx_path, sheet = 1, min_genes = 3) {
  tb <- readxl::read_excel(xlsx_path, sheet = sheet)
  
  nm <- toupper(trimws(colnames(tb)))
  nm <- make.unique(nm, sep = "_DUP")
  colnames(tb) <- nm
  
  lst <- lapply(tb, function(col) {
    x <- toupper(trimws(as.character(col)))
    x <- x[!is.na(x) & nzchar(x)]
    unique(x)
  })
  
  lens <- vapply(lst, length, integer(1))
  lst[lens >= min_genes]
}

pathways_input <- read_pathways_wide_xlsx_sym(
  xlsx_path = xlsx_path,
  min_genes = min_genes_per_pathway
)

if (length(pathways_input) == 0) {
  stop("No pathway in T_cell_signalling.xlsx contains enough genes.")
}

# =========================================================
# 07 GSVA/ssGSEA compatibility wrapper
# Follows the GSVA 2.x Param-interface logic used in the prior scRNA script.
# =========================================================
run_ssgsea <- function(expr_mat, gene_sets, min_size = 3, kcdf = "Gaussian") {
  build_args <- function(fun) {
    fn <- get(fun, asNamespace("GSVA"))
    a <- names(formals(fn))
    args <- list()
    
    if ("exprData" %in% a) args$exprData <- expr_mat
    if ("geneSets" %in% a) args$geneSets <- gene_sets
    if ("minSize" %in% a) args$minSize <- min_size
    if ("maxSize" %in% a) args$maxSize <- Inf
    if ("kcdf" %in% a) args$kcdf <- kcdf
    if ("norm" %in% a) args$norm <- TRUE
    if ("ssgseaNorm" %in% a) args$ssgseaNorm <- TRUE
    if ("maxDiff" %in% a) args$maxDiff <- TRUE
    if ("mx.diff" %in% a) args$mx.diff <- TRUE
    
    args
  }
  
  param <- NULL
  
  args1 <- build_args("ssgseaParam")
  if (length(args1) > 0) {
    param <- try(do.call(GSVA::ssgseaParam, args1), silent = TRUE)
  }
  
  if (inherits(param, "try-error") || is.null(param)) {
    args2 <- build_args("gsvaParam")
    if ("method" %in% names(formals(get("gsvaParam", asNamespace("GSVA"))))) {
      args2$method <- "ssgsea"
    }
    param <- try(do.call(GSVA::gsvaParam, args2), silent = TRUE)
  }
  
  if (inherits(param, "try-error") || is.null(param)) {
    stop("Unable to construct GSVA ssGSEA parameter object.")
  }
  
  es <- GSVA::gsva(param = param, verbose = FALSE)
  
  if (inherits(es, "SummarizedExperiment")) {
    return(SummarizedExperiment::assay(es))
  }
  
  if (is.matrix(es)) return(es)
  
  stop("Unknown GSVA output class: ", paste(class(es), collapse = ", "))
}

# =========================================================
# 08 load metadata and expression
# =========================================================
meta <- as.data.table(readRDS(meta_file))
norm <- readRDS(norm_file)

need_meta <- c(
  "cell_id",
  "FOV",
  "Run_Tissue_name",
  "BC6k_refined",
  "immunotype_BW_top_40"
)

miss_meta <- setdiff(need_meta, names(meta))
if (length(miss_meta) > 0) {
  stop("Missing metadata columns: ", paste(miss_meta, collapse = ", "))
}

meta <- meta[
  !is.na(cell_id) &
    !is.na(FOV) &
    !is.na(Run_Tissue_name) &
    !is.na(BC6k_refined) &
    !is.na(immunotype_BW_top_40)
]

meta[, cell_id := as.character(cell_id)]
meta[, FOV := as.character(FOV)]
meta[, Run_Tissue_name := as.character(Run_Tissue_name)]
meta[, BC6k_refined := as.character(BC6k_refined)]
meta[, immunotype_BW_top_40 := as.character(immunotype_BW_top_40)]
meta <- meta[immunotype_BW_top_40 %in% immunotype_order]

get_expr <- function(norm) {
  if (is.matrix(norm) || inherits(norm, "Matrix")) return(norm)
  
  if (inherits(norm, "Seurat")) {
    assay <- if ("RNA" %in% names(norm@assays)) {
      "RNA"
    } else {
      Seurat::DefaultAssay(norm)
    }
    
    return(
      tryCatch(
        Seurat::GetAssayData(norm, assay = assay, slot = "data"),
        error = function(e) {
          Seurat::GetAssayData(norm, assay = assay, slot = "counts")
        }
      )
    )
  }
  
  if (is.list(norm)) {
    for (nm in c("data", "expr", "expression", "mat", "matrix", "norm", "logcounts", "counts")) {
      if (
        !is.null(norm[[nm]]) &&
        (is.matrix(norm[[nm]]) || inherits(norm[[nm]], "Matrix"))
      ) {
        return(norm[[nm]])
      }
    }
  }
  
  stop("Cannot extract expression matrix from norm.")
}

expr <- get_expr(norm)

# Current normalized expression is expected as genes x cells.
# Transpose only if cell IDs are found on rows rather than columns.
cell_hit_col <- sum(meta$cell_id %in% colnames(expr))
cell_hit_row <- sum(meta$cell_id %in% rownames(expr))

if (cell_hit_col == 0 && cell_hit_row > 0) {
  expr <- t(expr)
}

if (is.null(rownames(expr)) || is.null(colnames(expr))) {
  stop("Expression matrix requires both row names and column names.")
}

common_cells <- intersect(meta$cell_id, colnames(expr))
common_cells <- common_cells[match(common_cells, meta$cell_id)]
common_cells <- common_cells[!is.na(common_cells)]

if (length(common_cells) == 0) {
  stop("No shared cells between metadata and normalized expression.")
}

expr <- expr[, common_cells, drop = FALSE]
meta <- meta[match(common_cells, cell_id)]
stopifnot(identical(colnames(expr), meta$cell_id))

# =========================================================
# 09 pathway coverage in the CosMx 6k panel
# =========================================================
expr_gene_upper <- toupper(rownames(expr))

pathway_qc <- rbindlist(
  lapply(names(pathways_input), function(pw) {
    genes_in <- pathways_input[[pw]]
    genes_hit <- intersect(genes_in, unique(expr_gene_upper))
    
    data.table(
      pathway = pw,
      n_genes_input = length(genes_in),
      n_genes_detected_6k = length(genes_hit),
      usable = length(genes_hit) >= min_genes_per_pathway,
      genes_detected_6k = paste(genes_hit, collapse = ";"),
      genes_missing_6k = paste(setdiff(genes_in, genes_hit), collapse = ";")
    )
  })
)

save_dat(
  pathway_qc,
  "03_4_3_Tcell_signalling_pathway_coverage_6k",
  "csv"
)

pathways_use <- pathways_input[pathway_qc[usable == TRUE, pathway]]

if (length(pathways_use) == 0) {
  stop("No T-cell signalling pathway retains >=3 genes in the CosMx 6k panel.")
}

# =========================================================
# 10 FOV pseudobulk
# =========================================================
collapse_uppercase_gene_rows <- function(mat) {
  gene_upper <- toupper(rownames(mat))
  
  if (!anyDuplicated(gene_upper)) {
    rownames(mat) <- gene_upper
    return(mat)
  }
  
  mat_dense <- as.matrix(mat)
  sums <- rowsum(mat_dense, group = gene_upper, reorder = FALSE)
  n_dup <- table(gene_upper)
  sweep(sums, 1, as.numeric(n_dup[rownames(sums)]), "/")
}

build_FOV_pseudobulk <- function(meta_subset, min_cells_per_FOV, population_name) {
  fov_count <- meta_subset[
    ,
    .(n_cells = .N),
    by = .(FOV, Run_Tissue_name, immunotype_BW_top_40)
  ]
  
  fov_count[, passed := n_cells >= min_cells_per_FOV]
  
  save_dat(
    fov_count,
    paste0("03_4_3_", safe_tag(population_name), "_cell_counts_per_FOV"),
    "csv"
  )
  
  use_dt <- meta_subset[FOV %in% fov_count[passed == TRUE, FOV]]
  
  if (nrow(use_dt) == 0) {
    stop("No FOV passed the cell-count threshold for ", population_name)
  }
  
  expr_sel <- expr[, use_dt$cell_id, drop = FALSE]
  stopifnot(identical(colnames(expr_sel), use_dt$cell_id))
  
  fov_anno <- unique(
    use_dt[, .(FOV, Run_Tissue_name, immunotype_BW_top_40)]
  )
  
  conflict <- fov_anno[, .N, by = FOV][N > 1]
  if (nrow(conflict) > 0) {
    stop("Conflicting FOV annotations detected for ", population_name)
  }
  
  fov_list <- split(use_dt$cell_id, use_dt$FOV)
  
  pb_list <- lapply(names(fov_list), function(fv) {
    m <- expr_sel[, fov_list[[fv]], drop = FALSE]
    
    if (ncol(m) == 1) {
      as.numeric(m[, 1])
    } else {
      Matrix::rowMeans(m)
    }
  })
  
  pb_mat <- do.call(cbind, pb_list)
  rownames(pb_mat) <- rownames(expr_sel)
  colnames(pb_mat) <- names(fov_list)
  pb_mat <- collapse_uppercase_gene_rows(pb_mat)
  
  fov_anno <- fov_anno[match(colnames(pb_mat), FOV)]
  stopifnot(identical(fov_anno$FOV, colnames(pb_mat)))
  
  fov_n_cells <- fov_count[passed == TRUE, .(FOV, n_cells)]
  fov_anno <- merge(fov_anno, fov_n_cells, by = "FOV", all.x = TRUE, sort = FALSE)
  fov_anno <- fov_anno[match(colnames(pb_mat), FOV)]
  
  list(
    pb_mat = pb_mat,
    fov_anno = fov_anno,
    fov_count = fov_count
  )
}

# =========================================================
# 11 calculate the three score systems
# =========================================================
score_FOV_ssGSEA <- function(pb_mat) {
  pathways_this <- lapply(
    pathways_use,
    function(gs) intersect(gs, rownames(pb_mat))
  )
  
  pathways_this <- pathways_this[
    vapply(pathways_this, length, integer(1)) >= min_genes_per_pathway
  ]
  
  if (length(pathways_this) == 0) {
    stop("No pathway retains enough genes after FOV pseudobulk preparation.")
  }
  
  # Score 1: normalized ssGSEA score returned by GSVA/ssGSEA.
  normalized_score <- run_ssgsea(
    expr_mat = as.matrix(pb_mat),
    gene_sets = pathways_this,
    min_size = min_genes_per_pathway,
    kcdf = kcdf_mode
  )
  
  # Score 2: pathway-wise robust z-score.
  # For each pathway independently, standardize across all FOVs:
  #   (score - pathway median across FOVs) / pathway MAD across FOVs.
  pathway_robust_z <- normalized_score
  
  for (i in seq_len(nrow(normalized_score))) {
    v <- as.numeric(normalized_score[i, ])
    med_i <- median(v, na.rm = TRUE)
    mad_i <- mad(v, center = med_i, constant = 1.4826, na.rm = TRUE)
    
    # Rare fallback if MAD is zero.
    if (!is.finite(mad_i) || mad_i <= 0) {
      mad_i <- sd(v, na.rm = TRUE)
    }
    if (!is.finite(mad_i) || mad_i <= 0) {
      mad_i <- 1
    }
    
    pathway_robust_z[i, ] <- (v - med_i) / mad_i
  }
  
  # Score 3: FOV-centered ssGSEA.
  # Within each FOV, center the pathway-score distribution using median/MAD.
  fov_med <- matrixStats::colMedians(normalized_score, na.rm = TRUE)
  fov_mad <- matrixStats::colMads(normalized_score, na.rm = TRUE)
  fov_mad[!is.finite(fov_mad) | fov_mad <= 0] <- 1
  
  FOV_centered <- sweep(normalized_score, 2, fov_med, "-")
  FOV_centered <- sweep(FOV_centered, 2, fov_mad, "/")
  
  list(
    normalized_ssGSEA = normalized_score,
    pathway_robust_z = pathway_robust_z,
    FOV_centered = FOV_centered,
    pathways_used = pathways_this
  )
}

score_matrix_to_long <- function(score_mat, fov_anno, population_name, score_mode) {
  d <- as.data.table(as.table(score_mat))
  setnames(d, c("pathway", "FOV", "score"))
  
  d[, pathway := as.character(pathway)]
  d[, FOV := as.character(FOV)]
  d[, score := as.numeric(score)]
  
  d <- merge(d, fov_anno, by = "FOV", all.x = TRUE, sort = FALSE)
  
  d[, population := population_name]
  d[, score_mode := score_mode]
  d[, immunotype_BW_top_40 := factor(immunotype_BW_top_40, levels = immunotype_order)]
  
  d
}

# =========================================================
# 12 statistical framework A: linear model + Tukey
# one patient = one FOV, so FOVs are independent biological replicates.
# Run_Tissue_name is treated as a fixed technical/batch covariate when estimable.
# =========================================================
run_stats_LM_Tukey <- function(score_long, population_name, score_mode_name) {
  group_count <- unique(
    score_long[
      ,
      .(n_FOV = uniqueN(FOV)),
      by = immunotype_BW_top_40
    ]
  )
  
  group_count <- merge(
    data.table(
      immunotype_BW_top_40 = factor(immunotype_order, levels = immunotype_order)
    ),
    group_count,
    by = "immunotype_BW_top_40",
    all.x = TRUE
  )
  group_count[is.na(n_FOV), n_FOV := 0L]
  
  if (!all(group_count$n_FOV >= min_FOV_per_group_for_stats)) {
    warning(
      population_name, " / ", score_mode_name,
      ": insufficient FOVs in one or more immunotypes for LM/Tukey statistics."
    )
    return(list(overall = data.table(), pairwise = data.table(), group_count = group_count))
  }
  
  pathways <- unique(score_long$pathway)
  overall_list <- list()
  pairwise_list <- list()
  
  for (pw in pathways) {
    d <- copy(score_long[pathway == pw & is.finite(score)])
    if (nrow(d) == 0 || uniqueN(d$immunotype_BW_top_40) < 3) next
    
    d[, immunotype_BW_top_40 := factor(immunotype_BW_top_40, levels = immunotype_order)]
    d[, Run_Tissue_name := factor(Run_Tissue_name)]
    
    # Prefer the slide-adjusted fixed-effect model.
    fit_full <- tryCatch(
      lm(score ~ immunotype_BW_top_40 + Run_Tissue_name, data = d),
      error = function(e) NULL
    )
    
    model_used <- "score ~ immunotype + Run_Tissue_name"
    slide_adjusted <- TRUE
    
    # If slide causes rank deficiency/confounding, fall back to one-way LM.
    if (
      is.null(fit_full) ||
      qr(model.matrix(fit_full))$rank < ncol(model.matrix(fit_full))
    ) {
      fit_full <- lm(score ~ immunotype_BW_top_40, data = d)
      model_used <- "score ~ immunotype"
      slide_adjusted <- FALSE
    }
    
    # Overall partial F-test for the immunotype term.
    if (slide_adjusted) {
      fit_reduced <- lm(score ~ Run_Tissue_name, data = d)
      av <- anova(fit_reduced, fit_full)
      overall_p <- as.numeric(av$`Pr(>F)`[2])
    } else {
      av <- anova(fit_full)
      overall_p <- as.numeric(av$`Pr(>F)`[1])
    }
    
    sm <- d[
      ,
      .(
        n_FOV = uniqueN(FOV),
        mean_score = mean(score, na.rm = TRUE),
        median_score = median(score, na.rm = TRUE),
        sd_score = sd(score, na.rm = TRUE)
      ),
      by = immunotype_BW_top_40
    ]
    
    get_group_value <- function(grp, col) {
      z <- sm[immunotype_BW_top_40 == grp, get(col)]
      if (length(z) == 0) return(NA_real_)
      as.numeric(z[1])
    }
    
    overall_list[[length(overall_list) + 1]] <- data.table(
      population = population_name,
      score_mode = score_mode_name,
      stat_framework = "LM_Tukey",
      pathway = pw,
      model_used = model_used,
      slide_adjusted = slide_adjusted,
      overall_p = overall_p,
      n_ID = get_group_value("ID", "n_FOV"),
      n_IE = get_group_value("IE", "n_FOV"),
      n_IN = get_group_value("IN", "n_FOV"),
      mean_ID = get_group_value("ID", "mean_score"),
      mean_IE = get_group_value("IE", "mean_score"),
      mean_IN = get_group_value("IN", "mean_score"),
      median_ID = get_group_value("ID", "median_score"),
      median_IE = get_group_value("IE", "median_score"),
      median_IN = get_group_value("IN", "median_score")
    )
    
    em <- tryCatch(
      emmeans::emmeans(fit_full, specs = ~ immunotype_BW_top_40),
      error = function(e) NULL
    )
    
    if (!is.null(em)) {
      pr <- tryCatch(
        as.data.table(
          summary(
            emmeans::contrast(
              em,
              method = "pairwise",
              adjust = "tukey"
            )
          )
        ),
        error = function(e) data.table()
      )
      
      if (nrow(pr) > 0) {
        pr[, contrast_std := gsub(" - ", "_vs_", contrast, fixed = TRUE)]
        pr[, population := population_name]
        pr[, score_mode := score_mode_name]
        pr[, stat_framework := "LM_Tukey"]
        pr[, pathway := pw]
        pr[, model_used := model_used]
        pr[, slide_adjusted := slide_adjusted]
        pr[, p_display := p.value]
        
        pairwise_list[[length(pairwise_list) + 1]] <- pr
      }
    }
  }
  
  overall_dt <- if (length(overall_list) > 0) {
    rbindlist(overall_list, fill = TRUE)
  } else {
    data.table()
  }
  
  pairwise_dt <- if (length(pairwise_list) > 0) {
    rbindlist(pairwise_list, fill = TRUE)
  } else {
    data.table()
  }
  
  if (nrow(overall_dt) > 0) {
    overall_dt[, overall_FDR_across_pathways := p.adjust(overall_p, method = "BH")]
  }
  
  list(
    overall = overall_dt,
    pairwise = pairwise_dt,
    group_count = group_count
  )
}

# =========================================================
# 13 statistical framework B: Kruskal-Wallis + pairwise Wilcoxon BH
# This is a nonparametric sensitivity analysis.
# It does not include Run_Tissue_name as a covariate.
# =========================================================
run_stats_KW_WilcoxonBH <- function(score_long, population_name, score_mode_name) {
  group_count <- unique(
    score_long[
      ,
      .(n_FOV = uniqueN(FOV)),
      by = immunotype_BW_top_40
    ]
  )
  
  group_count <- merge(
    data.table(
      immunotype_BW_top_40 = factor(immunotype_order, levels = immunotype_order)
    ),
    group_count,
    by = "immunotype_BW_top_40",
    all.x = TRUE
  )
  group_count[is.na(n_FOV), n_FOV := 0L]
  
  if (!all(group_count$n_FOV >= min_FOV_per_group_for_stats)) {
    warning(
      population_name, " / ", score_mode_name,
      ": insufficient FOVs in one or more immunotypes for KW/Wilcoxon statistics."
    )
    return(list(overall = data.table(), pairwise = data.table(), group_count = group_count))
  }
  
  pathways <- unique(score_long$pathway)
  overall_list <- list()
  pairwise_list <- list()
  
  pair_defs <- data.table(
    group1 = c("ID", "ID", "IE"),
    group2 = c("IE", "IN", "IN")
  )
  
  for (pw in pathways) {
    d <- copy(score_long[pathway == pw & is.finite(score)])
    if (nrow(d) == 0 || uniqueN(d$immunotype_BW_top_40) < 3) next
    
    d[, immunotype_BW_top_40 := factor(immunotype_BW_top_40, levels = immunotype_order)]
    
    kw <- tryCatch(
      kruskal.test(score ~ immunotype_BW_top_40, data = d),
      error = function(e) NULL
    )
    
    overall_p <- if (is.null(kw)) NA_real_ else as.numeric(kw$p.value)
    
    sm <- d[
      ,
      .(
        n_FOV = uniqueN(FOV),
        mean_score = mean(score, na.rm = TRUE),
        median_score = median(score, na.rm = TRUE),
        sd_score = sd(score, na.rm = TRUE)
      ),
      by = immunotype_BW_top_40
    ]
    
    get_group_value <- function(grp, col) {
      z <- sm[immunotype_BW_top_40 == grp, get(col)]
      if (length(z) == 0) return(NA_real_)
      as.numeric(z[1])
    }
    
    overall_list[[length(overall_list) + 1]] <- data.table(
      population = population_name,
      score_mode = score_mode_name,
      stat_framework = "KW_WilcoxonBH",
      pathway = pw,
      overall_p = overall_p,
      n_ID = get_group_value("ID", "n_FOV"),
      n_IE = get_group_value("IE", "n_FOV"),
      n_IN = get_group_value("IN", "n_FOV"),
      mean_ID = get_group_value("ID", "mean_score"),
      mean_IE = get_group_value("IE", "mean_score"),
      mean_IN = get_group_value("IN", "mean_score"),
      median_ID = get_group_value("ID", "median_score"),
      median_IE = get_group_value("IE", "median_score"),
      median_IN = get_group_value("IN", "median_score")
    )
    
    pr <- copy(pair_defs)
    pr[, p_raw := vapply(seq_len(.N), function(j) {
      x <- d[immunotype_BW_top_40 == group1[j], score]
      y <- d[immunotype_BW_top_40 == group2[j], score]
      
      tryCatch(
        wilcox.test(x, y, exact = FALSE)$p.value,
        error = function(e) NA_real_
      )
    }, numeric(1))]
    
    # Three pre-defined pairwise comparisons within each pathway.
    pr[, p_BH3 := p.adjust(p_raw, method = "BH")]
    pr[, p_display := p_BH3]
    pr[, contrast_std := paste0(group1, "_vs_", group2)]
    pr[, population := population_name]
    pr[, score_mode := score_mode_name]
    pr[, stat_framework := "KW_WilcoxonBH"]
    pr[, pathway := pw]
    
    pairwise_list[[length(pairwise_list) + 1]] <- pr
  }
  
  overall_dt <- if (length(overall_list) > 0) {
    rbindlist(overall_list, fill = TRUE)
  } else {
    data.table()
  }
  
  pairwise_dt <- if (length(pairwise_list) > 0) {
    rbindlist(pairwise_list, fill = TRUE)
  } else {
    data.table()
  }
  
  if (nrow(overall_dt) > 0) {
    overall_dt[, overall_FDR_across_pathways := p.adjust(overall_p, method = "BH")]
  }
  
  list(
    overall = overall_dt,
    pairwise = pairwise_dt,
    group_count = group_count
  )
}

# =========================================================
# 14 statistical framework C: independent pairwise Welch t-tests
# No three-group omnibus test is used in this framework.
# Figure stars/p values use the raw two-sided Welch p value.
# BH3 and across-pathway FDR are also saved in the CSV for reference.
# =========================================================
run_stats_Pairwise_Welch <- function(score_long, population_name, score_mode_name) {
  group_count <- unique(
    score_long[, .(n_FOV = uniqueN(FOV)), by = immunotype_BW_top_40]
  )
  group_count <- merge(
    data.table(immunotype_BW_top_40 = factor(immunotype_order, levels = immunotype_order)),
    group_count,
    by = "immunotype_BW_top_40",
    all.x = TRUE
  )
  group_count[is.na(n_FOV), n_FOV := 0L]
  
  pair_defs <- data.table(
    group1 = c("ID", "ID", "IE"),
    group2 = c("IE", "IN", "IN")
  )
  
  pathways <- unique(score_long$pathway)
  pairwise_list <- list()
  
  for (pw in pathways) {
    d <- copy(score_long[pathway == pw & is.finite(score)])
    if (nrow(d) == 0) next
    d[, immunotype_BW_top_40 := factor(immunotype_BW_top_40, levels = immunotype_order)]
    
    pr <- copy(pair_defs)
    pr[, n_group1 := vapply(group1, function(g) length(d[immunotype_BW_top_40 == g, score]), integer(1))]
    pr[, n_group2 := vapply(group2, function(g) length(d[immunotype_BW_top_40 == g, score]), integer(1))]
    pr[, mean_group1 := vapply(group1, function(g) mean(d[immunotype_BW_top_40 == g, score], na.rm = TRUE), numeric(1))]
    pr[, mean_group2 := vapply(group2, function(g) mean(d[immunotype_BW_top_40 == g, score], na.rm = TRUE), numeric(1))]
    pr[, median_group1 := vapply(group1, function(g) median(d[immunotype_BW_top_40 == g, score], na.rm = TRUE), numeric(1))]
    pr[, median_group2 := vapply(group2, function(g) median(d[immunotype_BW_top_40 == g, score], na.rm = TRUE), numeric(1))]
    pr[, mean_difference_group1_minus_group2 := mean_group1 - mean_group2]
    
    pr[, p_raw := vapply(seq_len(.N), function(j) {
      x <- d[immunotype_BW_top_40 == group1[j], score]
      y <- d[immunotype_BW_top_40 == group2[j], score]
      if (length(x) < min_FOV_per_group_for_stats || length(y) < min_FOV_per_group_for_stats) {
        return(NA_real_)
      }
      tryCatch(
        stats::t.test(x = x, y = y, alternative = "two.sided", var.equal = FALSE)$p.value,
        error = function(e) NA_real_
      )
    }, numeric(1))]
    
    pr[, p_BH3 := p.adjust(p_raw, method = "BH")]
    pr[, p_display := p_raw]
    pr[, contrast_std := paste0(group1, "_vs_", group2)]
    pr[, population := population_name]
    pr[, score_mode := score_mode_name]
    pr[, stat_framework := "Pairwise_Welch_rawP"]
    pr[, pathway := pw]
    pairwise_list[[length(pairwise_list) + 1]] <- pr
  }
  
  pairwise_dt <- if (length(pairwise_list) > 0) rbindlist(pairwise_list, fill = TRUE) else data.table()
  if (nrow(pairwise_dt) > 0) {
    pairwise_dt[, p_FDR_across_pathways_same_contrast := p.adjust(p_raw, method = "BH"), by = contrast_std]
  }
  
  list(
    overall = data.table(),
    pairwise = pairwise_dt,
    group_count = group_count
  )
}

# =========================================================
# 15 annotation helpers for all statistical frameworks
# =========================================================
get_pair_p <- function(pair_dt, pathway_name, a, b) {
  key1 <- paste0(a, "_vs_", b)
  key2 <- paste0(b, "_vs_", a)
  
  z <- pair_dt[pathway == pathway_name & contrast_std == key1, p_display]
  if (length(z) > 0) return(as.numeric(z[1]))
  
  z <- pair_dt[pathway == pathway_name & contrast_std == key2, p_display]
  if (length(z) > 0) return(as.numeric(z[1]))
  
  NA_real_
}

build_plot_annotations <- function(score_long, stats_obj) {
  if (nrow(stats_obj$pairwise) == 0) {
    return(list(overall = data.table(), pairwise = data.table()))
  }
  
  heights <- score_long[
    ,
    .(
      ymin = min(score, na.rm = TRUE),
      ymax = max(score, na.rm = TRUE)
    ),
    by = pathway
  ]
  
  heights[, yrange := pmax(ymax - ymin, 1e-6)]
  heights[, base_y := ymax + 0.10 * yrange]
  heights[, step_y := 0.14 * yrange]
  heights[, tip_y := 0.035 * yrange]
  
  overall_ann <- data.table()
  if (nrow(stats_obj$overall) > 0) {
    overall_ann <- merge(
      stats_obj$overall[, .(pathway, overall_p)],
      heights[, .(pathway, yrange, base_y, step_y)],
      by = "pathway",
      all.x = TRUE
    )
    overall_ann[, x := 2]
    overall_ann[, y := base_y + 3.15 * step_y]
    overall_ann[, label := paste0("overall p=", vapply(overall_p, format_p, character(1)))]
  }
  
  # Bracket order is chosen to minimize overlap:
  # ID-IE lower, IE-IN middle, ID-IN upper.
  pair_defs <- data.table(
    group1 = c("ID", "IE", "ID"),
    group2 = c("IE", "IN", "IN"),
    x = c(1, 2, 1),
    xend = c(2, 3, 3),
    idx = c(1, 2, 3)
  )
  
  pw_list <- lapply(seq_len(nrow(heights)), function(i) {
    pw <- heights$pathway[i]
    h <- heights[i]
    
    tmp <- copy(pair_defs)
    tmp[, pathway := pw]
    
    tmp[, p_value := vapply(seq_len(.N), function(j) {
      get_pair_p(stats_obj$pairwise, pw, group1[j], group2[j])
    }, numeric(1))]
    
    tmp[, y := h$base_y + (idx - 1) * h$step_y]
    tmp[, y_label := y + 0.018 * h$yrange]
    tmp[, xmid := (x + xend) / 2]
    tmp[, tip_y := h$tip_y]
    tmp[, label := paste0(
      "p=",
      vapply(p_value, format_p, character(1)),
      " ",
      vapply(p_value, p_to_star, character(1))
    )]
    
    tmp
  })
  
  pairwise_ann <- rbindlist(pw_list, fill = TRUE)
  
  list(
    overall = overall_ann,
    pairwise = pairwise_ann
  )
}

# =========================================================
# 15 one population page series inside an already-open PDF
# =========================================================
print_violin_pages <- function(
    score_long,
    stats_obj,
    population_name,
    score_mode_name,
    stat_framework_name
) {
  if (nrow(score_long) == 0) return(invisible(NULL))
  
  pathway_order <- if (nrow(stats_obj$overall) > 0) {
    stats_obj$overall[order(overall_p), pathway]
  } else if (nrow(stats_obj$pairwise) > 0) {
    ord_dt <- stats_obj$pairwise[
      ,
      .(min_pairwise_p = if (all(!is.finite(p_display))) Inf else min(p_display[is.finite(p_display)], na.rm = TRUE)),
      by = pathway
    ]
    ord_dt[order(min_pairwise_p, pathway), pathway]
  } else {
    sort(unique(score_long$pathway))
  }
  
  ann_obj <- build_plot_annotations(score_long, stats_obj)
  n_pages <- ceiling(length(pathway_order) / violin_per_page)
  
  for (page_i in seq_len(n_pages)) {
    idx1 <- (page_i - 1) * violin_per_page + 1
    idx2 <- min(page_i * violin_per_page, length(pathway_order))
    page_pw <- pathway_order[idx1:idx2]
    
    plot_dt <- copy(score_long[pathway %in% page_pw])
    plot_dt[, pathway := factor(pathway, levels = page_pw)]
    plot_dt[, x_num := match(as.character(immunotype_BW_top_40), immunotype_order)]
    
    y_label <- switch(
      score_mode_name,
      "normalized_ssGSEA" = "Normalized ssGSEA score (FOV-level)",
      "pathway_robust_z" = "Pathway-wise robust z-score",
      "FOV_centered" = "FOV-centered ssGSEA score",
      score_mode_name
    )
    
    stat_subtitle <- switch(
      stat_framework_name,
      "LM_Tukey" = "Overall partial F-test; pairwise Tukey-adjusted estimated marginal means",
      "KW_WilcoxonBH" = "Overall Kruskal-Wallis; pairwise Wilcoxon with BH correction across 3 comparisons",
      "Pairwise_Welch_rawP" = "Independent two-group Welch t-tests; raw two-sided pairwise p values; no omnibus three-group test",
      stat_framework_name
    )
    
    p <- ggplot(
      plot_dt,
      aes(
        x = x_num,
        y = score,
        fill = immunotype_BW_top_40,
        group = immunotype_BW_top_40
      )
    ) +
      geom_violin(
        scale = "width",
        adjust = 1,
        trim = FALSE,
        color = NA,
        alpha = 0.86
      ) +
      geom_boxplot(
        width = 0.18,
        outlier.shape = NA,
        alpha = 0.30,
        color = "black",
        linewidth = 0.35
      ) +
      stat_summary(
        fun = median,
        geom = "point",
        size = 1.35,
        shape = 21,
        fill = "white",
        color = "black",
        stroke = 0.35
      ) +
      geom_jitter(
        width = jitter_width,
        size = jitter_size,
        alpha = 0.56,
        color = "black"
      ) +
      scale_fill_manual(
        values = immunotype_cols,
        breaks = immunotype_order,
        drop = FALSE,
        guide = "none"
      ) +
      scale_x_continuous(
        breaks = 1:3,
        labels = immunotype_order,
        limits = c(0.5, 3.5),
        expand = expansion(mult = c(0.02, 0.02))
      ) +
      facet_wrap(
        ~ pathway,
        scales = "free_y",
        ncol = violin_ncol
      ) +
      labs(
        x = NULL,
        y = y_label,
        title = paste0(
          population_name,
          " T-cell signalling — ",
          score_mode_name,
          " — ",
          stat_framework_name,
          " — page ", page_i, "/", n_pages
        ),
        subtitle = stat_subtitle
      ) +
      theme_bw(base_size = 13) +
      theme(
        strip.background = element_rect(fill = "grey92", color = NA),
        strip.text = element_text(face = "bold", size = 9),
        axis.text.x = element_text(size = 10, face = "bold"),
        axis.text.y = element_text(size = 8),
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 9.5),
        plot.margin = margin(8, 22, 8, 6),
        panel.spacing.y = grid::unit(1.05, "lines")
      )
    
    if (nrow(ann_obj$pairwise) > 0) {
      seg_page <- copy(ann_obj$pairwise[pathway %in% page_pw])
      seg_page[, pathway := factor(pathway, levels = page_pw)]
      
      p <- p +
        geom_segment(
          data = seg_page,
          aes(x = x, xend = xend, y = y, yend = y),
          inherit.aes = FALSE,
          linewidth = 0.34
        ) +
        geom_segment(
          data = seg_page,
          aes(x = x, xend = x, y = y - tip_y, yend = y),
          inherit.aes = FALSE,
          linewidth = 0.34
        ) +
        geom_segment(
          data = seg_page,
          aes(x = xend, xend = xend, y = y - tip_y, yend = y),
          inherit.aes = FALSE,
          linewidth = 0.34
        ) +
        geom_text(
          data = seg_page,
          aes(x = xmid, y = y_label, label = label),
          inherit.aes = FALSE,
          size = 2.25,
          vjust = 0,
          fontface = "bold"
        )
    }
    
    if (nrow(ann_obj$overall) > 0) {
      ov_page <- copy(ann_obj$overall[pathway %in% page_pw])
      ov_page[, pathway := factor(pathway, levels = page_pw)]
      
      p <- p +
        geom_text(
          data = ov_page,
          aes(x = x, y = y, label = label),
          inherit.aes = FALSE,
          size = 2.35,
          vjust = 0,
          fontface = "bold"
        )
    }
    
    p <- p + coord_cartesian(clip = "off")
    
    print(p)
  }
  
  invisible(NULL)
}

# =========================================================
# 17 group-mean heatmap
# Columns = ID / IE / IN; rows = pathways; each tile starts from the group mean.
# Rows are z-scored across the three group means for visualization.
# =========================================================
make_group_mean_heatmap <- function(score_long, stats_obj_for_order, population_name, score_mode_name) {
  mean_dt <- score_long[
    ,
    .(
      mean_score = mean(score, na.rm = TRUE),
      median_score = median(score, na.rm = TRUE),
      n_FOV = uniqueN(FOV)
    ),
    by = .(pathway, immunotype_BW_top_40)
  ]
  
  wide <- dcast(mean_dt, pathway ~ immunotype_BW_top_40, value.var = "mean_score")
  for (grp in immunotype_order) {
    if (!grp %in% names(wide)) wide[, (grp) := NA_real_]
  }
  
  mat <- as.matrix(wide[, ..immunotype_order])
  rownames(mat) <- wide$pathway
  z <- t(scale(t(mat)))
  z[!is.finite(z)] <- 0
  
  pathway_order <- if (!is.null(stats_obj_for_order) && nrow(stats_obj_for_order$overall) > 0) {
    stats_obj_for_order$overall[order(overall_p), pathway]
  } else {
    rownames(z)
  }
  pathway_order <- pathway_order[pathway_order %in% rownames(z)]
  
  heat_dt <- as.data.table(as.table(z))
  setnames(heat_dt, c("pathway", "immunotype", "z"))
  heat_dt[, immunotype := factor(as.character(immunotype), levels = immunotype_order)]
  heat_dt[, pathway := factor(as.character(pathway), levels = rev(pathway_order))]
  
  p <- ggplot(heat_dt, aes(x = immunotype, y = pathway, fill = z)) +
    geom_tile(color = "white", linewidth = 0.45) +
    scale_fill_gradient2(
      low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
      midpoint = 0, name = "row z-score"
    ) +
    scale_x_discrete(limits = immunotype_order, drop = FALSE) +
    labs(
      x = NULL, y = NULL,
      title = paste0(population_name, " T-cell signalling group-mean heatmap"),
      subtitle = paste0(score_mode_name, "; mean score by immunotype; pathway rows z-scaled for display")
    ) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(size = 10, face = "bold"),
      axis.text.y = element_text(size = 8),
      plot.title = element_text(face = "bold", size = 11),
      plot.subtitle = element_text(size = 9)
    )
  
  save_pdf(
    p,
    paste0("03_4_3_", safe_tag(population_name), "_", safe_tag(score_mode_name), "_mean_heatmap_ID_IE_IN"),
    w = mean_heatmap_width,
    h = max(mean_heatmap_min_height, mean_heatmap_height_per_pathway * uniqueN(mean_dt$pathway))
  )
  
  save_dat(
    mean_dt,
    paste0("03_4_3_", safe_tag(population_name), "_", safe_tag(score_mode_name), "_group_mean_scores"),
    "csv"
  )
  
  invisible(list(mean_table = mean_dt, heatmap_table = heat_dt))
}

# =========================================================
# 18 run one population: pseudobulk + 3 scores + 3 statistical frameworks
# =========================================================
run_population <- function(population_name, selector_fun, min_cells_per_FOV) {
  message("")
  message("=========================================================")
  message("03_4_3 population: ", population_name)
  message("=========================================================")
  
  pop_meta <- copy(meta[selector_fun(BC6k_refined)])
  
  if (nrow(pop_meta) == 0) {
    stop("No cells found for population: ", population_name)
  }
  
  pb_obj <- build_FOV_pseudobulk(
    meta_subset = pop_meta,
    min_cells_per_FOV = min_cells_per_FOV,
    population_name = population_name
  )
  
  save_dat(
    pb_obj$fov_anno,
    paste0("03_4_3_", safe_tag(population_name), "_FOV_annotation"),
    "csv"
  )
  
  score_obj <- score_FOV_ssGSEA(pb_obj$pb_mat)
  
  pathway_used_dt <- rbindlist(
    lapply(names(score_obj$pathways_used), function(pw) {
      data.table(
        population = population_name,
        pathway = pw,
        n_genes_used = length(score_obj$pathways_used[[pw]]),
        genes_used = paste(score_obj$pathways_used[[pw]], collapse = ";")
      )
    })
  )
  
  save_dat(
    pathway_used_dt,
    paste0("03_4_3_", safe_tag(population_name), "_pathway_genes_used"),
    "csv"
  )
  
  score_modes <- c("normalized_ssGSEA", "pathway_robust_z", "FOV_centered")
  
  score_long_list <- list()
  LM_stats_list <- list()
  KW_stats_list <- list()
  Welch_stats_list <- list()
  
  for (score_mode_i in score_modes) {
    score_mat_i <- score_obj[[score_mode_i]]
    
    save_dat(
      score_mat_i,
      paste0("03_4_3_", safe_tag(population_name), "_matrix_", score_mode_i),
      "rds"
    )
    
    long_i <- score_matrix_to_long(
      score_mat = score_mat_i,
      fov_anno = pb_obj$fov_anno,
      population_name = population_name,
      score_mode = score_mode_i
    )
    
    save_dat(
      long_i,
      paste0("03_4_3_", safe_tag(population_name), "_long_", score_mode_i),
      "csv"
    )
    
    save_dat(
      long_i,
      paste0("03_4_3_", safe_tag(population_name), "_long_", score_mode_i),
      "rds"
    )
    
    LM_i <- run_stats_LM_Tukey(
      score_long = long_i,
      population_name = population_name,
      score_mode_name = score_mode_i
    )
    
    KW_i <- run_stats_KW_WilcoxonBH(
      score_long = long_i,
      population_name = population_name,
      score_mode_name = score_mode_i
    )
    
    Welch_i <- run_stats_Pairwise_Welch(
      score_long = long_i,
      population_name = population_name,
      score_mode_name = score_mode_i
    )
    
    if (nrow(LM_i$overall) > 0) {
      save_dat(
        LM_i$overall,
        paste0("03_4_3_", safe_tag(population_name), "_", score_mode_i, "_LM_Tukey_overall"),
        "csv"
      )
    }
    
    if (nrow(LM_i$pairwise) > 0) {
      save_dat(
        LM_i$pairwise,
        paste0("03_4_3_", safe_tag(population_name), "_", score_mode_i, "_LM_Tukey_pairwise"),
        "csv"
      )
    }
    
    if (nrow(KW_i$overall) > 0) {
      save_dat(
        KW_i$overall,
        paste0("03_4_3_", safe_tag(population_name), "_", score_mode_i, "_KW_WilcoxonBH_overall"),
        "csv"
      )
    }
    
    if (nrow(KW_i$pairwise) > 0) {
      save_dat(
        KW_i$pairwise,
        paste0("03_4_3_", safe_tag(population_name), "_", score_mode_i, "_KW_WilcoxonBH_pairwise"),
        "csv"
      )
    }
    
    if (nrow(Welch_i$pairwise) > 0) {
      save_dat(
        Welch_i$pairwise,
        paste0("03_4_3_", safe_tag(population_name), "_", score_mode_i, "_Pairwise_Welch_rawP_pairwise"),
        "csv"
      )
    }
    
    score_long_list[[score_mode_i]] <- long_i
    LM_stats_list[[score_mode_i]] <- LM_i
    KW_stats_list[[score_mode_i]] <- KW_i
    Welch_stats_list[[score_mode_i]] <- Welch_i
  }
  
  list(
    score_long = score_long_list,
    LM = LM_stats_list,
    KW = KW_stats_list,
    Welch = Welch_stats_list,
    fov_anno = pb_obj$fov_anno,
    pathway_used = pathway_used_dt
  )
}

# =========================================================
# 19 run Total T and CD8 T
# =========================================================
res_Total_T <- run_population(
  population_name = "Total_T",
  selector_fun = function(x) x %in% total_T_labels,
  min_cells_per_FOV = min_cells_per_FOV_total_T
)

res_CD8_T <- run_population(
  population_name = "CD8_T",
  selector_fun = function(x) x == CD8_label,
  min_cells_per_FOV = min_cells_per_FOV_CD8
)

# =========================================================
# 20 create nine violin figure sets
# 3 score systems x 3 statistical frameworks.
# Each PDF contains Total_T page(s), then CD8_T page(s).
# =========================================================
score_modes <- c("normalized_ssGSEA", "pathway_robust_z", "FOV_centered")
stat_modes <- c("LM_Tukey", "KW_WilcoxonBH", "Pairwise_Welch_rawP")

for (score_mode_i in score_modes) {
  for (stat_mode_i in stat_modes) {
    pdf_file <- file.path(
      fig_dir,
      paste0(
        "03_4_3_SCORE_", score_mode_i, "__STAT_", stat_mode_i,
        "__TotalT_and_CD8_violin.pdf"
      )
    )
    if (file.exists(pdf_file)) file.remove(pdf_file)
    grDevices::cairo_pdf(filename = pdf_file, width = violin_pdf_width, height = violin_pdf_height)
    
    if (stat_mode_i == "LM_Tukey") {
      total_stats <- res_Total_T$LM[[score_mode_i]]
      cd8_stats <- res_CD8_T$LM[[score_mode_i]]
    } else if (stat_mode_i == "KW_WilcoxonBH") {
      total_stats <- res_Total_T$KW[[score_mode_i]]
      cd8_stats <- res_CD8_T$KW[[score_mode_i]]
    } else {
      total_stats <- res_Total_T$Welch[[score_mode_i]]
      cd8_stats <- res_CD8_T$Welch[[score_mode_i]]
    }
    
    print_violin_pages(
      score_long = res_Total_T$score_long[[score_mode_i]],
      stats_obj = total_stats,
      population_name = "Total_T",
      score_mode_name = score_mode_i,
      stat_framework_name = stat_mode_i
    )
    print_violin_pages(
      score_long = res_CD8_T$score_long[[score_mode_i]],
      stats_obj = cd8_stats,
      population_name = "CD8_T",
      score_mode_name = score_mode_i,
      stat_framework_name = stat_mode_i
    )
    
    grDevices::dev.off()
    message("Saved figure set: ", pdf_file)
  }
}

# =========================================================
# 21 restored mean heatmaps: 3 scores x 2 populations = 6 PDFs
# =========================================================
for (score_mode_i in score_modes) {
  make_group_mean_heatmap(
    res_Total_T$score_long[[score_mode_i]],
    res_Total_T$LM[[score_mode_i]],
    "Total_T",
    score_mode_i
  )
  make_group_mean_heatmap(
    res_CD8_T$score_long[[score_mode_i]],
    res_CD8_T$LM[[score_mode_i]],
    "CD8_T",
    score_mode_i
  )
}

# =========================================================
# 22 combined master statistics
# =========================================================
collect_stat_tables <- function(which_table = c("overall", "pairwise")) {
  which_table <- match.arg(which_table)
  
  out <- list()
  
  for (pop_name in c("Total_T", "CD8_T")) {
    res_pop <- if (pop_name == "Total_T") res_Total_T else res_CD8_T
    
    for (score_mode_i in score_modes) {
      for (stat_mode_i in stat_modes) {
        obj <- if (stat_mode_i == "LM_Tukey") {
          res_pop$LM[[score_mode_i]]
        } else if (stat_mode_i == "KW_WilcoxonBH") {
          res_pop$KW[[score_mode_i]]
        } else {
          res_pop$Welch[[score_mode_i]]
        }
        
        tab <- obj[[which_table]]
        if (nrow(tab) > 0) {
          out[[length(out) + 1]] <- tab
        }
      }
    }
  }
  
  if (length(out) == 0) return(data.table())
  rbindlist(out, fill = TRUE)
}

master_overall <- collect_stat_tables("overall")
master_pairwise <- collect_stat_tables("pairwise")

if (nrow(master_overall) > 0) {
  save_dat(
    master_overall,
    "03_4_3_ALL_scores_ALL_stats_overall_MASTER",
    "csv"
  )
  save_dat(
    master_overall,
    "03_4_3_ALL_scores_ALL_stats_overall_MASTER",
    "rds"
  )
}

if (nrow(master_pairwise) > 0) {
  save_dat(
    master_pairwise,
    "03_4_3_ALL_scores_ALL_stats_pairwise_MASTER",
    "csv"
  )
  save_dat(
    master_pairwise,
    "03_4_3_ALL_scores_ALL_stats_pairwise_MASTER",
    "rds"
  )
}

# =========================================================
# 23 compact comparison table
# Lets the user compare whether pathway conclusions are stable across
# all 3 score systems and both statistical frameworks.
# =========================================================
compact_overall <- copy(master_overall)

if (nrow(compact_overall) > 0) {
  compact_overall[, rank_high_to_low := apply(
    cbind(mean_ID, mean_IE, mean_IN),
    1,
    function(x) {
      names(x) <- c("ID", "IE", "IN")
      x <- x[is.finite(x)]
      if (length(x) == 0) return(NA_character_)
      paste(names(sort(x, decreasing = TRUE)), collapse = " > ")
    }
  )]
  
  setorder(
    compact_overall,
    population,
    pathway,
    score_mode,
    stat_framework
  )
  
  save_dat(
    compact_overall,
    "03_4_3_score_stat_framework_comparison_summary",
    "csv"
  )
}

# =========================================================
# 24 parameter record
# =========================================================
parameter_dt <- data.table(
  parameter = c(
    "biological_replicate", "metadata_file", "expression_file", "pathway_file",
    "immunotype_order", "immunotype_colors", "total_T_definition", "CD8_definition",
    "min_cells_per_FOV_total_T", "min_cells_per_FOV_CD8",
    "score_1", "score_2", "score_3",
    "stat_framework_A", "stat_framework_A_pairwise",
    "stat_framework_B", "stat_framework_B_pairwise",
    "stat_framework_C", "stat_framework_C_pairwise", "stat_framework_C_figure_p",
    "violin_figure_sets", "mean_heatmaps", "figure_height_inches"
  ),
  value = c(
    "one patient = one FOV", meta_file, norm_file, xlsx_path,
    paste(immunotype_order, collapse = " -> "),
    "ID=#000000; IE=#0000FF; IN=#FF0000",
    paste(total_T_labels, collapse = ";"), CD8_label,
    as.character(min_cells_per_FOV_total_T), as.character(min_cells_per_FOV_CD8),
    "normalized ssGSEA returned by GSVA ssGSEA with normalization enabled",
    "within each pathway across FOVs: (score - median) / MAD; SD fallback if MAD=0",
    "within each FOV across pathways: (score - FOV median) / FOV MAD",
    "linear model; slide as fixed covariate when estimable; overall partial F-test",
    "estimated marginal means; Tukey-adjusted pairwise p values",
    "Kruskal-Wallis overall test",
    "pairwise Wilcoxon rank-sum; BH across 3 pre-defined comparisons within each pathway",
    "independent pairwise-only analysis; no omnibus three-group test",
    "two-sided Welch t-test for ID-vs-IE, ID-vs-IN and IE-vs-IN separately",
    "raw Welch p displayed; BH3 and across-pathway FDR also saved in CSV",
    "9 PDFs = 3 scores x 3 statistical frameworks; each PDF contains Total_T and CD8_T",
    "6 PDFs = 3 scores x 2 populations; ID/IE/IN group means with row z-score display",
    as.character(violin_pdf_height)
  )
)

save_dat(parameter_dt, "03_4_3_parameters_v4", "csv")

message("")
message("03_4_3 v4 finished successfully.")
message("Three score systems: normalized_ssGSEA; pathway_robust_z; FOV_centered")
message("Three statistical frameworks: LM_Tukey; KW_WilcoxonBH; Pairwise_Welch_rawP")
message("Pairwise_Welch_rawP uses separate two-sided Welch t-tests and no omnibus three-group test.")
message("Nine violin figure-set PDFs are written to: ", fig_dir)
message("Six mean-heatmap PDFs are also written to: ", fig_dir)
message("Violin PDF size: ", violin_pdf_width, " x ", violin_pdf_height, " inches")
message("Colors are fixed: ID pure black, IE pure blue, IN pure red.")

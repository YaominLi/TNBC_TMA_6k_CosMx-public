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
raw_dir <- file.path(proj, "raw")
out_dir <- file.path(proj, "results", "03_3", "data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

save_dat_overwrite <- function(x, filename, ext = c("csv", "rds")) {
  ext <- match.arg(ext)
  out <- file.path(out_dir, paste0(filename, ".", ext))
  if (file.exists(out)) file.remove(out)
  if (ext == "csv") write.csv(x, out, row.names = FALSE)
  if (ext == "rds") saveRDS(x, out)
  message("Saved: ", out)
  invisible(out)
}

# =========================================================
# 02 packages
# =========================================================
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
})

# =========================================================
# 03 user input: optional custom gene-set score
# if score_genes is empty, no extra score column will be added
# =========================================================
score_name <- "Neuroactive ligand-receptor interaction"   
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

# scoring cells:
# "tumor_only" = only tumor cells within each FOV
# "all_cells"  = all cells within each FOV
score_basis <- "tumor_only"

# minimum number of cells required in a FOV for scoring
min_cells_for_scoring <- 10

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
  "BC6k_refined", "immunotype_BW_top_40"
)
miss <- setdiff(need_cols, names(dt))
if (length(miss) > 0) {
  stop("Missing required columns: ", paste(miss, collapse = ", "))
}

dt <- dt[
  !is.na(FOV) &
    !is.na(Run_Tissue_name) &
    !is.na(BC6k_refined) &
    !is.na(immunotype_BW_top_40)
]

dt[, FOV := as.character(FOV)]
dt[, Run_Tissue_name := as.character(Run_Tissue_name)]
dt[, BC6k_refined := as.character(BC6k_refined)]
dt[, immunotype_BW_top_40 := as.character(immunotype_BW_top_40)]

# keep only the current three-class FOV labels
dt <- dt[immunotype_BW_top_40 %in% c("ID", "IE", "IN")]

# =========================================================
# 06 one row per FOV basic annotation
# =========================================================
fov_anno <- unique(dt[, .(
  FOV,
  Run_Tissue_name,
  immunotype_BW_top_40
)])

chk <- fov_anno[, .N, by = FOV][N > 1]
if (nrow(chk) > 0) {
  stop("Some FOVs have conflicting Run_Tissue_name or immunotype_BW_top_40 assignments.")
}

# =========================================================
# 07 basic counts per FOV
# =========================================================
fov_counts_all <- dt[, .(n_cells_all = .N), by = FOV]
fov_counts_tumor <- dt[BC6k_refined == "tumor", .(n_cells_tumor = .N), by = FOV]
fov_counts_nontumor <- dt[BC6k_refined != "tumor", .(n_cells_nontumor = .N), by = FOV]

fov_summary <- Reduce(
  function(x, y) merge(x, y, by = "FOV", all = TRUE),
  list(fov_anno, fov_counts_all, fov_counts_tumor, fov_counts_nontumor)
)

for (cc in c("n_cells_all", "n_cells_tumor", "n_cells_nontumor")) {
  fov_summary[is.na(get(cc)), (cc) := 0L]
}

# =========================================================
# 08 all-cell composition by BC6k_refined
# =========================================================
all_long <- dt[, .N, by = .(FOV, BC6k_refined)]
all_long <- merge(
  all_long,
  fov_summary[, .(FOV, n_cells_all)],
  by = "FOV",
  all.x = TRUE
)
all_long[, pct_all := fifelse(n_cells_all > 0, N / n_cells_all, NA_real_)]

all_wide <- dcast(
  all_long,
  FOV ~ BC6k_refined,
  value.var = "pct_all",
  fill = 0
)
if (ncol(all_wide) > 1) {
  old_cols <- setdiff(names(all_wide), "FOV")
  setnames(all_wide, old_cols, paste0("pct_all__", old_cols))
}

# =========================================================
# 09 non-tumor composition by BC6k_refined
# =========================================================
nontumor_dt <- dt[BC6k_refined != "tumor"]

nontumor_long <- nontumor_dt[, .N, by = .(FOV, BC6k_refined)]
nontumor_long <- merge(
  nontumor_long,
  fov_summary[, .(FOV, n_cells_nontumor)],
  by = "FOV",
  all.x = TRUE
)
nontumor_long[, pct_nontumor := fifelse(n_cells_nontumor > 0, N / n_cells_nontumor, NA_real_)]

nontumor_wide <- dcast(
  nontumor_long,
  FOV ~ BC6k_refined,
  value.var = "pct_nontumor",
  fill = 0
)
if (ncol(nontumor_wide) > 1) {
  old_cols <- setdiff(names(nontumor_wide), "FOV")
  setnames(nontumor_wide, old_cols, paste0("pct_nontumor__", old_cols))
}

# =========================================================
# 10 merge one-row-per-FOV summary
# =========================================================
fov_summary <- merge(fov_summary, all_wide, by = "FOV", all.x = TRUE)
fov_summary <- merge(fov_summary, nontumor_wide, by = "FOV", all.x = TRUE)

# =========================================================
# 11 optional custom gene-set score
# scoring logic:
#   1) tumor-only or all-cells within each FOV
#   2) FOV-level pseudobulk = rowMeans across selected cells
#   3) for each gene, z-score within Run_Tissue_name
#   4) final FOV score = mean z-score across genes in the set
# =========================================================
if (length(score_genes) > 0) {
  
  norm_file <- file.path(raw_dir, "TNBC_6k_norm_Bogang.RDS")
  if (!file.exists(norm_file)) {
    stop("Missing file: ", norm_file)
  }
  norm <- readRDS(norm_file)
  
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
  
  # ensure genes x cells
  if (!any(dt$cell_id %in% colnames(expr)) && any(dt$cell_id %in% rownames(expr))) {
    expr <- t(expr)
  }
  
  common_cells <- intersect(dt$cell_id, colnames(expr))
  common_cells <- common_cells[match(common_cells, dt$cell_id)]
  common_cells <- common_cells[!is.na(common_cells)]
  
  expr <- expr[, common_cells, drop = FALSE]
  dt_expr <- dt[match(common_cells, cell_id)]
  
  if (!identical(colnames(expr), dt_expr$cell_id)) {
    stop("Expression matrix cell order does not match metadata.")
  }
  
  # choose scoring cells
  if (score_basis == "tumor_only") {
    score_dt <- dt_expr[BC6k_refined == "tumor"]
  } else if (score_basis == "all_cells") {
    score_dt <- copy(dt_expr)
  } else {
    stop("score_basis must be either 'tumor_only' or 'all_cells'.")
  }
  
  # minimum cells per FOV for scoring
  fov_n_score_cells <- score_dt[, .N, by = .(FOV, Run_Tissue_name)]
  valid_fov_score <- fov_n_score_cells[N >= min_cells_for_scoring, FOV]
  score_dt <- score_dt[FOV %in% valid_fov_score]
  
  if (nrow(score_dt) > 0) {
    
    expr_score <- expr[, score_dt$cell_id, drop = FALSE]
    if (!identical(colnames(expr_score), score_dt$cell_id)) {
      stop("Scoring expression matrix cell order mismatch.")
    }
    
    score_genes <- unique(as.character(score_genes))
    genes_present <- intersect(score_genes, rownames(expr_score))
    genes_missing <- setdiff(score_genes, genes_present)
    
    gene_use_dt <- data.table(
      score_name = score_name,
      input_gene = score_genes,
      in_expr = score_genes %in% rownames(expr_score)
    )
    save_dat_overwrite(as.data.frame(gene_use_dt), paste0("03_3_", score_name, "_gene_usage"), ext = "csv")
    
    if (length(genes_present) == 0) {
      warning("No input genes present in expression matrix. Score column will be NA.")
      fov_summary[, (score_name) := NA_real_]
    } else {
      
      expr_use <- expr_score[genes_present, , drop = FALSE]
      
      # FOV-level pseudobulk
      fov_list <- split(seq_len(nrow(score_dt)), score_dt$FOV)
      
      pb_list <- lapply(names(fov_list), function(fv) {
        idx <- fov_list[[fv]]
        cells_use <- score_dt$cell_id[idx]
        mat <- expr_use[, cells_use, drop = FALSE]
        
        if (ncol(mat) == 1) {
          val <- as.numeric(mat[, 1])
        } else {
          val <- Matrix::rowMeans(mat)
        }
        
        data.table(
          FOV = fv,
          gene = genes_present,
          expr_mean = val
        )
      })
      
      pb_long <- rbindlist(pb_list)
      pb_long <- merge(
        pb_long,
        unique(score_dt[, .(FOV, Run_Tissue_name)]),
        by = "FOV",
        all.x = TRUE
      )
      
      # z-score within Run_Tissue_name for each gene
      zscore_safe <- function(x) {
        s <- sd(x, na.rm = TRUE)
        m <- mean(x, na.rm = TRUE)
        if (is.na(s) || s == 0) return(rep(0, length(x)))
        (x - m) / s
      }
      
      pb_long[, gene_z_within_slide := zscore_safe(expr_mean), by = .(Run_Tissue_name, gene)]
      
      fov_score <- pb_long[, .(
        score = mean(gene_z_within_slide, na.rm = TRUE),
        n_genes_used = uniqueN(gene)
      ), by = .(FOV, Run_Tissue_name)]
      
      setnames(fov_score, "score", score_name)
      
      save_dat_overwrite(as.data.frame(pb_long), paste0("03_3_", score_name, "_pseudobulk_long"), ext = "csv")
      save_dat_overwrite(as.data.frame(fov_score), paste0("03_3_", score_name, "_FOV_score"), ext = "csv")
      
      fov_summary <- merge(
        fov_summary,
        fov_score[, .(FOV, Run_Tissue_name, n_genes_used, get(score_name))],
        by = c("FOV", "Run_Tissue_name"),
        all.x = TRUE
      )
      setnames(fov_summary, "V4", score_name)
      setnames(fov_summary, "n_genes_used", paste0(score_name, "__n_genes_used"))
    }
    
  } else {
    warning("No FOV passed the minimum cell threshold for scoring. Score column will be NA.")
    fov_summary[, (score_name) := NA_real_]
    fov_summary[, (paste0(score_name, "__n_genes_used")) := NA_integer_]
  }
}

# =========================================================
# 12 column order
# =========================================================
fixed_cols <- c(
  "FOV", "Run_Tissue_name", "immunotype_BW_top_40",
  "n_cells_all", "n_cells_tumor", "n_cells_nontumor"
)

other_cols <- setdiff(names(fov_summary), fixed_cols)
setcolorder(fov_summary, c(fixed_cols, other_cols))

setorder(fov_summary, immunotype_BW_top_40, Run_Tissue_name, FOV)

# =========================================================
# 13 save
# =========================================================
save_dat_overwrite(as.data.frame(fov_summary), "03_3_FOV_summary_table", ext = "csv")
save_dat_overwrite(as.data.frame(fov_summary), "03_3_FOV_summary_table", ext = "rds")

message("03_3 finished successfully.")
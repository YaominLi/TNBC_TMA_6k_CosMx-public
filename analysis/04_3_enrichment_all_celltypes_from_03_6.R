rm(list = ls())
gc()
graphics.off()

options(stringsAsFactors = FALSE)
set.seed(1)

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

raw <- file.path(proj, "raw")
gene_set_dir <- file.path(proj, "gene_sets")

out_dir <- file.path(proj, "results", "04_3")
dat_dir <- file.path(out_dir, "data")
checkpoint_dir <- file.path(dat_dir, "checkpoints")

dir.create(dat_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

save_dat_overwrite <- function(x, filename, ext = c("csv", "rds")) {
  ext <- match.arg(ext)
  out <- file.path(dat_dir, paste0(filename, ".", ext))
  
  if (file.exists(out)) {
    ok_rm <- tryCatch(
      file.remove(out),
      warning = function(e) FALSE,
      error = function(e) FALSE
    )
    if (!isTRUE(ok_rm)) stop("Cannot overwrite existing file: ", out)
  }
  
  if (ext == "csv") {
    write.csv(x, out, row.names = FALSE)
  } else {
    saveRDS(x, out)
  }
  
  message("Saved data: ", out)
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
  library(limma)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(fgsea)
})

# =========================================================
# 03 analysis settings
# =========================================================

# ---------------------------------------------------------
# Cell types
# ---------------------------------------------------------
# 04_1 / 04_2 already provide the formal tumor-only analysis.
# Therefore 04_3 excludes tumor by default and analyzes other cell types.
include_tumor_in_04_3 <- FALSE

broad_celltype_col <- "BC6k_refined"

# NULL = use every broad cell type present in the current metadata.
broad_celltypes_to_run <- NULL

# Yi updated fine CD8 annotation.
include_cd8_subtypes <- TRUE
fine_celltype_col <- "celltype_new"

cd8_subtypes_to_run <- c(
  "cd8_cytotoxic",
  "cd8_exhausted",
  "cd8_naive",
  "cd8_tcm",
  "cd8_tem"
)

# ---------------------------------------------------------
# FOV inclusion thresholds
# ---------------------------------------------------------
# Tumor used 30 cells/FOV in 04_1 and 04_2.
# For other broad cell types, 10 cells/FOV is used.
# For fine CD8 subtypes, 5 cells/FOV is used because these are less abundant.
min_cells_per_fov_broad <- 10L
min_cells_per_fov_fine_cd8 <- 5L
min_cells_per_fov_tumor <- 30L

# Require all three immunotypes to have enough FOVs.
# This is important because the same fitted model is used for both
# pairwise and true one-vs-rest contrasts.
min_fov_per_immunotype <- 3L

immunotype_levels <- c("ID", "IE", "IN")

# ---------------------------------------------------------
# DEG thresholds
# ---------------------------------------------------------
padj_cutoff <- 0.05
logfc_cutoff <- 0.25

# ---------------------------------------------------------
# Enrichment
# ---------------------------------------------------------
gsea_padj_query_cutoff <- 0.20

gmt_files <- list(
  H  = file.path(gene_set_dir, "h.all.v2025.1.Hs.symbols.gmt"),
  C2 = file.path(gene_set_dir, "c2.all.v2025.1.Hs.symbols.gmt"),
  C5 = file.path(gene_set_dir, "c5.all.v2025.1.Hs.symbols.gmt"),
  C7 = file.path(gene_set_dir, "c7.all.v2025.1.Hs.symbols.gmt")
)

# Original 04_1 / 04_2 also ran clusterProfiler KEGG online.
# Keep it optional because KEGG annotation can require internet access.
# FALSE still retains KEGG-related curated gene sets within local C2 GMT.
run_online_kegg <- FALSE

# ---------------------------------------------------------
# Computation
# ---------------------------------------------------------
# Keep the enrichment code itself the same style as the already-working
# 04_1 / 04_2. Do not introduce a new nested parallel backend here.
# data.table threading is safe and does not change the statistical analysis.
n_dt_threads <- parallel::detectCores(logical = TRUE)
if (is.na(n_dt_threads) || n_dt_threads < 1L) n_dt_threads <- 1L
data.table::setDTthreads(min(8L, n_dt_threads))

# Reuse completed cell-type checkpoints only when intentionally resuming the same code/parameters.
# FALSE is the reproducibility-safe default because checkpoint filenames do not encode all inputs/settings.
resume_existing <- FALSE

# =========================================================
# 04 load data
# =========================================================
norm_file <- file.path(raw, "TNBC_6k_norm_Bogang.RDS")
meta_new_file <- file.path(
  proj,
  "results",
  "03_2",
  "data",
  "03_2_meta_with_immunotype_BW_top_40.rds"
)

if (!file.exists(norm_file)) stop("Missing expression file: ", norm_file)
if (!file.exists(meta_new_file)) stop("Missing 03_2 metadata file: ", meta_new_file)

norm <- readRDS(norm_file)
meta_new <- as.data.table(readRDS(meta_new_file))

need_cols <- c(
  "cell_id",
  "FOV",
  "Run_Tissue_name",
  "BC6k_refined",
  "immunotype_BW_top_40"
)

miss_cols <- setdiff(need_cols, names(meta_new))
if (length(miss_cols) > 0) {
  stop("Missing required metadata columns: ", paste(miss_cols, collapse = ", "))
}

if (include_cd8_subtypes && !fine_celltype_col %in% names(meta_new)) {
  warning(
    "Fine CD8 column '", fine_celltype_col,
    "' not found. Fine CD8 analyses will be skipped."
  )
  include_cd8_subtypes <- FALSE
}

for (nm in names(gmt_files)) {
  if (!file.exists(gmt_files[[nm]])) {
    stop("Missing GMT file: ", gmt_files[[nm]])
  }
}

# =========================================================
# 05 extract expression matrix
# same logic as 04_1 / 04_2
# =========================================================
get_expr <- function(norm) {
  if (is.matrix(norm) || inherits(norm, "Matrix")) return(norm)
  
  if (inherits(norm, "Seurat")) {
    assay <- if ("RNA" %in% names(norm@assays)) {
      "RNA"
    } else {
      Seurat::DefaultAssay(norm)
    }
    
    m <- tryCatch(
      Seurat::GetAssayData(norm, assay = assay, slot = "data"),
      error = function(e) {
        Seurat::GetAssayData(norm, assay = assay, slot = "counts")
      }
    )
    return(m)
  }
  
  if (is.list(norm)) {
    for (nm in c(
      "data", "expr", "expression", "mat",
      "matrix", "norm", "logcounts", "counts"
    )) {
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

if (
  !any(meta_new$cell_id %in% colnames(expr)) &&
  any(meta_new$cell_id %in% rownames(expr))
) {
  expr <- t(expr)
}

common_cells <- intersect(meta_new$cell_id, colnames(expr))
common_cells <- common_cells[match(common_cells, meta_new$cell_id)]
common_cells <- common_cells[!is.na(common_cells)]

expr <- expr[, common_cells, drop = FALSE]
meta_new <- meta_new[match(common_cells, meta_new$cell_id)]

stopifnot(identical(colnames(expr), meta_new$cell_id))

if (anyDuplicated(rownames(expr)) > 0) {
  expr_dt <- data.table(
    gene = rownames(expr),
    as.matrix(expr)
  )
  expr_dt <- expr_dt[, lapply(.SD, mean), by = gene]
  expr <- as.matrix(expr_dt[, -1])
  rownames(expr) <- expr_dt$gene
}

meta_new[, cell_id := as.character(cell_id)]
meta_new[, FOV := as.character(FOV)]
meta_new[, Run_Tissue_name := as.character(Run_Tissue_name)]
meta_new[, BC6k_refined := as.character(BC6k_refined)]
meta_new[, immunotype_BW_top_40 := as.character(immunotype_BW_top_40)]

if (fine_celltype_col %in% names(meta_new)) {
  meta_new[, (fine_celltype_col) := as.character(get(fine_celltype_col))]
}

meta_new <- meta_new[
  !is.na(FOV) &
    !is.na(Run_Tissue_name) &
    !is.na(immunotype_BW_top_40) &
    immunotype_BW_top_40 %in% immunotype_levels
]

# =========================================================
# 06 define cell-type analysis manifest
# =========================================================
if (is.null(broad_celltypes_to_run)) {
  broad_celltypes_to_run <- sort(
    unique(
      meta_new[
        !is.na(get(broad_celltype_col)),
        get(broad_celltype_col)
      ]
    )
  )
}

if (!include_tumor_in_04_3) {
  broad_celltypes_to_run <- setdiff(
    broad_celltypes_to_run,
    "tumor"
  )
}

analysis_manifest_list <- list()

for (ct in broad_celltypes_to_run) {
  min_cells_this <- if (ct == "tumor") {
    min_cells_per_fov_tumor
  } else {
    min_cells_per_fov_broad
  }
  
  analysis_manifest_list[[length(analysis_manifest_list) + 1]] <- data.table(
    annotation_level = "broad",
    celltype_col = broad_celltype_col,
    celltype = ct,
    min_cells_per_fov = min_cells_this
  )
}

if (include_cd8_subtypes) {
  available_fine <- sort(
    unique(
      meta_new[
        !is.na(get(fine_celltype_col)),
        get(fine_celltype_col)
      ]
    )
  )
  
  valid_cd8 <- intersect(
    cd8_subtypes_to_run,
    available_fine
  )
  
  missing_cd8 <- setdiff(
    cd8_subtypes_to_run,
    valid_cd8
  )
  
  if (length(missing_cd8) > 0) {
    warning(
      "Fine CD8 subtypes not found and skipped: ",
      paste(missing_cd8, collapse = ", ")
    )
  }
  
  for (ct in valid_cd8) {
    analysis_manifest_list[[length(analysis_manifest_list) + 1]] <- data.table(
      annotation_level = "fine_CD8",
      celltype_col = fine_celltype_col,
      celltype = ct,
      min_cells_per_fov = min_cells_per_fov_fine_cd8
    )
  }
}

analysis_manifest <- rbindlist(
  analysis_manifest_list,
  fill = TRUE
)

setorder(
  analysis_manifest,
  annotation_level,
  celltype
)

save_dat_overwrite(
  as.data.frame(analysis_manifest),
  "04_3_analysis_manifest",
  ext = "csv"
)

# =========================================================
# 07 common enrichment helpers
# adapted directly from 04_1 / 04_2
# =========================================================
measured_genes <- rownames(expr)

bg_map <- tryCatch(
  clusterProfiler::bitr(
    measured_genes,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  ),
  error = function(e) data.frame()
)

bg_map <- as.data.table(bg_map)
bg_map <- unique(bg_map[!is.na(ENTREZID)])
bg_entrez <- unique(bg_map$ENTREZID)

save_dat_overwrite(
  as.data.frame(bg_map),
  "04_3_measured_genes_SYMBOL_to_ENTREZ",
  ext = "csv"
)

convert_symbols_to_entrez <- function(genes) {
  if (length(genes) == 0) {
    return(
      data.table(
        SYMBOL = character(),
        ENTREZID = character()
      )
    )
  }
  
  out <- tryCatch(
    clusterProfiler::bitr(
      genes,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db
    ),
    error = function(e) data.frame()
  )
  
  out <- as.data.table(out)
  unique(out[!is.na(ENTREZID)])
}

make_ranked_entrez <- function(tt, stat_col = "t") {
  rk <- tt[, .(gene, stat = get(stat_col))]
  rk <- rk[!is.na(stat)]
  
  map <- convert_symbols_to_entrez(rk$gene)
  rk <- merge(
    rk,
    map,
    by.x = "gene",
    by.y = "SYMBOL"
  )
  
  if (nrow(rk) == 0) return(numeric())
  
  rk <- rk[order(-abs(stat))]
  rk <- rk[!duplicated(ENTREZID)]
  
  geneList <- rk$stat
  names(geneList) <- rk$ENTREZID
  
  sort(
    geneList,
    decreasing = TRUE
  )
}

make_ranked_symbol <- function(tt, stat_col = "t") {
  rk <- tt[, .(gene, stat = get(stat_col))]
  rk <- rk[!is.na(stat)]
  rk <- rk[!duplicated(gene)]
  
  geneList <- rk$stat
  names(geneList) <- rk$gene
  
  sort(
    geneList,
    decreasing = TRUE
  )
}

run_ora_go <- function(
    genes_symbol,
    universe_entrez,
    ont = "BP"
) {
  gene_map <- convert_symbols_to_entrez(
    genes_symbol
  )
  
  if (nrow(gene_map) < 5) return(NULL)
  
  tryCatch(
    clusterProfiler::enrichGO(
      gene = unique(gene_map$ENTREZID),
      universe = universe_entrez,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = ont,
      pAdjustMethod = "BH",
      pvalueCutoff = 0.1,
      qvalueCutoff = 0.2,
      readable = TRUE
    ),
    error = function(e) NULL
  )
}

run_ora_kegg <- function(
    genes_symbol,
    universe_entrez
) {
  if (!run_online_kegg) return(NULL)
  
  gene_map <- convert_symbols_to_entrez(
    genes_symbol
  )
  
  if (nrow(gene_map) < 5) return(NULL)
  
  tryCatch(
    clusterProfiler::enrichKEGG(
      gene = unique(gene_map$ENTREZID),
      universe = universe_entrez,
      organism = "hsa",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.1,
      qvalueCutoff = 0.2
    ),
    error = function(e) NULL
  )
}

run_gsea_go <- function(
    geneList_entrez,
    ont = "BP"
) {
  if (length(geneList_entrez) < 10) return(NULL)
  
  tryCatch(
    clusterProfiler::gseGO(
      geneList = geneList_entrez,
      OrgDb = org.Hs.eg.db,
      ont = ont,
      keyType = "ENTREZID",
      minGSSize = 5,
      maxGSSize = 500,
      pvalueCutoff = 0.2,
      verbose = FALSE
    ),
    error = function(e) NULL
  )
}

run_gsea_kegg <- function(
    geneList_entrez
) {
  if (!run_online_kegg) return(NULL)
  if (length(geneList_entrez) < 10) return(NULL)
  
  tryCatch(
    clusterProfiler::gseKEGG(
      geneList = geneList_entrez,
      organism = "hsa",
      minGSSize = 5,
      maxGSSize = 500,
      pvalueCutoff = 0.2,
      verbose = FALSE
    ),
    error = function(e) NULL
  )
}

read_gmt_to_pathway_list <- function(gmt_file) {
  gmt_dt <- clusterProfiler::read.gmt(
    gmt_file
  )
  
  gmt_dt <- as.data.table(
    gmt_dt
  )
  
  setnames(
    gmt_dt,
    c("term", "gene")
  )
  
  pathways <- split(
    gmt_dt$gene,
    gmt_dt$term
  )
  
  lapply(
    pathways,
    unique
  )
}

run_fgsea_gmt <- function(
    geneList_symbol,
    pathways,
    minSize = 5,
    maxSize = 500
) {
  if (length(geneList_symbol) < 10) return(NULL)
  if (length(pathways) == 0) return(NULL)
  
  fg <- tryCatch(
    fgsea::fgseaMultilevel(
      pathways = pathways,
      stats = geneList_symbol,
      minSize = minSize,
      maxSize = maxSize
    ),
    error = function(e) NULL
  )
  
  if (
    is.null(fg) ||
    nrow(fg) == 0
  ) {
    return(NULL)
  }
  
  fg <- as.data.table(
    fg
  )
  
  if (
    "padj" %in% names(fg) &&
    !("p.adjust" %in% names(fg))
  ) {
    fg[, p.adjust := padj]
  }
  
  if (
    "pval" %in% names(fg) &&
    !("pvalue" %in% names(fg))
  ) {
    fg[, pvalue := pval]
  }
  
  fg[, absNES := abs(NES)]
  setorder(
    fg,
    p.adjust,
    -absNES
  )
  
  fg[, leadingEdge := vapply(
    leadingEdge,
    function(x) paste(x, collapse = "/"),
    FUN.VALUE = character(1)
  )]
  
  fg[, core_enrichment := leadingEdge]
  
  fg[, .(
    ID = pathway,
    Description = pathway,
    NES,
    pvalue = if ("pvalue" %in% names(fg)) {
      pvalue
    } else {
      NA_real_
    },
    p.adjust = if ("p.adjust" %in% names(fg)) {
      p.adjust
    } else {
      NA_real_
    },
    core_enrichment
  )]
}

gmt_pathways <- lapply(
  gmt_files,
  read_gmt_to_pathway_list
)

gmt_summary <- rbindlist(
  lapply(
    names(gmt_pathways),
    function(nm) {
      data.table(
        collection = nm,
        n_pathways = length(gmt_pathways[[nm]]),
        n_unique_genes = length(
          unique(
            unlist(
              gmt_pathways[[nm]]
            )
          )
        )
      )
    }
  )
)

save_dat_overwrite(
  as.data.frame(gmt_summary),
  "04_3_GMT_collection_summary",
  ext = "csv"
)

# =========================================================
# 08 harmonize enrichment result
# =========================================================
enrich_obj_to_dt <- function(
    obj,
    collection_name
) {
  if (is.null(obj)) return(data.table())
  
  df <- tryCatch(
    as.data.frame(obj),
    error = function(e) NULL
  )
  
  if (
    is.null(df) ||
    nrow(df) == 0
  ) {
    return(data.table())
  }
  
  x <- as.data.table(df)
  
  if (!"ID" %in% names(x)) return(data.table())
  
  data.table(
    collection = collection_name,
    pathway_id = as.character(x$ID),
    Description = if ("Description" %in% names(x)) {
      as.character(x$Description)
    } else {
      as.character(x$ID)
    },
    NES = if ("NES" %in% names(x)) {
      as.numeric(x$NES)
    } else {
      NA_real_
    },
    pvalue = if ("pvalue" %in% names(x)) {
      as.numeric(x$pvalue)
    } else {
      NA_real_
    },
    p.adjust = if ("p.adjust" %in% names(x)) {
      as.numeric(x$p.adjust)
    } else {
      NA_real_
    },
    core_enrichment = if ("core_enrichment" %in% names(x)) {
      as.character(x$core_enrichment)
    } else if ("geneID" %in% names(x)) {
      as.character(x$geneID)
    } else {
      NA_character_
    }
  )
}

fgsea_dt_to_dt <- function(
    x,
    collection_name
) {
  if (
    is.null(x) ||
    nrow(x) == 0
  ) {
    return(data.table())
  }
  
  x <- as.data.table(
    x
  )
  
  data.table(
    collection = collection_name,
    pathway_id = as.character(x$ID),
    Description = as.character(x$Description),
    NES = as.numeric(x$NES),
    pvalue = as.numeric(x$pvalue),
    p.adjust = as.numeric(x$p.adjust),
    core_enrichment = as.character(x$core_enrichment)
  )
}

# =========================================================
# 09 run enrichment for one ranked contrast
# =========================================================
run_gsea_all_collections <- function(tt) {
  out_list <- list()
  
  geneList_entrez <- make_ranked_entrez(
    tt,
    stat_col = "t"
  )
  
  geneList_symbol <- make_ranked_symbol(
    tt,
    stat_col = "t"
  )
  
  for (ont in c("BP", "CC", "MF")) {
    ggo <- run_gsea_go(
      geneList_entrez,
      ont = ont
    )
    
    if (!is.null(ggo)) {
      ggo <- tryCatch(
        setReadable(
          ggo,
          OrgDb = org.Hs.eg.db,
          keyType = "ENTREZID"
        ),
        error = function(e) ggo
      )
    }
    
    tmp <- enrich_obj_to_dt(
      ggo,
      paste0("GO_", ont)
    )
    
    if (nrow(tmp) > 0) {
      out_list[[length(out_list) + 1]] <- tmp
    }
  }
  
  if (run_online_kegg) {
    gkegg <- run_gsea_kegg(
      geneList_entrez
    )
    
    if (!is.null(gkegg)) {
      gkegg <- tryCatch(
        setReadable(
          gkegg,
          OrgDb = org.Hs.eg.db,
          keyType = "ENTREZID"
        ),
        error = function(e) gkegg
      )
    }
    
    tmp <- enrich_obj_to_dt(
      gkegg,
      "KEGG"
    )
    
    if (nrow(tmp) > 0) {
      out_list[[length(out_list) + 1]] <- tmp
    }
  }
  
  for (collection_nm in names(gmt_pathways)) {
    fg_dt <- run_fgsea_gmt(
      geneList_symbol,
      gmt_pathways[[collection_nm]],
      minSize = 5,
      maxSize = 500
    )
    
    tmp <- fgsea_dt_to_dt(
      fg_dt,
      collection_nm
    )
    
    if (nrow(tmp) > 0) {
      out_list[[length(out_list) + 1]] <- tmp
    }
  }
  
  if (length(out_list) == 0) {
    return(data.table())
  }
  
  rbindlist(
    out_list,
    fill = TRUE
  )
}

# =========================================================
# 10 run pairwise ORA for one contrast
# =========================================================
run_pairwise_ora <- function(tt) {
  result_list <- list()
  
  up_genes <- tt[
    adj.P.Val < padj_cutoff &
      logFC >= logfc_cutoff,
    gene
  ]
  
  down_genes <- tt[
    adj.P.Val < padj_cutoff &
      logFC <= -logfc_cutoff,
    gene
  ]
  
  for (direction_name in c(
    "up_first",
    "up_second"
  )) {
    genes_use <- if (direction_name == "up_first") {
      up_genes
    } else {
      down_genes
    }
    
    for (ont in c(
      "BP",
      "CC",
      "MF"
    )) {
      ego <- run_ora_go(
        genes_use,
        bg_entrez,
        ont = ont
      )
      
      tmp <- enrich_obj_to_dt(
        ego,
        paste0("GO_", ont)
      )
      
      if (nrow(tmp) > 0) {
        tmp[, direction := direction_name]
        result_list[[length(result_list) + 1]] <- tmp
      }
    }
    
    if (run_online_kegg) {
      ekegg <- run_ora_kegg(
        genes_use,
        bg_entrez
      )
      
      tmp <- enrich_obj_to_dt(
        ekegg,
        "KEGG"
      )
      
      if (nrow(tmp) > 0) {
        tmp[, direction := direction_name]
        result_list[[length(result_list) + 1]] <- tmp
      }
    }
  }
  
  if (length(result_list) == 0) {
    return(data.table())
  }
  
  rbindlist(
    result_list,
    fill = TRUE
  )
}

# =========================================================
# 11 helper: safe subtype score
# =========================================================
mean_available <- function(...) {
  x <- c(...)
  if (all(is.na(x))) {
    return(NA_real_)
  }
  mean(
    x,
    na.rm = TRUE
  )
}

# =========================================================
# 12 main cell-type loop
# =========================================================
all_pairwise_deg <- list()
all_vsrest_deg <- list()
all_pairwise_ora <- list()
all_pairwise_gsea <- list()
all_vsrest_gsea <- list()
all_qc <- list()
run_log <- list()

for (i in seq_len(nrow(analysis_manifest))) {
  annotation_level_i <- analysis_manifest$annotation_level[i]
  celltype_col_i <- analysis_manifest$celltype_col[i]
  celltype_i <- analysis_manifest$celltype[i]
  min_cells_i <- analysis_manifest$min_cells_per_fov[i]
  
  tag <- paste0(
    safe_tag(annotation_level_i),
    "__",
    safe_tag(celltype_i)
  )
  
  message("")
  message("=========================================================")
  message(
    "04_3 cell type ",
    i,
    "/",
    nrow(analysis_manifest),
    ": ",
    annotation_level_i,
    " / ",
    celltype_i
  )
  message("=========================================================")
  
  checkpoint_file <- file.path(
    checkpoint_dir,
    paste0(
      "04_3_",
      tag,
      "__complete.rds"
    )
  )
  
  if (
    resume_existing &&
    file.exists(checkpoint_file)
  ) {
    message(
      "Reusing completed checkpoint: ",
      checkpoint_file
    )
    
    result_i <- readRDS(
      checkpoint_file
    )
  } else {
    cell_dt <- meta_new[
      !is.na(get(celltype_col_i)) &
        get(celltype_col_i) == celltype_i &
        !is.na(FOV) &
        !is.na(Run_Tissue_name) &
        !is.na(immunotype_BW_top_40) &
        immunotype_BW_top_40 %in% immunotype_levels
    ]
    
    fov_counts <- cell_dt[
      ,
      .N,
      by = .(
        FOV,
        immunotype_BW_top_40,
        Run_Tissue_name
      )
    ][
      order(-N)
    ]
    
    valid_fov <- fov_counts[
      N >= min_cells_i,
      FOV
    ]
    
    cell_dt <- cell_dt[
      FOV %in% valid_fov
    ]
    
    fov_counts_valid <- fov_counts[
      FOV %in% valid_fov
    ]
    
    n_fov_group <- fov_counts_valid[
      ,
      .(
        n_FOV = uniqueN(FOV)
      ),
      by = immunotype_BW_top_40
    ]
    
    group_grid <- data.table(
      immunotype_BW_top_40 = immunotype_levels
    )
    
    n_fov_group <- merge(
      group_grid,
      n_fov_group,
      by = "immunotype_BW_top_40",
      all.x = TRUE
    )
    
    n_fov_group[
      is.na(n_FOV),
      n_FOV := 0L
    ]
    
    n_slide <- uniqueN(
      fov_counts_valid$Run_Tissue_name
    )
    
    qc_i <- data.table(
      annotation_level = annotation_level_i,
      celltype_col = celltype_col_i,
      celltype = celltype_i,
      min_cells_per_fov = min_cells_i,
      n_cells_total = nrow(cell_dt),
      n_FOV_total = uniqueN(cell_dt$FOV),
      n_slides = n_slide,
      n_FOV_ID = n_fov_group[
        immunotype_BW_top_40 == "ID",
        n_FOV
      ],
      n_FOV_IE = n_fov_group[
        immunotype_BW_top_40 == "IE",
        n_FOV
      ],
      n_FOV_IN = n_fov_group[
        immunotype_BW_top_40 == "IN",
        n_FOV
      ]
    )
    
    eligible <- all(
      n_fov_group$n_FOV >= min_fov_per_immunotype
    ) &&
      uniqueN(cell_dt$FOV) >= 6 &&
      n_slide >= 2
    
    qc_i[, eligible := eligible]
    
    if (!eligible) {
      qc_i[, status := "skipped_insufficient_FOV_or_slide"]
      
      result_i <- list(
        qc = qc_i,
        pairwise_deg = data.table(),
        vsrest_deg = data.table(),
        pairwise_ora = data.table(),
        pairwise_gsea = data.table(),
        vsrest_gsea = data.table()
      )
      
      saveRDS(
        result_i,
        checkpoint_file
      )
    } else {
      expr_cell <- expr[
        ,
        cell_dt$cell_id,
        drop = FALSE
      ]
      
      stopifnot(
        identical(
          colnames(expr_cell),
          cell_dt$cell_id
        )
      )
      
      fov_anno <- unique(
        cell_dt[
          ,
          .(
            FOV,
            immunotype_BW_top_40,
            Run_Tissue_name
          )
        ]
      )
      
      conflicting_fov <- fov_anno[
        ,
        .N,
        by = FOV
      ][
        N > 1
      ]
      
      if (nrow(conflicting_fov) > 0) {
        stop(
          "Some FOVs have conflicting annotations for cell type: ",
          celltype_i
        )
      }
      
      fov_list <- split(
        seq_len(nrow(cell_dt)),
        cell_dt$FOV
      )
      
      pb_list <- lapply(
        names(fov_list),
        function(fv) {
          idx <- fov_list[[fv]]
          
          cells_use <- cell_dt$cell_id[
            idx
          ]
          
          mat <- expr_cell[
            ,
            cells_use,
            drop = FALSE
          ]
          
          if (ncol(mat) == 1) {
            val <- as.numeric(
              mat[
                ,
                1
              ]
            )
          } else {
            val <- Matrix::rowMeans(
              mat
            )
          }
          
          data.table(
            FOV = fv,
            gene = rownames(expr_cell),
            expr_mean = val
          )
        }
      )
      
      pb_long <- rbindlist(
        pb_list
      )
      
      pb_wide <- dcast(
        pb_long,
        gene ~ FOV,
        value.var = "expr_mean"
      )
      
      pb_mat <- as.matrix(
        pb_wide[
          ,
          -1
        ]
      )
      
      rownames(pb_mat) <- pb_wide$gene
      
      fov_anno <- fov_anno[
        match(
          colnames(pb_mat),
          FOV
        )
      ]
      
      stopifnot(
        identical(
          fov_anno$FOV,
          colnames(pb_mat)
        )
      )
      
      fov_anno[
        ,
        immunotype_BW_top_40 := factor(
          immunotype_BW_top_40,
          levels = immunotype_levels
        )
      ]
      
      fov_anno[
        ,
        Run_Tissue_name := factor(
          Run_Tissue_name
        )
      ]
      
      design <- model.matrix(
        ~ 0 + immunotype_BW_top_40,
        data = fov_anno
      )
      
      colnames(design) <- make.names(
        colnames(design)
      )
      
      if (qr(design)$rank < ncol(design)) {
        stop(
          "Design matrix is not full rank for cell type: ",
          celltype_i
        )
      }
      
      corfit <- limma::duplicateCorrelation(
        pb_mat,
        design = design,
        block = fov_anno$Run_Tissue_name
      )
      
      fit <- limma::lmFit(
        pb_mat,
        design = design,
        block = fov_anno$Run_Tissue_name,
        correlation = corfit$consensus
      )
      
      # ===================================================
      # 12A true pairwise contrasts: 04_1 logic
      # ===================================================
      contrast_pairwise <- limma::makeContrasts(
        IE_vs_IN =
          immunotype_BW_top_40IE -
          immunotype_BW_top_40IN,
        
        ID_vs_IN =
          immunotype_BW_top_40ID -
          immunotype_BW_top_40IN,
        
        ID_vs_IE =
          immunotype_BW_top_40ID -
          immunotype_BW_top_40IE,
        
        levels = design
      )
      
      fit_pairwise <- limma::contrasts.fit(
        fit,
        contrast_pairwise
      )
      
      fit_pairwise <- limma::eBayes(
        fit_pairwise,
        trend = TRUE
      )
      
      pairwise_deg_list <- list()
      pairwise_ora_list <- list()
      pairwise_gsea_list <- list()
      
      for (cn in colnames(contrast_pairwise)) {
        tt <- limma::topTable(
          fit_pairwise,
          coef = cn,
          number = Inf,
          sort.by = "P"
        )
        
        tt$gene <- rownames(tt)
        tt <- as.data.table(tt)
        
        tt[
          ,
          is_sig :=
            adj.P.Val < padj_cutoff &
            abs(logFC) >= logfc_cutoff
        ]
        
        tt[
          ,
          direction := fifelse(
            logFC > 0,
            "up_in_first",
            "up_in_second"
          )
        ]
        
        tt[
          ,
          `:=`(
            annotation_level = annotation_level_i,
            celltype_col = celltype_col_i,
            celltype = celltype_i,
            contrast = cn,
            analysis_mode = "pairwise"
          )
        ]
        
        pairwise_deg_list[[length(pairwise_deg_list) + 1]] <- tt
        
        ora_dt <- run_pairwise_ora(
          tt
        )
        
        if (nrow(ora_dt) > 0) {
          ora_dt[
            ,
            `:=`(
              annotation_level = annotation_level_i,
              celltype_col = celltype_col_i,
              celltype = celltype_i,
              contrast = cn,
              analysis_mode = "pairwise_ORA"
            )
          ]
          
          pairwise_ora_list[[length(pairwise_ora_list) + 1]] <- ora_dt
        }
        
        gsea_dt <- run_gsea_all_collections(
          tt
        )
        
        if (nrow(gsea_dt) > 0) {
          gsea_dt[
            ,
            `:=`(
              annotation_level = annotation_level_i,
              celltype_col = celltype_col_i,
              celltype = celltype_i,
              contrast = cn,
              subtype = NA_character_,
              analysis_mode = "pairwise_GSEA"
            )
          ]
          
          pairwise_gsea_list[[length(pairwise_gsea_list) + 1]] <- gsea_dt
        }
      }
      
      pairwise_deg_i <- rbindlist(
        pairwise_deg_list,
        fill = TRUE
      )
      
      pairwise_ora_i <- if (length(pairwise_ora_list) > 0) {
        rbindlist(
          pairwise_ora_list,
          fill = TRUE
        )
      } else {
        data.table()
      }
      
      pairwise_gsea_i <- if (length(pairwise_gsea_list) > 0) {
        rbindlist(
          pairwise_gsea_list,
          fill = TRUE
        )
      } else {
        data.table()
      }
      
      # ===================================================
      # 12B true one-vs-rest contrasts: 04_2 logic
      # ===================================================
      contrast_vsrest <- limma::makeContrasts(
        ID_vs_rest =
          immunotype_BW_top_40ID -
          (
            immunotype_BW_top_40IE +
              immunotype_BW_top_40IN
          ) / 2,
        
        IE_vs_rest =
          immunotype_BW_top_40IE -
          (
            immunotype_BW_top_40ID +
              immunotype_BW_top_40IN
          ) / 2,
        
        IN_vs_rest =
          immunotype_BW_top_40IN -
          (
            immunotype_BW_top_40ID +
              immunotype_BW_top_40IE
          ) / 2,
        
        levels = design
      )
      
      fit_vsrest <- limma::contrasts.fit(
        fit,
        contrast_vsrest
      )
      
      fit_vsrest <- limma::eBayes(
        fit_vsrest,
        trend = TRUE
      )
      
      vsrest_deg_list <- list()
      vsrest_gsea_list <- list()
      
      for (cn in colnames(contrast_vsrest)) {
        tt <- limma::topTable(
          fit_vsrest,
          coef = cn,
          number = Inf,
          sort.by = "P"
        )
        
        tt$gene <- rownames(tt)
        tt <- as.data.table(tt)
        
        tt[
          ,
          is_sig :=
            adj.P.Val < padj_cutoff &
            abs(logFC) >= logfc_cutoff
        ]
        
        subtype_i <- sub(
          "_vs_rest$",
          "",
          cn
        )
        
        tt[
          ,
          `:=`(
            annotation_level = annotation_level_i,
            celltype_col = celltype_col_i,
            celltype = celltype_i,
            contrast = cn,
            subtype = subtype_i,
            analysis_mode = "one_vs_rest"
          )
        ]
        
        vsrest_deg_list[[length(vsrest_deg_list) + 1]] <- tt
        
        gsea_dt <- run_gsea_all_collections(
          tt
        )
        
        if (nrow(gsea_dt) > 0) {
          gsea_dt[
            ,
            `:=`(
              annotation_level = annotation_level_i,
              celltype_col = celltype_col_i,
              celltype = celltype_i,
              contrast = cn,
              subtype = subtype_i,
              analysis_mode = "one_vs_rest_GSEA"
            )
          ]
          
          vsrest_gsea_list[[length(vsrest_gsea_list) + 1]] <- gsea_dt
        }
      }
      
      vsrest_deg_i <- rbindlist(
        vsrest_deg_list,
        fill = TRUE
      )
      
      vsrest_gsea_i <- if (length(vsrest_gsea_list) > 0) {
        rbindlist(
          vsrest_gsea_list,
          fill = TRUE
        )
      } else {
        data.table()
      }
      
      qc_i[
        ,
        `:=`(
          status = "completed",
          consensus_correlation = corfit$consensus
        )
      ]
      
      result_i <- list(
        qc = qc_i,
        pairwise_deg = pairwise_deg_i,
        vsrest_deg = vsrest_deg_i,
        pairwise_ora = pairwise_ora_i,
        pairwise_gsea = pairwise_gsea_i,
        vsrest_gsea = vsrest_gsea_i
      )
      
      saveRDS(
        result_i,
        checkpoint_file
      )
    }
  }
  
  all_qc[[length(all_qc) + 1]] <- result_i$qc
  
  if (nrow(result_i$pairwise_deg) > 0) {
    all_pairwise_deg[[length(all_pairwise_deg) + 1]] <- result_i$pairwise_deg
  }
  
  if (nrow(result_i$vsrest_deg) > 0) {
    all_vsrest_deg[[length(all_vsrest_deg) + 1]] <- result_i$vsrest_deg
  }
  
  if (nrow(result_i$pairwise_ora) > 0) {
    all_pairwise_ora[[length(all_pairwise_ora) + 1]] <- result_i$pairwise_ora
  }
  
  if (nrow(result_i$pairwise_gsea) > 0) {
    all_pairwise_gsea[[length(all_pairwise_gsea) + 1]] <- result_i$pairwise_gsea
  }
  
  if (nrow(result_i$vsrest_gsea) > 0) {
    all_vsrest_gsea[[length(all_vsrest_gsea) + 1]] <- result_i$vsrest_gsea
  }
  
  run_log[[length(run_log) + 1]] <- result_i$qc[
    ,
    .(
      annotation_level,
      celltype,
      status,
      eligible,
      n_FOV_ID,
      n_FOV_IE,
      n_FOV_IN,
      n_slides
    )
  ]
  
  gc(
    verbose = FALSE
  )
}

# =========================================================
# 13 combine master results
# =========================================================
qc_master <- rbindlist(
  all_qc,
  fill = TRUE
)

pairwise_deg_master <- if (length(all_pairwise_deg) > 0) {
  rbindlist(
    all_pairwise_deg,
    fill = TRUE
  )
} else {
  data.table()
}

vsrest_deg_master <- if (length(all_vsrest_deg) > 0) {
  rbindlist(
    all_vsrest_deg,
    fill = TRUE
  )
} else {
  data.table()
}

pairwise_ora_master <- if (length(all_pairwise_ora) > 0) {
  rbindlist(
    all_pairwise_ora,
    fill = TRUE
  )
} else {
  data.table()
}

pairwise_gsea_master <- if (length(all_pairwise_gsea) > 0) {
  rbindlist(
    all_pairwise_gsea,
    fill = TRUE
  )
} else {
  data.table()
}

vsrest_gsea_master <- if (length(all_vsrest_gsea) > 0) {
  rbindlist(
    all_vsrest_gsea,
    fill = TRUE
  )
} else {
  data.table()
}

run_log_master <- rbindlist(
  run_log,
  fill = TRUE
)

save_dat_overwrite(
  as.data.frame(qc_master),
  "04_3_celltype_QC_summary",
  ext = "csv"
)

save_dat_overwrite(
  as.data.frame(run_log_master),
  "04_3_run_log",
  ext = "csv"
)

if (nrow(pairwise_deg_master) > 0) {
  save_dat_overwrite(
    as.data.frame(pairwise_deg_master),
    "04_3_PAIRWISE_DEG_ALL_celltypes_ALL_genes",
    ext = "csv"
  )
  
  save_dat_overwrite(
    pairwise_deg_master,
    "04_3_PAIRWISE_DEG_ALL_celltypes_ALL_genes",
    ext = "rds"
  )
}

if (nrow(vsrest_deg_master) > 0) {
  save_dat_overwrite(
    as.data.frame(vsrest_deg_master),
    "04_3_vsREST_DEG_ALL_celltypes_ALL_genes",
    ext = "csv"
  )
  
  save_dat_overwrite(
    vsrest_deg_master,
    "04_3_vsREST_DEG_ALL_celltypes_ALL_genes",
    ext = "rds"
  )
}

if (nrow(pairwise_ora_master) > 0) {
  save_dat_overwrite(
    as.data.frame(pairwise_ora_master),
    "04_3_PAIRWISE_ORA_complete_long_table",
    ext = "csv"
  )
  
  save_dat_overwrite(
    pairwise_ora_master,
    "04_3_PAIRWISE_ORA_complete_long_table",
    ext = "rds"
  )
}

if (nrow(pairwise_gsea_master) > 0) {
  save_dat_overwrite(
    as.data.frame(pairwise_gsea_master),
    "04_3_PAIRWISE_GSEA_complete_long_table",
    ext = "csv"
  )
  
  save_dat_overwrite(
    pairwise_gsea_master,
    "04_3_PAIRWISE_GSEA_complete_long_table",
    ext = "rds"
  )
}

if (nrow(vsrest_gsea_master) > 0) {
  save_dat_overwrite(
    as.data.frame(vsrest_gsea_master),
    "04_3_vsREST_GSEA_complete_long_table",
    ext = "csv"
  )
  
  save_dat_overwrite(
    vsrest_gsea_master,
    "04_3_vsREST_GSEA_complete_long_table",
    ext = "rds"
  )
}

# =========================================================
# 14 pairwise integrated GSEA table
# same 3-column score logic as 04_1
# =========================================================
if (nrow(pairwise_gsea_master) > 0) {
  pairwise_wide <- dcast(
    pairwise_gsea_master,
    annotation_level +
      celltype_col +
      celltype +
      collection +
      pathway_id +
      Description ~ contrast,
    value.var = c(
      "NES",
      "p.adjust",
      "core_enrichment"
    )
  )
  
  required_nes <- c(
    "NES_IE_vs_IN",
    "NES_ID_vs_IN",
    "NES_ID_vs_IE"
  )
  
  required_p <- c(
    "p.adjust_IE_vs_IN",
    "p.adjust_ID_vs_IN",
    "p.adjust_ID_vs_IE"
  )
  
  required_core <- c(
    "core_enrichment_IE_vs_IN",
    "core_enrichment_ID_vs_IN",
    "core_enrichment_ID_vs_IE"
  )
  
  for (nm in required_nes) {
    if (!nm %in% names(pairwise_wide)) {
      pairwise_wide[, (nm) := NA_real_]
    }
  }
  
  for (nm in required_p) {
    if (!nm %in% names(pairwise_wide)) {
      pairwise_wide[, (nm) := NA_real_]
    }
  }
  
  for (nm in required_core) {
    if (!nm %in% names(pairwise_wide)) {
      pairwise_wide[, (nm) := NA_character_]
    }
  }
  
  pairwise_wide[
    ,
    score_ID := mapply(
      mean_available,
      NES_ID_vs_IN,
      NES_ID_vs_IE
    )
  ]
  
  pairwise_wide[
    ,
    score_IE := mapply(
      mean_available,
      NES_IE_vs_IN,
      -NES_ID_vs_IE
    )
  ]
  
  pairwise_wide[
    ,
    score_IN := mapply(
      mean_available,
      -NES_IE_vs_IN,
      -NES_ID_vs_IN
    )
  ]
  
  score_mat <- as.matrix(
    pairwise_wide[
      ,
      .(
        score_ID,
        score_IE,
        score_IN
      )
    ]
  )
  
  pairwise_wide[
    ,
    best_subtype_by_score := NA_character_
  ]
  
  valid_idx <- rowSums(
    !is.na(
      score_mat
    )
  ) > 0
  
  if (any(valid_idx)) {
    tmp <- score_mat[
      valid_idx,
      ,
      drop = FALSE
    ]
    
    tmp[
      is.na(tmp)
    ] <- -Inf
    
    pairwise_wide[
      valid_idx,
      best_subtype_by_score :=
        c(
          "ID",
          "IE",
          "IN"
        )[
          max.col(
            tmp,
            ties.method = "first"
          )
        ]
    ]
  }
  
  save_dat_overwrite(
    as.data.frame(pairwise_wide),
    "04_3_PAIRWISE_GSEA_complete_integrated_table",
    ext = "csv"
  )
  
  save_dat_overwrite(
    pairwise_wide,
    "04_3_PAIRWISE_GSEA_complete_integrated_table",
    ext = "rds"
  )
}

# =========================================================
# 15 true one-vs-rest integrated GSEA table
# same logic as 04_2
# =========================================================
if (nrow(vsrest_gsea_master) > 0) {
  vsrest_wide <- dcast(
    vsrest_gsea_master,
    annotation_level +
      celltype_col +
      celltype +
      collection +
      pathway_id +
      Description ~ subtype,
    value.var = c(
      "NES",
      "p.adjust",
      "core_enrichment"
    )
  )
  
  for (subtype_name in c(
    "ID",
    "IE",
    "IN"
  )) {
    nm_nes <- paste0(
      "NES_",
      subtype_name
    )
    
    nm_p <- paste0(
      "p.adjust_",
      subtype_name
    )
    
    nm_core <- paste0(
      "core_enrichment_",
      subtype_name
    )
    
    if (!nm_nes %in% names(vsrest_wide)) {
      vsrest_wide[, (nm_nes) := NA_real_]
    }
    
    if (!nm_p %in% names(vsrest_wide)) {
      vsrest_wide[, (nm_p) := NA_real_]
    }
    
    if (!nm_core %in% names(vsrest_wide)) {
      vsrest_wide[, (nm_core) := NA_character_]
    }
  }
  
  score_mat <- as.matrix(
    vsrest_wide[
      ,
      .(
        NES_ID,
        NES_IE,
        NES_IN
      )
    ]
  )
  
  vsrest_wide[
    ,
    best_subtype_by_score := NA_character_
  ]
  
  valid_idx <- rowSums(
    !is.na(
      score_mat
    )
  ) > 0
  
  if (any(valid_idx)) {
    tmp <- score_mat[
      valid_idx,
      ,
      drop = FALSE
    ]
    
    tmp[
      is.na(tmp)
    ] <- -Inf
    
    vsrest_wide[
      valid_idx,
      best_subtype_by_score :=
        c(
          "ID",
          "IE",
          "IN"
        )[
          max.col(
            tmp,
            ties.method = "first"
          )
        ]
    ]
  }
  
  save_dat_overwrite(
    as.data.frame(vsrest_wide),
    "04_3_vsREST_GSEA_complete_integrated_table",
    ext = "csv"
  )
  
  save_dat_overwrite(
    vsrest_wide,
    "04_3_vsREST_GSEA_complete_integrated_table",
    ext = "rds"
  )
}

# =========================================================
# 16 universal GSEA query table
# =========================================================
gsea_query_master <- rbindlist(
  list(
    pairwise_gsea_master,
    vsrest_gsea_master
  ),
  fill = TRUE
)

if (nrow(gsea_query_master) > 0) {
  save_dat_overwrite(
    as.data.frame(gsea_query_master),
    "04_3_GSEA_QUERY_MASTER_long",
    ext = "csv"
  )
  
  save_dat_overwrite(
    gsea_query_master,
    "04_3_GSEA_QUERY_MASTER_long",
    ext = "rds"
  )
}

# =========================================================
# 17 convenience significant tables
# =========================================================
if (nrow(pairwise_gsea_master) > 0) {
  tmp <- pairwise_gsea_master[
    !is.na(p.adjust) &
      p.adjust < gsea_padj_query_cutoff
  ]
  
  tmp[, absNES_order := abs(NES)]
  
  setorder(
    tmp,
    annotation_level,
    celltype,
    contrast,
    collection,
    p.adjust,
    -absNES_order
  )
  
  tmp[, absNES_order := NULL]
  
  save_dat_overwrite(
    as.data.frame(tmp),
    "04_3_PAIRWISE_GSEA_padj_lt_0_20",
    ext = "csv"
  )
}

if (nrow(vsrest_gsea_master) > 0) {
  tmp <- vsrest_gsea_master[
    !is.na(p.adjust) &
      p.adjust < gsea_padj_query_cutoff
  ]
  
  tmp[, absNES_order := abs(NES)]
  
  setorder(
    tmp,
    annotation_level,
    celltype,
    subtype,
    collection,
    p.adjust,
    -absNES_order
  )
  
  tmp[, absNES_order := NULL]
  
  save_dat_overwrite(
    as.data.frame(tmp),
    "04_3_vsREST_GSEA_padj_lt_0_20",
    ext = "csv"
  )
}

if (nrow(pairwise_ora_master) > 0) {
  tmp <- pairwise_ora_master[
    !is.na(p.adjust) &
      p.adjust < 0.05
  ]
  
  setorder(
    tmp,
    annotation_level,
    celltype,
    contrast,
    direction,
    collection,
    p.adjust
  )
  
  save_dat_overwrite(
    as.data.frame(tmp),
    "04_3_PAIRWISE_ORA_padj_lt_0_05",
    ext = "csv"
  )
}

# =========================================================
# 18 parameter record
# =========================================================
parameter_dt <- data.table(
  parameter = c(
    "formal_input",
    "include_tumor_in_04_3",
    "include_cd8_subtypes",
    "min_cells_per_fov_broad",
    "min_cells_per_fov_fine_cd8",
    "min_cells_per_fov_tumor",
    "min_fov_per_immunotype",
    "padj_cutoff",
    "logfc_cutoff",
    "run_online_kegg",
    "GMT_collections",
    "blocking_method",
    "pairwise_method",
    "vsrest_method",
    "outer_parallelization"
  ),
  value = c(
    "03_2 metadata + TNBC_6k_norm_Bogang.RDS; FOV pseudobulk recomputed within each cell type",
    as.character(include_tumor_in_04_3),
    as.character(include_cd8_subtypes),
    as.character(min_cells_per_fov_broad),
    as.character(min_cells_per_fov_fine_cd8),
    as.character(min_cells_per_fov_tumor),
    as.character(min_fov_per_immunotype),
    as.character(padj_cutoff),
    as.character(logfc_cutoff),
    as.character(run_online_kegg),
    "H,C2,C5,C7",
    "limma duplicateCorrelation(block = Run_Tissue_name)",
    "true limma contrasts: IE_vs_IN, ID_vs_IN, ID_vs_IE",
    "true limma contrasts: ID_vs_rest, IE_vs_rest, IN_vs_rest",
    "sequential; data.table threads = 8"
  )
)

save_dat_overwrite(
  as.data.frame(parameter_dt),
  "04_3_parameters",
  ext = "csv"
)

message("")
message("04_3 finished successfully.")
message("")
message("Formal analysis:")
message("  cell-type-specific FOV pseudobulk")
message("  duplicateCorrelation(block = Run_Tissue_name)")
message("  true 04_1 pairwise contrasts")
message("  true 04_2 one-vs-rest contrasts")
message("")
message("Primary query files:")
message("  04_3_GSEA_QUERY_MASTER_long.csv")
message("  04_3_PAIRWISE_GSEA_complete_integrated_table.csv")
message("  04_3_vsREST_GSEA_complete_integrated_table.csv")
message("  04_3_PAIRWISE_ORA_complete_long_table.csv")
message("  04_3_celltype_QC_summary.csv")
message("")
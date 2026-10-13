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

out_dir <- file.path(proj, "results", "04_2")
fig_dir <- file.path(out_dir, "figures")
dat_dir <- file.path(out_dir, "data")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dat_dir, recursive = TRUE, showWarnings = FALSE)

save_dat_overwrite <- function(x, filename, ext = c("csv", "rds")) {
  ext <- match.arg(ext)
  out <- file.path(dat_dir, paste0(filename, ".", ext))
  
  if (file.exists(out)) {
    ok_rm <- tryCatch(file.remove(out), warning = function(e) FALSE, error = function(e) FALSE)
    if (!isTRUE(ok_rm)) {
      stop("Cannot overwrite existing file: ", out)
    }
  }
  
  if (ext == "csv") {
    tryCatch(
      write.csv(x, out, row.names = FALSE),
      error = function(e) stop("Failed to write csv: ", out, "\n", conditionMessage(e))
    )
  }
  
  if (ext == "rds") {
    tryCatch(
      saveRDS(x, out),
      error = function(e) stop("Failed to write rds: ", out, "\n", conditionMessage(e))
    )
  }
  
  message("Saved data: ", out)
  invisible(out)
}

save_plot_pdf <- function(plot_obj, filename, w = 7, h = 6) {
  out <- file.path(fig_dir, paste0(filename, ".pdf"))
  
  if (file.exists(out)) {
    ok_rm <- tryCatch(file.remove(out), warning = function(e) FALSE, error = function(e) FALSE)
    if (!isTRUE(ok_rm)) {
      stop("Cannot overwrite existing pdf: ", out)
    }
  }
  
  tryCatch(
    ggplot2::ggsave(out, plot = plot_obj, width = w, height = h, device = cairo_pdf),
    error = function(e) stop("Failed to write pdf: ", out, "\n", conditionMessage(e))
  )
  
  message("Saved figure: ", out)
  invisible(out)
}

# =========================================================
# 02 packages
# =========================================================
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(limma)
  library(ggplot2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(fgsea)
})

# =========================================================
# 03 load data
# =========================================================
norm <- readRDS(file.path(raw, "TNBC_6k_norm_Bogang.RDS"))

meta_new_file <- file.path(proj, "results", "03_2", "data", "03_2_meta_with_immunotype_BW_top_40.rds")
if (!file.exists(meta_new_file)) {
  stop("Cannot find 03_2 output: ", meta_new_file)
}
meta_new <- readRDS(meta_new_file)
meta_new <- as.data.table(meta_new)

gmt_files <- list(
  H  = file.path(gene_set_dir, "h.all.v2025.1.Hs.symbols.gmt"),
  C2 = file.path(gene_set_dir, "c2.all.v2025.1.Hs.symbols.gmt"),
  C5 = file.path(gene_set_dir, "c5.all.v2025.1.Hs.symbols.gmt"),
  C7 = file.path(gene_set_dir, "c7.all.v2025.1.Hs.symbols.gmt")
)
for (nm in names(gmt_files)) {
  if (!file.exists(gmt_files[[nm]])) stop("Missing GMT file: ", gmt_files[[nm]])
}

# =========================================================
# 04 checks
# =========================================================
need_cols <- c("cell_id", "FOV", "BC6k_refined", "Run_Tissue_name", "immunotype_BW_top_40")
stopifnot(all(need_cols %in% names(meta_new)))

# =========================================================
# 05 extract expression matrix
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

if (!any(meta_new$cell_id %in% colnames(expr)) && any(meta_new$cell_id %in% rownames(expr))) {
  expr <- t(expr)
}

common_cells <- intersect(meta_new$cell_id, colnames(expr))
common_cells <- common_cells[match(common_cells, meta_new$cell_id)]
common_cells <- common_cells[!is.na(common_cells)]

expr <- expr[, common_cells, drop = FALSE]
meta_new <- meta_new[match(common_cells, meta_new$cell_id)]

stopifnot(identical(colnames(expr), meta_new$cell_id))

if (anyDuplicated(rownames(expr)) > 0) {
  expr_dt <- data.table(gene = rownames(expr), as.matrix(expr))
  expr_dt <- expr_dt[, lapply(.SD, mean), by = gene]
  expr <- as.matrix(expr_dt[, -1])
  rownames(expr) <- expr_dt$gene
}

# =========================================================
# 06 tumor-only pseudobulk by FOV
# =========================================================
tumor_dt <- meta_new[
  BC6k_refined == "tumor" &
    !is.na(FOV) &
    !is.na(immunotype_BW_top_40) &
    immunotype_BW_top_40 %in% c("ID", "IE", "IN") &
    !is.na(Run_Tissue_name)
]

min_tumor_cells <- 30
tumor_fov_counts <- tumor_dt[, .N, by = .(FOV, immunotype_BW_top_40, Run_Tissue_name)][order(-N)]
save_dat_overwrite(as.data.frame(tumor_fov_counts), "04_2_tumor_cell_counts_per_FOV", ext = "csv")

valid_fov <- tumor_fov_counts[N >= min_tumor_cells, FOV]
tumor_dt <- tumor_dt[FOV %in% valid_fov]

if (nrow(tumor_dt) == 0) stop("No FOV passed min_tumor_cells threshold.")

expr_tumor <- expr[, tumor_dt$cell_id, drop = FALSE]
stopifnot(identical(colnames(expr_tumor), tumor_dt$cell_id))

fov_anno <- unique(tumor_dt[, .(FOV, immunotype_BW_top_40, Run_Tissue_name)])
chk1 <- fov_anno[, .N, by = FOV][N > 1]
if (nrow(chk1) > 0) stop("Some FOVs have conflicting annotations.")

fov_list <- split(seq_len(nrow(tumor_dt)), tumor_dt$FOV)

pb_list <- lapply(names(fov_list), function(fv) {
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

pb_long <- rbindlist(pb_list)
save_dat_overwrite(as.data.frame(pb_long), "04_2_tumor_pseudobulk_long", ext = "csv")

pb_wide <- dcast(pb_long, gene ~ FOV, value.var = "expr_mean")
pb_mat <- as.matrix(pb_wide[, -1])
rownames(pb_mat) <- pb_wide$gene

fov_anno <- fov_anno[match(colnames(pb_mat), FOV)]
stopifnot(identical(fov_anno$FOV, colnames(pb_mat)))

save_dat_overwrite(as.data.frame(fov_anno), "04_2_FOV_annotation_for_GSEA", ext = "csv")

slide_by_group <- fov_anno[, .N, by = .(Run_Tissue_name, immunotype_BW_top_40)]
save_dat_overwrite(as.data.frame(slide_by_group), "04_2_Run_Tissue_by_immunotype_table", ext = "csv")

# =========================================================
# 07 one-vs-rest limma model with duplicateCorrelation
# =========================================================
fov_anno[, immunotype_BW_top_40 := factor(immunotype_BW_top_40, levels = c("ID", "IE", "IN"))]
fov_anno[, Run_Tissue_name := factor(Run_Tissue_name)]

design <- model.matrix(~ 0 + immunotype_BW_top_40, data = fov_anno)
colnames(design) <- make.names(colnames(design))

corfit <- limma::duplicateCorrelation(
  pb_mat,
  design = design,
  block = fov_anno$Run_Tissue_name
)

save_dat_overwrite(
  data.frame(consensus.correlation = corfit$consensus),
  "04_2_duplicateCorrelation_summary",
  ext = "csv"
)

fit <- limma::lmFit(
  pb_mat,
  design = design,
  block = fov_anno$Run_Tissue_name,
  correlation = corfit$consensus
)

contrast_mat <- limma::makeContrasts(
  ID_vs_rest = immunotype_BW_top_40ID - (immunotype_BW_top_40IE + immunotype_BW_top_40IN) / 2,
  IE_vs_rest = immunotype_BW_top_40IE - (immunotype_BW_top_40ID + immunotype_BW_top_40IN) / 2,
  IN_vs_rest = immunotype_BW_top_40IN - (immunotype_BW_top_40ID + immunotype_BW_top_40IE) / 2,
  levels = design
)

fit2 <- limma::contrasts.fit(fit, contrast_mat)
fit2 <- limma::eBayes(fit2, trend = TRUE)

save_dat_overwrite(as.data.frame(contrast_mat), "04_2_contrast_matrix", ext = "csv")

contrast_names <- colnames(contrast_mat)

# optional DEG export
padj_cutoff <- 0.05
logfc_cutoff <- 0.25

deg_summary_list <- list()
deg_tables <- list()

for (cn in contrast_names) {
  tt <- limma::topTable(fit2, coef = cn, number = Inf, sort.by = "P")
  tt$gene <- rownames(tt)
  tt <- as.data.table(tt)
  tt[, is_sig := adj.P.Val < padj_cutoff & abs(logFC) >= logfc_cutoff]
  deg_tables[[cn]] <- tt
  
  save_dat_overwrite(as.data.frame(tt), paste0("04_2_DEG_", cn), ext = "csv")
  
  deg_summary_list[[cn]] <- data.table(
    contrast = cn,
    n_sig = sum(tt$is_sig, na.rm = TRUE),
    n_up = sum(tt$adj.P.Val < padj_cutoff & tt$logFC >  logfc_cutoff, na.rm = TRUE),
    n_down = sum(tt$adj.P.Val < padj_cutoff & tt$logFC < -logfc_cutoff, na.rm = TRUE)
  )
}
save_dat_overwrite(as.data.frame(rbindlist(deg_summary_list)), "04_2_DEG_summary", ext = "csv")

# =========================================================
# 08 helpers for GSEA
# =========================================================
convert_symbols_to_entrez <- function(genes) {
  if (length(genes) == 0) return(data.table(SYMBOL = character(), ENTREZID = character()))
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
  out <- unique(out[!is.na(ENTREZID)])
  out
}

make_ranked_entrez <- function(tt, stat_col = "t") {
  rk <- tt[, .(gene, stat = get(stat_col))]
  rk <- rk[!is.na(stat)]
  map <- convert_symbols_to_entrez(rk$gene)
  rk <- merge(rk, map, by.x = "gene", by.y = "SYMBOL")
  if (nrow(rk) == 0) return(numeric())
  rk <- rk[order(-abs(stat))]
  rk <- rk[!duplicated(ENTREZID)]
  geneList <- rk$stat
  names(geneList) <- rk$ENTREZID
  sort(geneList, decreasing = TRUE)
}

make_ranked_symbol <- function(tt, stat_col = "t") {
  rk <- tt[, .(gene, stat = get(stat_col))]
  rk <- rk[!is.na(stat)]
  rk <- rk[!duplicated(gene)]
  geneList <- rk$stat
  names(geneList) <- rk$gene
  sort(geneList, decreasing = TRUE)
}

run_gsea_go <- function(geneList_entrez, ont = "BP") {
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

run_gsea_kegg <- function(geneList_entrez) {
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
  gmt_dt <- clusterProfiler::read.gmt(gmt_file)
  gmt_dt <- as.data.table(gmt_dt)
  setnames(gmt_dt, c("term", "gene"))
  pathways <- split(gmt_dt$gene, gmt_dt$term)
  pathways <- lapply(pathways, unique)
  pathways
}

run_fgsea_gmt <- function(geneList_symbol, pathways, minSize = 5, maxSize = 500) {
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
  if (is.null(fg) || nrow(fg) == 0) return(NULL)
  
  fg <- as.data.table(fg)
  if ("padj" %in% names(fg) && !("p.adjust" %in% names(fg))) fg[, p.adjust := padj]
  if ("pval" %in% names(fg) && !("pvalue" %in% names(fg))) fg[, pvalue := pval]
  
  fg[, absNES := abs(NES)]
  setorder(fg, p.adjust, -absNES)
  
  fg[, leadingEdge := vapply(leadingEdge, function(x) paste(x, collapse = "/"), FUN.VALUE = character(1))]
  fg[, core_enrichment := leadingEdge]
  
  fg[, .(
    ID = pathway,
    Description = pathway,
    NES,
    pvalue = if ("pvalue" %in% names(fg)) pvalue else NA_real_,
    p.adjust = if ("p.adjust" %in% names(fg)) p.adjust else NA_real_,
    core_enrichment
  )]
}

save_enrich_results <- function(obj, prefix) {
  if (is.null(obj)) return(invisible(NULL))
  
  df <- tryCatch(as.data.frame(obj), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) {
    message("No enrichment result for: ", prefix)
    return(invisible(NULL))
  }
  
  save_dat_overwrite(df, prefix, ext = "csv")
  save_dat_overwrite(df, prefix, ext = "rds")
  
  p1 <- tryCatch(
    clusterProfiler::dotplot(obj, showCategory = min(15, nrow(df))) +
      ggplot2::ggtitle(prefix),
    error = function(e) NULL
  )
  if (!is.null(p1)) save_plot_pdf(p1, prefix, w = 8, h = 6)
  
  invisible(df)
}

save_fgsea_results <- function(dt, prefix) {
  if (is.null(dt) || nrow(dt) == 0) {
    message("No fgsea result for: ", prefix)
    return(invisible(NULL))
  }
  
  save_dat_overwrite(as.data.frame(dt), prefix, ext = "csv")
  save_dat_overwrite(as.data.frame(dt), prefix, ext = "rds")
  
  plot_dt <- copy(dt)[p.adjust < 0.2]
  if (nrow(plot_dt) == 0) return(invisible(NULL))
  
  plot_dt <- plot_dt[order(NES)]
  if (nrow(plot_dt) > 20) {
    idx <- unique(c(seq_len(min(10, nrow(plot_dt))), pmax(nrow(plot_dt) - 9, 1):nrow(plot_dt)))
    plot_dt <- plot_dt[idx]
  }
  
  plot_dt[, Description := factor(Description, levels = Description)]
  
  p <- ggplot(plot_dt, aes(x = NES, y = Description, color = p.adjust)) +
    geom_point(size = 2.2) +
    scale_color_gradient(low = "#B2182B", high = "#2166AC", trans = "reverse") +
    labs(
      title = prefix,
      x = "NES",
      y = NULL,
      color = "adj.P"
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 10),
      axis.text.y = element_text(size = 7)
    )
  
  save_plot_pdf(p, prefix, w = 7.2, h = max(4.2, 0.22 * nrow(plot_dt)))
  invisible(dt)
}

# =========================================================
# 09 prepare GMT pathways
# =========================================================
gmt_pathways <- lapply(gmt_files, read_gmt_to_pathway_list)

gmt_sizes <- rbindlist(lapply(names(gmt_pathways), function(nm) {
  data.table(
    collection = nm,
    n_pathways = length(gmt_pathways[[nm]]),
    n_unique_genes = length(unique(unlist(gmt_pathways[[nm]])))
  )
}))
save_dat_overwrite(as.data.frame(gmt_sizes), "04_2_GMT_collection_summary", ext = "csv")

# =========================================================
# 10 run one-vs-rest GSEA
# =========================================================
gsea_quick_list <- list()

for (cn in contrast_names) {
  tt <- deg_tables[[cn]]
  
  geneList_entrez <- make_ranked_entrez(tt, stat_col = "t")
  geneList_symbol <- make_ranked_symbol(tt, stat_col = "t")
  
  # GO
  for (go_ont in c("BP", "CC", "MF")) {
    ggo <- run_gsea_go(geneList_entrez, ont = go_ont)
    if (!is.null(ggo)) {
      ggo <- tryCatch(
        setReadable(ggo, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"),
        error = function(e) ggo
      )
    }
    save_enrich_results(ggo, paste0("04_2_GSEA_GO_", go_ont, "_", cn))
  }
  
  # KEGG
  gkegg <- run_gsea_kegg(geneList_entrez)
  if (!is.null(gkegg)) {
    gkegg <- tryCatch(
      setReadable(gkegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"),
      error = function(e) gkegg
    )
  }
  save_enrich_results(gkegg, paste0("04_2_GSEA_KEGG_", cn))
  
  # GMT collections
  for (collection_nm in names(gmt_pathways)) {
    fg_dt <- run_fgsea_gmt(geneList_symbol, gmt_pathways[[collection_nm]], minSize = 5, maxSize = 500)
    save_fgsea_results(fg_dt, paste0("04_2_GSEA_", collection_nm, "_", cn))
  }
  
  gsea_quick_list[[cn]] <- data.table(contrast = cn)
}
save_dat_overwrite(as.data.frame(rbindlist(gsea_quick_list, fill = TRUE)), "04_2_GSEA_run_summary", ext = "csv")

# =========================================================
# 11 build complete integrated GSEA table
# =========================================================
collections <- c("GO_BP", "GO_CC", "GO_MF", "KEGG", "H", "C2", "C5", "C7")
contrast_order <- c("ID_vs_rest", "IE_vs_rest", "IN_vs_rest")

gsea_file_map <- list(
  GO_BP = c("04_2_GSEA_GO_BP_ID_vs_rest.rds", "04_2_GSEA_GO_BP_IE_vs_rest.rds", "04_2_GSEA_GO_BP_IN_vs_rest.rds"),
  GO_CC = c("04_2_GSEA_GO_CC_ID_vs_rest.rds", "04_2_GSEA_GO_CC_IE_vs_rest.rds", "04_2_GSEA_GO_CC_IN_vs_rest.rds"),
  GO_MF = c("04_2_GSEA_GO_MF_ID_vs_rest.rds", "04_2_GSEA_GO_MF_IE_vs_rest.rds", "04_2_GSEA_GO_MF_IN_vs_rest.rds"),
  KEGG  = c("04_2_GSEA_KEGG_ID_vs_rest.rds",  "04_2_GSEA_KEGG_IE_vs_rest.rds",  "04_2_GSEA_KEGG_IN_vs_rest.rds"),
  H     = c("04_2_GSEA_H_ID_vs_rest.rds",     "04_2_GSEA_H_IE_vs_rest.rds",     "04_2_GSEA_H_IN_vs_rest.rds"),
  C2    = c("04_2_GSEA_C2_ID_vs_rest.rds",    "04_2_GSEA_C2_IE_vs_rest.rds",    "04_2_GSEA_C2_IN_vs_rest.rds"),
  C5    = c("04_2_GSEA_C5_ID_vs_rest.rds",    "04_2_GSEA_C5_IE_vs_rest.rds",    "04_2_GSEA_C5_IN_vs_rest.rds"),
  C7    = c("04_2_GSEA_C7_ID_vs_rest.rds",    "04_2_GSEA_C7_IE_vs_rest.rds",    "04_2_GSEA_C7_IN_vs_rest.rds")
)

read_gsea_result <- function(filepath, collection_name, contrast_name) {
  if (!file.exists(filepath)) return(NULL)
  
  obj <- tryCatch(readRDS(filepath), error = function(e) NULL)
  if (is.null(obj)) return(NULL)
  
  df <- tryCatch(as.data.frame(obj), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) return(NULL)
  
  dt <- as.data.table(df)
  
  if (!"ID" %in% names(dt) && "pathway" %in% names(dt)) dt[, ID := pathway]
  if (!"Description" %in% names(dt) && "pathway" %in% names(dt)) dt[, Description := pathway]
  if (!"p.adjust" %in% names(dt) && "padj" %in% names(dt)) dt[, p.adjust := padj]
  if (!"core_enrichment" %in% names(dt) && "leadingEdge" %in% names(dt)) dt[, core_enrichment := leadingEdge]
  
  need_cols <- c("ID", "Description", "NES", "p.adjust")
  miss <- setdiff(need_cols, names(dt))
  if (length(miss) > 0) return(NULL)
  
  dt <- dt[, .(
    collection = collection_name,
    subtype = gsub("_vs_rest", "", contrast_name),
    pathway_id = as.character(ID),
    Description = as.character(Description),
    NES = as.numeric(NES),
    p.adjust = as.numeric(p.adjust),
    core_enrichment = if ("core_enrichment" %in% names(dt)) as.character(core_enrichment) else NA_character_
  )]
  dt
}

gsea_list <- list()
for (set_nm in names(gsea_file_map)) {
  file_vec <- gsea_file_map[[set_nm]]
  for (i in seq_along(file_vec)) {
    fp <- file.path(dat_dir, file_vec[i])
    cn <- contrast_order[i]
    gsea_list[[paste(set_nm, cn, sep = "__")]] <- read_gsea_result(fp, set_nm, cn)
  }
}

gsea_long <- rbindlist(gsea_list, fill = TRUE)
if (nrow(gsea_long) == 0) stop("No GSEA results found for section 11.")

gsea_long[, subtype := factor(subtype, levels = c("ID", "IE", "IN"))]
save_dat_overwrite(as.data.frame(gsea_long), "04_2_GSEA_complete_long_table", ext = "csv")
save_dat_overwrite(as.data.frame(gsea_long), "04_2_GSEA_complete_long_table", ext = "rds")

gsea_wide <- dcast(
  gsea_long,
  collection + pathway_id + Description ~ subtype,
  value.var = c("NES", "p.adjust", "core_enrichment")
)

gsea_wide[, best_subtype_by_score := c("ID", "IE", "IN")[max.col(.SD, ties.method = "first")],
          .SDcols = c("NES_ID", "NES_IE", "NES_IN")]

gsea_wide[, sig_ID := !is.na(p.adjust_ID) & p.adjust_ID < 0.20 & NES_ID > 0]
gsea_wide[, sig_IE := !is.na(p.adjust_IE) & p.adjust_IE < 0.20 & NES_IE > 0]
gsea_wide[, sig_IN := !is.na(p.adjust_IN) & p.adjust_IN < 0.20 & NES_IN > 0]

gsea_wide[, n_sig_subtypes := rowSums(cbind(sig_ID, sig_IE, sig_IN), na.rm = TRUE)]

gsea_wide[, assigned_subtype := NA_character_]
gsea_wide[sig_ID == TRUE & (is.na(assigned_subtype)), assigned_subtype := "ID"]
gsea_wide[sig_IE == TRUE & (is.na(assigned_subtype)), assigned_subtype := "IE"]
gsea_wide[sig_IN == TRUE & (is.na(assigned_subtype)), assigned_subtype := "IN"]

# if >1 subtype significant, assign by largest NES
multi_idx <- gsea_wide[n_sig_subtypes > 1, which = TRUE]
if (length(multi_idx) > 0) {
  tmp <- as.matrix(gsea_wide[multi_idx, .(NES_ID, NES_IE, NES_IN)])
  tmp[is.na(tmp)] <- -Inf
  best <- apply(tmp, 1, function(x) c("ID", "IE", "IN")[which.max(x)])
  gsea_wide[multi_idx, assigned_subtype := best]
}

gsea_complete <- gsea_wide[, .(
  collection,
  pathway_id,
  Description,
  assigned_subtype,
  best_subtype_by_score,
  n_sig_subtypes,
  NES_ID,
  p.adjust_ID,
  core_enrichment_ID,
  NES_IE,
  p.adjust_IE,
  core_enrichment_IE,
  NES_IN,
  p.adjust_IN,
  core_enrichment_IN
)]

save_dat_overwrite(as.data.frame(gsea_complete), "04_2_GSEA_complete_integrated_table", ext = "csv")
save_dat_overwrite(as.data.frame(gsea_complete), "04_2_GSEA_complete_integrated_table", ext = "rds")

# =========================================================
# 12 automatic one-vs-rest GSEA heatmap
# =========================================================
n_show <- 15

plot_pool <- copy(gsea_complete)[!is.na(assigned_subtype)]
if (nrow(plot_pool) > 0) {
  plot_pool[, rank_score := fifelse(
    assigned_subtype == "ID", NES_ID,
    fifelse(assigned_subtype == "IE", NES_IE, NES_IN)
  )]
  
  plot_pool <- plot_pool[order(assigned_subtype, -rank_score)]
  plot_top <- plot_pool[, head(.SD, n_show), by = assigned_subtype]
  
  plot_top[, row_label := paste0("[", assigned_subtype, "] ", collection, " | ", Description)]
  
  row_order_dt <- unique(
    plot_top[, .(row_label, assigned_subtype, rank_score)][
      order(factor(assigned_subtype, levels = c("ID", "IE", "IN")), rank_score)
    ]
  )
  row_order <- rev(row_order_dt$row_label)
  
  heat_dt <- rbindlist(list(
    plot_top[, .(row_label, subtype = "ID", value = NES_ID)],
    plot_top[, .(row_label, subtype = "IE", value = NES_IE)],
    plot_top[, .(row_label, subtype = "IN", value = NES_IN)]
  ), fill = TRUE)
  
  heat_dt[, subtype := factor(subtype, levels = c("ID", "IE", "IN"))]
  heat_dt[, row_label := factor(row_label, levels = row_order)]
  
  p_main <- ggplot(heat_dt, aes(x = subtype, y = row_label, fill = value)) +
    geom_tile(color = "white", linewidth = 0.35) +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      name = "NES"
    ) +
    labs(
      x = NULL,
      y = NULL,
      title = "Subtype-specific one-vs-rest GSEA pathways",
      subtitle = paste0("Top ", n_show, " pathways shown per subtype")
    ) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(size = 10, face = "bold"),
      axis.text.y = element_text(size = 7),
      plot.title = element_text(face = "bold", size = 11),
      plot.subtitle = element_text(size = 9),
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8)
    )
  
  save_plot_pdf(
    p_main,
    "04_2_GSEA_subtype_specific_unified_heatmap",
    w = 7.2,
    h = max(4.5, 0.20 * length(unique(heat_dt$row_label)))
  )
}

# =========================================================
# 13 custom selected pathway heatmap by pathway_id
# =========================================================
# fill pathway_id and simplified display_name manually before running this section
selected_pathways <- data.table(
  pathway_id = c(
    # "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    # "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"
  ),
  display_name = c(
    # "IFNG response",
    # "EMT"
  )
)

if (nrow(selected_pathways) > 0) {
  selected_pathways[, input_order := seq_len(.N)]
  
  plot_dt <- merge(
    selected_pathways,
    gsea_complete,
    by = "pathway_id",
    all.x = TRUE
  )
  
  plot_dt <- plot_dt[!is.na(collection)]
  if (nrow(plot_dt) > 0) {
    plot_dt[is.na(assigned_subtype), assigned_subtype := best_subtype_by_score]
    
    dup_name <- plot_dt[, .N, by = display_name][N > 1, display_name]
    plot_dt[, row_label := ifelse(
      display_name %in% dup_name,
      paste0(display_name, " [", collection, "]"),
      display_name
    )]
    
    row_order_dt <- unique(
      plot_dt[, .(row_label, assigned_subtype, input_order)][
        order(factor(assigned_subtype, levels = c("ID", "IE", "IN")), input_order)
      ]
    )
    row_order <- rev(row_order_dt$row_label)
    
    heat_custom <- rbindlist(list(
      plot_dt[, .(row_label, subtype = "ID", value = NES_ID)],
      plot_dt[, .(row_label, subtype = "IE", value = NES_IE)],
      plot_dt[, .(row_label, subtype = "IN", value = NES_IN)]
    ), fill = TRUE)
    
    heat_custom[, subtype := factor(subtype, levels = c("ID", "IE", "IN"))]
    heat_custom[, row_label := factor(row_label, levels = row_order)]
    
    out_selected <- unique(plot_dt[, .(
      pathway_id,
      display_name,
      collection,
      Description,
      assigned_subtype,
      best_subtype_by_score,
      NES_ID,
      p.adjust_ID,
      core_enrichment_ID,
      NES_IE,
      p.adjust_IE,
      core_enrichment_IE,
      NES_IN,
      p.adjust_IN,
      core_enrichment_IN
    )])
    
    save_dat_overwrite(as.data.frame(out_selected), "04_2_GSEA_custom_selected_pathways_table", ext = "csv")
    save_dat_overwrite(as.data.frame(out_selected), "04_2_GSEA_custom_selected_pathways_table", ext = "rds")
    
    p_custom <- ggplot(heat_custom, aes(x = subtype, y = row_label, fill = value)) +
      geom_tile(color = "white", linewidth = 0.35) +
      scale_fill_gradient2(
        low = "#2166AC",
        mid = "white",
        high = "#B2182B",
        midpoint = 0,
        name = "NES"
      ) +
      labs(
        x = NULL,
        y = NULL,
        title = "Selected GSEA pathways",
        subtitle = "Custom pathway selection by pathway_id"
      ) +
      theme_bw(base_size = 10) +
      theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 10, face = "bold"),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 9),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8)
      )
    
    save_plot_pdf(
      p_custom,
      "04_2_GSEA_custom_selected_pathways_heatmap",
      w = 6.6,
      h = max(4.2, 0.22 * length(unique(heat_custom$row_label)))
    )
  }
}

message("04_2 one-vs-rest GSEA pipeline finished successfully.")
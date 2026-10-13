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

out_dir <- file.path(proj, "results", "04_1")
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

save_plot_pdf <- function(plot_obj, filename, w = 7, h = 6) {
  out <- file.path(fig_dir, paste0(filename, ".pdf"))
  ggplot2::ggsave(out, plot = plot_obj, width = w, height = h, device = cairo_pdf)
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

# GMT files
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
    immunotype_BW_top_40 %in% c("IN", "IE", "ID") &
    !is.na(Run_Tissue_name)
]

min_tumor_cells <- 30
tumor_fov_counts <- tumor_dt[, .N, by = .(FOV, immunotype_BW_top_40, Run_Tissue_name)][order(-N)]
save_dat_overwrite(as.data.frame(tumor_fov_counts), "04_1_tumor_cell_counts_per_FOV", ext = "csv")

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
save_dat_overwrite(as.data.frame(pb_long), "04_1_tumor_pseudobulk_long", ext = "csv")

pb_wide <- dcast(pb_long, gene ~ FOV, value.var = "expr_mean")
pb_mat <- as.matrix(pb_wide[, -1])
rownames(pb_mat) <- pb_wide$gene

fov_anno <- fov_anno[match(colnames(pb_mat), FOV)]
stopifnot(identical(fov_anno$FOV, colnames(pb_mat)))

save_dat_overwrite(as.data.frame(fov_anno), "04_1_FOV_annotation_for_DE", ext = "csv")

slide_by_group <- fov_anno[, .N, by = .(Run_Tissue_name, immunotype_BW_top_40)]
save_dat_overwrite(as.data.frame(slide_by_group), "04_1_Run_Tissue_by_immunotype_table", ext = "csv")

# =========================================================
# 07 DE analysis with duplicateCorrelation
# =========================================================
fov_anno[, immunotype_BW_top_40 := factor(immunotype_BW_top_40, levels = c("IN", "IE", "ID"))]
fov_anno[, Run_Tissue_name := factor(Run_Tissue_name)]

design <- model.matrix(~ 0 + immunotype_BW_top_40, data = fov_anno)
colnames(design) <- make.names(colnames(design))

corfit <- limma::duplicateCorrelation(
  pb_mat,
  design = design,
  block = fov_anno$Run_Tissue_name
)

cor_dt <- data.table(consensus.correlation = corfit$consensus)
save_dat_overwrite(as.data.frame(cor_dt), "04_1_duplicateCorrelation_summary", ext = "csv")

fit <- limma::lmFit(
  pb_mat,
  design = design,
  block = fov_anno$Run_Tissue_name,
  correlation = corfit$consensus
)

contrast_mat <- limma::makeContrasts(
  IE_vs_IN = immunotype_BW_top_40IE - immunotype_BW_top_40IN,
  ID_vs_IN = immunotype_BW_top_40ID - immunotype_BW_top_40IN,
  ID_vs_IE = immunotype_BW_top_40ID - immunotype_BW_top_40IE,
  levels = design
)

fit2 <- limma::contrasts.fit(fit, contrast_mat)
fit2 <- limma::eBayes(fit2, trend = TRUE)

save_dat_overwrite(as.data.frame(contrast_mat), "04_1_contrast_matrix", ext = "csv")

padj_cutoff <- 0.05
logfc_cutoff <- 0.25
contrast_names <- colnames(contrast_mat)

deg_tables <- list()
deg_summary_list <- list()

for (cn in contrast_names) {
  tt <- limma::topTable(fit2, coef = cn, number = Inf, sort.by = "P")
  tt$gene <- rownames(tt)
  tt <- as.data.table(tt)
  
  tt[, is_sig := adj.P.Val < padj_cutoff & abs(logFC) >= logfc_cutoff]
  tt[, direction := fifelse(logFC > 0, "up_in_first", "up_in_second")]
  deg_tables[[cn]] <- tt
  
  save_dat_overwrite(as.data.frame(tt), paste0("04_1_DEG_", cn), ext = "csv")
  save_dat_overwrite(as.data.frame(tt), paste0("04_1_DEG_", cn), ext = "rds")
  
  up_genes <- tt[adj.P.Val < padj_cutoff & logFC >= logfc_cutoff, gene]
  down_genes <- tt[adj.P.Val < padj_cutoff & logFC <= -logfc_cutoff, gene]
  
  save_dat_overwrite(data.frame(gene = up_genes),   paste0("04_1_DEG_", cn, "_up_genes"), ext = "csv")
  save_dat_overwrite(data.frame(gene = down_genes), paste0("04_1_DEG_", cn, "_down_genes"), ext = "csv")
  
  deg_summary_list[[cn]] <- data.table(
    contrast = cn,
    n_sig = sum(tt$is_sig, na.rm = TRUE),
    n_up_first = sum(tt$adj.P.Val < padj_cutoff & tt$logFC >= logfc_cutoff, na.rm = TRUE),
    n_up_second = sum(tt$adj.P.Val < padj_cutoff & tt$logFC <= -logfc_cutoff, na.rm = TRUE)
  )
  
  tt[, neglog10padj := -log10(adj.P.Val + 1e-300)]
  p_volcano <- ggplot(tt, aes(x = logFC, y = neglog10padj)) +
    geom_point(alpha = 0.6, size = 1) +
    geom_vline(xintercept = c(-logfc_cutoff, logfc_cutoff), linetype = 2, color = "grey50") +
    geom_hline(yintercept = -log10(padj_cutoff), linetype = 2, color = "grey50") +
    labs(
      title = paste0("Volcano: ", cn, " (duplicateCorrelation)"),
      x = "logFC",
      y = "-log10(adj.P.Val)"
    ) +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(face = "bold"))
  save_plot_pdf(p_volcano, paste0("04_1_volcano_", cn), w = 6, h = 5)
}

deg_summary <- rbindlist(deg_summary_list, fill = TRUE)
save_dat_overwrite(as.data.frame(deg_summary), "04_1_DEG_summary", ext = "csv")

# =========================================================
# 08 ID conversion helpers
# =========================================================
measured_genes <- rownames(pb_mat)

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
  geneList <- sort(geneList, decreasing = TRUE)
  geneList
}

make_ranked_symbol <- function(tt, stat_col = "t") {
  rk <- tt[, .(gene, stat = get(stat_col))]
  rk <- rk[!is.na(stat)]
  rk <- rk[!duplicated(gene)]
  geneList <- rk$stat
  names(geneList) <- rk$gene
  geneList <- sort(geneList, decreasing = TRUE)
  geneList
}

# =========================================================
# 09 enrichment functions
# =========================================================
run_ora_go <- function(genes_symbol, universe_entrez, ont = "BP") {
  gene_map <- convert_symbols_to_entrez(genes_symbol)
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

run_ora_kegg <- function(genes_symbol, universe_entrez) {
  gene_map <- convert_symbols_to_entrez(genes_symbol)
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
  
  # fgsea 返回列通常包括 pathway, pval, padj, ES, NES, size, leadingEdge
  # 这里统一命名，避免后面混乱
  if ("padj" %in% names(fg) && !("p.adjust" %in% names(fg))) {
    fg[, p.adjust := padj]
  }
  if ("pval" %in% names(fg) && !("pvalue" %in% names(fg))) {
    fg[, pvalue := pval]
  }
  
  fg[, absNES := abs(NES)]
  setorder(fg, p.adjust, -absNES)
  
  fg[, leadingEdge := vapply(leadingEdge, function(x) paste(x, collapse = "/"), FUN.VALUE = character(1))]
  fg[, core_enrichment := leadingEdge]
  
  out <- fg[, .(
    ID = pathway,
    Description = pathway,
    setSize = if ("size" %in% names(fg)) size else NA_integer_,
    enrichmentScore = if ("ES" %in% names(fg)) ES else NA_real_,
    NES,
    pvalue = if ("pvalue" %in% names(fg)) pvalue else NA_real_,
    p.adjust = if ("p.adjust" %in% names(fg)) p.adjust else NA_real_,
    padj = if ("p.adjust" %in% names(fg)) p.adjust else NA_real_,
    size = if ("size" %in% names(fg)) size else NA_integer_,
    leadingEdge,
    core_enrichment
  )]
  
  out
}

save_enrich_results <- function(obj, prefix) {
  if (is.null(obj)) return(invisible(NULL))
  df <- tryCatch(as.data.frame(obj), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) return(invisible(NULL))
  
  save_dat_overwrite(df, prefix, ext = "csv")
  save_dat_overwrite(df, prefix, ext = "rds")
  
  p1 <- tryCatch(
    clusterProfiler::dotplot(obj, showCategory = min(15, nrow(df))) +
      ggplot2::ggtitle(prefix),
    error = function(e) NULL
  )
  if (!is.null(p1)) save_plot_pdf(p1, prefix, w = 8, h = 6)
}

save_fgsea_results <- function(dt, prefix) {
  if (is.null(dt) || nrow(dt) == 0) return(invisible(NULL))
  
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
    geom_point(size = 2.5) +
    scale_color_gradient(low = "#B2182B", high = "#2166AC", trans = "reverse") +
    labs(
      title = prefix,
      x = "NES",
      y = NULL,
      color = "adj.P"
    ) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))
  save_plot_pdf(p, prefix, w = 8, h = max(4.5, 0.22 * nrow(plot_dt)))
}

# =========================================================
# 10 prepare GMT pathway lists
# =========================================================
gmt_pathways <- lapply(gmt_files, read_gmt_to_pathway_list)

gmt_sizes <- rbindlist(lapply(names(gmt_pathways), function(nm) {
  data.table(
    collection = nm,
    n_pathways = length(gmt_pathways[[nm]]),
    n_unique_genes = length(unique(unlist(gmt_pathways[[nm]])))
  )
}))
save_dat_overwrite(as.data.frame(gmt_sizes), "04_1_GMT_collection_summary", ext = "csv")

# =========================================================
# 11 run GO / KEGG / GSEA / GMT-GSEA
# =========================================================
for (cn in contrast_names) {
  tt <- deg_tables[[cn]]
  
  up_genes <- tt[adj.P.Val < padj_cutoff & logFC >= logfc_cutoff, gene]
  down_genes <- tt[adj.P.Val < padj_cutoff & logFC <= -logfc_cutoff, gene]
  
  # ---------- ORA GO: BP / CC / MF ----------
  for (go_ont in c("BP", "CC", "MF")) {
    ego_up <- run_ora_go(up_genes, bg_entrez, ont = go_ont)
    ego_down <- run_ora_go(down_genes, bg_entrez, ont = go_ont)
    
    save_enrich_results(ego_up,   paste0("04_1_GO_", go_ont, "_", cn, "_up"))
    save_enrich_results(ego_down, paste0("04_1_GO_", go_ont, "_", cn, "_down"))
  }
  
  # ---------- ORA KEGG ----------
  ekegg_up <- run_ora_kegg(up_genes, bg_entrez)
  ekegg_down <- run_ora_kegg(down_genes, bg_entrez)
  
  if (!is.null(ekegg_up)) {
    ekegg_up <- tryCatch(
      setReadable(ekegg_up, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"),
      error = function(e) ekegg_up
    )
  }
  if (!is.null(ekegg_down)) {
    ekegg_down <- tryCatch(
      setReadable(ekegg_down, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"),
      error = function(e) ekegg_down
    )
  }
  
  save_enrich_results(ekegg_up,   paste0("04_1_KEGG_", cn, "_up"))
  save_enrich_results(ekegg_down, paste0("04_1_KEGG_", cn, "_down"))
  
  # ---------- clusterProfiler GSEA: GO / KEGG ----------
  geneList_entrez <- make_ranked_entrez(tt, stat_col = "t")
  
  gsea_go_quick <- list()
  
  for (go_ont in c("BP", "CC", "MF")) {
    ggo <- run_gsea_go(geneList_entrez, ont = go_ont)
    if (!is.null(ggo)) {
      ggo <- tryCatch(
        setReadable(ggo, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"),
        error = function(e) ggo
      )
    }
    save_enrich_results(ggo, paste0("04_1_GSEA_GO_", go_ont, "_", cn))
    
    top_go <- tryCatch(as.data.frame(ggo), error = function(e) data.frame())
    gsea_go_quick[[go_ont]] <- if (nrow(top_go) > 0) paste(head(top_go$Description, 10), collapse = " | ") else NA_character_
  }
  
  gkegg <- run_gsea_kegg(geneList_entrez)
  if (!is.null(gkegg)) {
    gkegg <- tryCatch(
      setReadable(gkegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"),
      error = function(e) gkegg
    )
  }
  save_enrich_results(gkegg, paste0("04_1_GSEA_KEGG_", cn))
  
  top_kegg <- tryCatch(as.data.frame(gkegg), error = function(e) data.frame())
  
  # ---------- GMT-based GSEA: H / C2 / C5 / C7 ----------
  geneList_symbol <- make_ranked_symbol(tt, stat_col = "t")
  
  gmt_quick <- list()
  for (collection_nm in names(gmt_pathways)) {
    fg_dt <- run_fgsea_gmt(geneList_symbol, gmt_pathways[[collection_nm]], minSize = 5, maxSize = 500)
    save_fgsea_results(fg_dt, paste0("04_1_GSEA_", collection_nm, "_", cn))
    
    top_dt <- if (!is.null(fg_dt)) fg_dt[p.adjust < 0.2] else data.table()
    gmt_quick[[collection_nm]] <- if (nrow(top_dt) > 0) paste(head(top_dt$Description, 10), collapse = " | ") else NA_character_
  }
  
  quick_dt <- data.table(
    contrast = cn,
    top_GSEA_GO_BP_terms = gsea_go_quick[["BP"]],
    top_GSEA_GO_CC_terms = gsea_go_quick[["CC"]],
    top_GSEA_GO_MF_terms = gsea_go_quick[["MF"]],
    top_GSEA_KEGG_terms = if (nrow(top_kegg) > 0) paste(head(top_kegg$Description, 10), collapse = " | ") else NA_character_,
    top_GSEA_H_terms = gmt_quick[["H"]],
    top_GSEA_C2_terms = gmt_quick[["C2"]],
    top_GSEA_C5_terms = gmt_quick[["C5"]],
    top_GSEA_C7_terms = gmt_quick[["C7"]]
  )
  save_dat_overwrite(as.data.frame(quick_dt), paste0("04_1_enrichment_quick_summary_", cn), ext = "csv")
}


# =========================================================
# 12 unified subtype-specific GSEA heatmap after updated 04_1
# collections included:
# GO_BP / GO_CC / GO_MF / KEGG / H / C2 / C5 / C7
# =========================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ---------- parameters ----------
gsea_padj_cutoff <- 0.20
n_show <- 15
subtype_order <- c("ID", "IE", "IN")
contrast_order <- c("IE_vs_IN", "ID_vs_IN", "ID_vs_IE")

# ---------- GSEA result files generated by updated 04_1 ----------
gsea_sets <- list(
  GO_BP = c("04_1_GSEA_GO_BP_IE_vs_IN.rds", "04_1_GSEA_GO_BP_ID_vs_IN.rds", "04_1_GSEA_GO_BP_ID_vs_IE.rds"),
  GO_CC = c("04_1_GSEA_GO_CC_IE_vs_IN.rds", "04_1_GSEA_GO_CC_ID_vs_IN.rds", "04_1_GSEA_GO_CC_ID_vs_IE.rds"),
  GO_MF = c("04_1_GSEA_GO_MF_IE_vs_IN.rds", "04_1_GSEA_GO_MF_ID_vs_IN.rds", "04_1_GSEA_GO_MF_ID_vs_IE.rds"),
  KEGG  = c("04_1_GSEA_KEGG_IE_vs_IN.rds",  "04_1_GSEA_KEGG_ID_vs_IN.rds",  "04_1_GSEA_KEGG_ID_vs_IE.rds"),
  H     = c("04_1_GSEA_H_IE_vs_IN.rds",     "04_1_GSEA_H_ID_vs_IN.rds",     "04_1_GSEA_H_ID_vs_IE.rds"),
  C2    = c("04_1_GSEA_C2_IE_vs_IN.rds",    "04_1_GSEA_C2_ID_vs_IN.rds",    "04_1_GSEA_C2_ID_vs_IE.rds"),
  C5    = c("04_1_GSEA_C5_IE_vs_IN.rds",    "04_1_GSEA_C5_ID_vs_IN.rds",    "04_1_GSEA_C5_ID_vs_IE.rds"),
  C7    = c("04_1_GSEA_C7_IE_vs_IN.rds",    "04_1_GSEA_C7_ID_vs_IN.rds",    "04_1_GSEA_C7_ID_vs_IE.rds")
)

# ---------- helper: read one GSEA result ----------
read_gsea_result <- function(filepath, collection_name, contrast_name) {
  if (!file.exists(filepath)) return(NULL)
  
  obj <- tryCatch(readRDS(filepath), error = function(e) NULL)
  if (is.null(obj)) return(NULL)
  
  df <- tryCatch(as.data.frame(obj), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) return(NULL)
  
  dt <- as.data.table(df)
  
  # harmonize columns between clusterProfiler GSEA and fgsea GMT-GSEA
  if (!"ID" %in% names(dt) && "pathway" %in% names(dt)) dt[, ID := pathway]
  if (!"Description" %in% names(dt) && "pathway" %in% names(dt)) dt[, Description := pathway]
  if (!"p.adjust" %in% names(dt) && "padj" %in% names(dt)) dt[, p.adjust := padj]
  if (!"core_enrichment" %in% names(dt) && "leadingEdge" %in% names(dt)) dt[, core_enrichment := leadingEdge]
  
  need_cols <- c("ID", "Description", "NES", "p.adjust")
  miss <- setdiff(need_cols, names(dt))
  if (length(miss) > 0) return(NULL)
  
  dt <- dt[, .(
    ID = as.character(ID),
    Description = as.character(Description),
    NES = as.numeric(NES),
    p.adjust = as.numeric(p.adjust),
    core_enrichment = if ("core_enrichment" %in% names(dt)) as.character(core_enrichment) else NA_character_
  )]
  
  dt[, collection := collection_name]
  dt[, contrast := contrast_name]
  dt
}

# ---------- load all GSEA results into one unified pool ----------
gsea_list <- list()
for (set_nm in names(gsea_sets)) {
  file_vec <- gsea_sets[[set_nm]]
  for (i in seq_along(file_vec)) {
    fp <- file.path(dat_dir, file_vec[i])
    cn <- contrast_order[i]
    gsea_list[[paste(set_nm, cn, sep = "__")]] <- read_gsea_result(fp, set_nm, cn)
  }
}

gsea_dt <- rbindlist(gsea_list, fill = TRUE)
if (nrow(gsea_dt) == 0) stop("No GSEA result files found for section 11.")

# ---------- reshape wide ----------
gsea_wide <- dcast(
  gsea_dt,
  collection + ID + Description ~ contrast,
  value.var = c("NES", "p.adjust", "core_enrichment")
)

# ---------- subtype-level score reconstruction ----------
# contrast definitions:
# IE_vs_IN = IE - IN
# ID_vs_IN = ID - IN
# ID_vs_IE = ID - IE

gsea_wide[, score_ID := rowMeans(cbind(NES_ID_vs_IN, NES_ID_vs_IE), na.rm = TRUE)]
gsea_wide[, score_IE := rowMeans(cbind(NES_IE_vs_IN, -NES_ID_vs_IE), na.rm = TRUE)]
gsea_wide[, score_IN := rowMeans(cbind(-NES_IE_vs_IN, -NES_ID_vs_IN), na.rm = TRUE)]

# conservative rank score
gsea_wide[, rank_ID := pmin(NES_ID_vs_IN, NES_ID_vs_IE, na.rm = TRUE)]
gsea_wide[, rank_IE := pmin(NES_IE_vs_IN, -NES_ID_vs_IE, na.rm = TRUE)]
gsea_wide[, rank_IN := pmin(-NES_IE_vs_IN, -NES_ID_vs_IN, na.rm = TRUE)]

# ---------- define subtype-specific pathways ----------
# ID-high: ID > IN and ID > IE
gsea_wide[, sig_ID :=
            !is.na(p.adjust_ID_vs_IN) & !is.na(p.adjust_ID_vs_IE) &
            p.adjust_ID_vs_IN < gsea_padj_cutoff &
            p.adjust_ID_vs_IE < gsea_padj_cutoff &
            NES_ID_vs_IN > 0 &
            NES_ID_vs_IE > 0
]

# IE-high: IE > IN and IE > ID
gsea_wide[, sig_IE :=
            !is.na(p.adjust_IE_vs_IN) & !is.na(p.adjust_ID_vs_IE) &
            p.adjust_IE_vs_IN < gsea_padj_cutoff &
            p.adjust_ID_vs_IE < gsea_padj_cutoff &
            NES_IE_vs_IN > 0 &
            NES_ID_vs_IE < 0
]

# IN-high: IN > IE and IN > ID
gsea_wide[, sig_IN :=
            !is.na(p.adjust_IE_vs_IN) & !is.na(p.adjust_ID_vs_IN) &
            p.adjust_IE_vs_IN < gsea_padj_cutoff &
            p.adjust_ID_vs_IN < gsea_padj_cutoff &
            NES_IE_vs_IN < 0 &
            NES_ID_vs_IN < 0
]

gsea_wide[, assigned_subtype := NA_character_]
gsea_wide[sig_ID == TRUE, assigned_subtype := "ID"]
gsea_wide[sig_IE == TRUE, assigned_subtype := "IE"]
gsea_wide[sig_IN == TRUE, assigned_subtype := "IN"]

# resolve occasional overlap by highest subtype score
overlap_idx <- gsea_wide[
  rowSums(cbind(sig_ID, sig_IE, sig_IN), na.rm = TRUE) > 1,
  which = TRUE
]
if (length(overlap_idx) > 0) {
  tmp <- gsea_wide[overlap_idx, .(score_ID, score_IE, score_IN)]
  best <- apply(as.matrix(tmp), 1, function(x) c("ID", "IE", "IN")[which.max(x)])
  gsea_wide[overlap_idx, assigned_subtype := best]
}

gsea_sel <- gsea_wide[!is.na(assigned_subtype)]
if (nrow(gsea_sel) == 0) stop("No subtype-specific GSEA pathways found under current cutoff.")

gsea_sel[, subtype_rank := fifelse(
  assigned_subtype == "ID", rank_ID,
  fifelse(assigned_subtype == "IE", rank_IE, rank_IN)
)]

# unified top-n selection across ALL collections together
gsea_sel <- gsea_sel[order(assigned_subtype, -subtype_rank)]
gsea_sel_top <- gsea_sel[, head(.SD, n_show), by = assigned_subtype]

# ---------- export unified master table ----------
out_master <- gsea_sel_top[, .(
  collection,
  pathway_id = ID,
  Description,
  assigned_subtype,
  score_ID,
  score_IE,
  score_IN,
  rank_ID,
  rank_IE,
  rank_IN,
  NES_IE_vs_IN,
  NES_ID_vs_IN,
  NES_ID_vs_IE,
  p.adjust_IE_vs_IN,
  p.adjust_ID_vs_IN,
  p.adjust_ID_vs_IE,
  core_enrichment_IE_vs_IN,
  core_enrichment_ID_vs_IN,
  core_enrichment_ID_vs_IE
)]

save_dat_overwrite(as.data.frame(out_master), "04_1_GSEA_subtype_specific_pathways_unified_table", ext = "csv")
save_dat_overwrite(as.data.frame(out_master), "04_1_GSEA_subtype_specific_pathways_unified_table", ext = "rds")

# ---------- build unified heatmap ----------
gsea_sel_top[, row_label := paste0("[", assigned_subtype, "] ", collection, " | ", Description)]

row_order_dt <- unique(
  gsea_sel_top[, .(row_label, assigned_subtype, subtype_rank)][
    order(factor(assigned_subtype, levels = subtype_order), subtype_rank)
  ]
)
row_order <- rev(row_order_dt$row_label)

heat_dt <- rbindlist(list(
  gsea_sel_top[, .(row_label, subtype = "ID", value = score_ID)],
  gsea_sel_top[, .(row_label, subtype = "IE", value = score_IE)],
  gsea_sel_top[, .(row_label, subtype = "IN", value = score_IN)]
), fill = TRUE)

heat_dt[, subtype := factor(subtype, levels = subtype_order)]
heat_dt[, row_label := factor(row_label, levels = row_order)]

p_main <- ggplot(heat_dt, aes(x = subtype, y = row_label, fill = value)) +
  geom_tile(color = "white", linewidth = 0.35) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    name = "Subtype\nscore"
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "Subtype-specific GSEA pathways",
    subtitle = paste0(
      "Unified GSEA pool: GO_BP / GO_CC / GO_MF / KEGG / H / C2 / C5 / C7; top ",
      n_show,
      " pathways shown per subtype"
    )
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(face = "bold")
  )

save_plot_pdf(
  p_main,
  "04_1_GSEA_subtype_specific_pathways_unified_heatmap",
  w = 7.8,
  h = max(5, 0.20 * length(unique(heat_dt$row_label)))
)

# ---------- export and plot collection subsets from the SAME unified selected table ----------
for (collection_name in names(gsea_sets)) {
  sub_tab <- out_master[collection == collection_name]
  if (nrow(sub_tab) == 0) next
  
  save_dat_overwrite(
    as.data.frame(sub_tab),
    paste0("04_1_GSEA_subtype_specific_pathways_unified_", collection_name, "_subset"),
    ext = "csv"
  )
  
  sub_sel <- gsea_sel_top[collection == collection_name]
  if (nrow(sub_sel) == 0) next
  
  sub_sel[, row_label_sub := paste0("[", assigned_subtype, "] ", Description)]
  
  row_order_sub_dt <- unique(
    sub_sel[, .(row_label_sub, assigned_subtype, subtype_rank)][
      order(factor(assigned_subtype, levels = subtype_order), subtype_rank)
    ]
  )
  row_order_sub <- rev(row_order_sub_dt$row_label_sub)
  
  heat_sub <- rbindlist(list(
    sub_sel[, .(row_label_sub, subtype = "ID", value = score_ID)],
    sub_sel[, .(row_label_sub, subtype = "IE", value = score_IE)],
    sub_sel[, .(row_label_sub, subtype = "IN", value = score_IN)]
  ), fill = TRUE)
  
  heat_sub[, subtype := factor(subtype, levels = subtype_order)]
  heat_sub[, row_label_sub := factor(row_label_sub, levels = row_order_sub)]
  
  p_sub <- ggplot(heat_sub, aes(x = subtype, y = row_label_sub, fill = value)) +
    geom_tile(color = "white", linewidth = 0.35) +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      name = "Subtype\nscore"
    ) +
    labs(
      x = NULL,
      y = NULL,
      title = paste0(collection_name, " pathways from unified selected GSEA set"),
      subtitle = "Subset extracted from the unified top pathways"
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(size = 12, face = "bold"),
      axis.text.y = element_text(size = 8),
      plot.title = element_text(face = "bold")
    )
  
  save_plot_pdf(
    p_sub,
    paste0("04_1_GSEA_subtype_specific_pathways_unified_", collection_name, "_subset_heatmap"),
    w = 6.6,
    h = max(4.5, 0.22 * length(unique(heat_sub$row_label_sub)))
  )
}

message("Section 12 finished successfully.")

# =========================================================
# 12 build one complete integrated GSEA table
# =========================================================

suppressPackageStartupMessages({
  library(data.table)
})

# ---------- parameters ----------
subtype_order <- c("ID", "IE", "IN")
contrast_order <- c("IE_vs_IN", "ID_vs_IN", "ID_vs_IE")
gsea_padj_cutoff <- 0.20

# ---------- GSEA result files ----------
gsea_sets <- list(
  GO_BP = c("04_1_GSEA_GO_BP_IE_vs_IN.rds", "04_1_GSEA_GO_BP_ID_vs_IN.rds", "04_1_GSEA_GO_BP_ID_vs_IE.rds"),
  GO_CC = c("04_1_GSEA_GO_CC_IE_vs_IN.rds", "04_1_GSEA_GO_CC_ID_vs_IN.rds", "04_1_GSEA_GO_CC_ID_vs_IE.rds"),
  GO_MF = c("04_1_GSEA_GO_MF_IE_vs_IN.rds", "04_1_GSEA_GO_MF_ID_vs_IN.rds", "04_1_GSEA_GO_MF_ID_vs_IE.rds"),
  KEGG  = c("04_1_GSEA_KEGG_IE_vs_IN.rds",  "04_1_GSEA_KEGG_ID_vs_IN.rds",  "04_1_GSEA_KEGG_ID_vs_IE.rds"),
  H     = c("04_1_GSEA_H_IE_vs_IN.rds",     "04_1_GSEA_H_ID_vs_IN.rds",     "04_1_GSEA_H_ID_vs_IE.rds"),
  C2    = c("04_1_GSEA_C2_IE_vs_IN.rds",    "04_1_GSEA_C2_ID_vs_IN.rds",    "04_1_GSEA_C2_ID_vs_IE.rds"),
  C5    = c("04_1_GSEA_C5_IE_vs_IN.rds",    "04_1_GSEA_C5_ID_vs_IN.rds",    "04_1_GSEA_C5_ID_vs_IE.rds"),
  C7    = c("04_1_GSEA_C7_IE_vs_IN.rds",    "04_1_GSEA_C7_ID_vs_IN.rds",    "04_1_GSEA_C7_ID_vs_IE.rds")
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
    pathway_id = as.character(ID),
    Description = as.character(Description),
    NES = as.numeric(NES),
    p.adjust = as.numeric(p.adjust),
    core_enrichment = if ("core_enrichment" %in% names(dt)) as.character(core_enrichment) else NA_character_
  )]
  dt[, collection := collection_name]
  dt[, contrast := contrast_name]
  dt
}

# ---------- load all ----------
gsea_list <- list()
for (set_nm in names(gsea_sets)) {
  file_vec <- gsea_sets[[set_nm]]
  for (i in seq_along(file_vec)) {
    fp <- file.path(dat_dir, file_vec[i])
    cn <- contrast_order[i]
    gsea_list[[paste(set_nm, cn, sep = "__")]] <- read_gsea_result(fp, set_nm, cn)
  }
}
gsea_long <- rbindlist(gsea_list, fill = TRUE)
if (nrow(gsea_long) == 0) stop("No GSEA files found for section 12.")

save_dat_overwrite(as.data.frame(gsea_long), "04_1_GSEA_complete_long_table", ext = "csv")
save_dat_overwrite(as.data.frame(gsea_long), "04_1_GSEA_complete_long_table", ext = "rds")

# ---------- wide ----------
gsea_wide <- dcast(
  gsea_long,
  collection + pathway_id + Description ~ contrast,
  value.var = c("NES", "p.adjust", "core_enrichment")
)

# subtype-level scores
gsea_wide[, score_ID := rowMeans(cbind(NES_ID_vs_IN, NES_ID_vs_IE), na.rm = TRUE)]
gsea_wide[, score_IE := rowMeans(cbind(NES_IE_vs_IN, -NES_ID_vs_IE), na.rm = TRUE)]
gsea_wide[, score_IN := rowMeans(cbind(-NES_IE_vs_IN, -NES_ID_vs_IN), na.rm = TRUE)]

gsea_wide[, rank_ID := pmin(NES_ID_vs_IN, NES_ID_vs_IE, na.rm = TRUE)]
gsea_wide[, rank_IE := pmin(NES_IE_vs_IN, -NES_ID_vs_IE, na.rm = TRUE)]
gsea_wide[, rank_IN := pmin(-NES_IE_vs_IN, -NES_ID_vs_IN, na.rm = TRUE)]

# subtype-specific label
gsea_wide[, sig_ID :=
            !is.na(p.adjust_ID_vs_IN) & !is.na(p.adjust_ID_vs_IE) &
            p.adjust_ID_vs_IN < gsea_padj_cutoff &
            p.adjust_ID_vs_IE < gsea_padj_cutoff &
            NES_ID_vs_IN > 0 & NES_ID_vs_IE > 0
]
gsea_wide[, sig_IE :=
            !is.na(p.adjust_IE_vs_IN) & !is.na(p.adjust_ID_vs_IE) &
            p.adjust_IE_vs_IN < gsea_padj_cutoff &
            p.adjust_ID_vs_IE < gsea_padj_cutoff &
            NES_IE_vs_IN > 0 & NES_ID_vs_IE < 0
]
gsea_wide[, sig_IN :=
            !is.na(p.adjust_IE_vs_IN) & !is.na(p.adjust_ID_vs_IN) &
            p.adjust_IE_vs_IN < gsea_padj_cutoff &
            p.adjust_ID_vs_IN < gsea_padj_cutoff &
            NES_IE_vs_IN < 0 & NES_ID_vs_IN < 0
]

gsea_wide[, assigned_subtype := NA_character_]
gsea_wide[sig_ID == TRUE, assigned_subtype := "ID"]
gsea_wide[sig_IE == TRUE, assigned_subtype := "IE"]
gsea_wide[sig_IN == TRUE, assigned_subtype := "IN"]

overlap_idx <- gsea_wide[rowSums(cbind(sig_ID, sig_IE, sig_IN), na.rm = TRUE) > 1, which = TRUE]
if (length(overlap_idx) > 0) {
  tmp <- gsea_wide[overlap_idx, .(score_ID, score_IE, score_IN)]
  best <- apply(as.matrix(tmp), 1, function(x) c("ID", "IE", "IN")[which.max(x)])
  gsea_wide[overlap_idx, assigned_subtype := best]
}

# best subtype by score, independent of significance rule
gsea_wide[, best_subtype_by_score := c("ID", "IE", "IN")[max.col(.SD, ties.method = "first")],
          .SDcols = c("score_ID", "score_IE", "score_IN")]

# complete integrated table
gsea_complete <- copy(gsea_wide)[, .(
  collection,
  pathway_id,
  Description,
  assigned_subtype,
  best_subtype_by_score,
  score_ID,
  score_IE,
  score_IN,
  rank_ID,
  rank_IE,
  rank_IN,
  NES_IE_vs_IN,
  p.adjust_IE_vs_IN,
  core_enrichment_IE_vs_IN,
  NES_ID_vs_IN,
  p.adjust_ID_vs_IN,
  core_enrichment_ID_vs_IN,
  NES_ID_vs_IE,
  p.adjust_ID_vs_IE,
  core_enrichment_ID_vs_IE
)]

save_dat_overwrite(as.data.frame(gsea_complete), "04_1_GSEA_complete_integrated_table", ext = "csv")
save_dat_overwrite(as.data.frame(gsea_complete), "04_1_GSEA_complete_integrated_table", ext = "rds")

message("Section 12 finished successfully.")


# =========================================================
# 13 custom selected pathway heatmap from
# 04_1_GSEA_complete_integrated_table
# 3-column layout: ID / IE / IN
# score definition:
#   score_ID = mean(NES_ID_vs_IN, NES_ID_vs_IE)
#   score_IE = mean(NES_IE_vs_IN, -NES_ID_vs_IE)
#   score_IN = mean(-NES_IE_vs_IN, -NES_ID_vs_IN)
# =========================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ---------- explicit path ----------
proj <- tryCatch(normalizePath(getwd(), winslash = "/"), error = function(e) NA)
if (is.na(proj) || !file.exists(file.path(proj, "results", "04_1", "data"))) {
  proj <- find_project_root()
}
dat_dir_04_1 <- file.path(proj, "results", "04_1", "data")
fig_dir_04_1 <- file.path(proj, "results", "04_1", "figures")

dir.create(fig_dir_04_1, recursive = TRUE, showWarnings = FALSE)

save_dat_overwrite_local <- function(x, filename, ext = c("csv", "rds")) {
  ext <- match.arg(ext)
  out <- file.path(dat_dir_04_1, paste0(filename, ".", ext))
  
  if (file.exists(out)) {
    ok_rm <- tryCatch(file.remove(out), warning = function(e) FALSE, error = function(e) FALSE)
    if (!isTRUE(ok_rm)) {
      stop("Cannot overwrite existing file: ", out)
    }
  }
  
  if (ext == "csv") write.csv(x, out, row.names = FALSE)
  if (ext == "rds") saveRDS(x, out)
  
  message("Saved data: ", out)
  invisible(out)
}

save_plot_pdf_local <- function(plot_obj, filename, w = 7, h = 6) {
  out <- file.path(fig_dir_04_1, paste0(filename, ".pdf"))
  
  if (file.exists(out)) {
    ok_rm <- tryCatch(file.remove(out), warning = function(e) FALSE, error = function(e) FALSE)
    if (!isTRUE(ok_rm)) {
      stop("Cannot overwrite existing pdf: ", out)
    }
  }
  
  ggplot2::ggsave(out, plot = plot_obj, width = w, height = h, device = cairo_pdf)
  message("Saved figure: ", out)
  invisible(out)
}

# ---------- read integrated table ----------
gsea_complete_file_rds <- file.path(dat_dir_04_1, "04_1_GSEA_complete_integrated_table.rds")
gsea_complete_file_csv <- file.path(dat_dir_04_1, "04_1_GSEA_complete_integrated_table.csv")

if (file.exists(gsea_complete_file_rds)) {
  gsea_complete <- as.data.table(readRDS(gsea_complete_file_rds))
} else if (file.exists(gsea_complete_file_csv)) {
  gsea_complete <- fread(gsea_complete_file_csv)
} else {
  stop("Missing integrated GSEA table in: ", dat_dir_04_1)
}

# ---------- required columns ----------
need_cols <- c(
  "collection", "pathway_id", "Description",
  "NES_IE_vs_IN", "p.adjust_IE_vs_IN", "core_enrichment_IE_vs_IN",
  "NES_ID_vs_IN", "p.adjust_ID_vs_IN", "core_enrichment_ID_vs_IN",
  "NES_ID_vs_IE", "p.adjust_ID_vs_IE", "core_enrichment_ID_vs_IE"
)
miss <- setdiff(need_cols, names(gsea_complete))
if (length(miss) > 0) {
  stop("Missing columns in 04_1_GSEA_complete_integrated_table: ", paste(miss, collapse = ", "))
}

# ---------- user input ----------
selected_pathways <- data.table(
  pathway_id = c(
    "GOCC_FIBRILLAR_COLLAGEN_TRIMER",
    "WP_TGFB_SMAD_SIGNALING",
    "GOMF_INSULIN_LIKE_GROWTH_FACTOR_BINDING",
    "NIKOLSKY_BREAST_CANCER_15Q26_AMPLICON",
    "GOBP_POSITIVE_REGULATION_OF_EXTRACELLULAR_MATRIX_ASSEMBLY",
    
    "hsa04080",
    "REACTOME_GPCR_LIGAND_BINDING",
    "GOMF_NEUROPEPTIDE_RECEPTOR_ACTIVITY",
    "GO:0030594",
    "GOMF_FIBROBLAST_GROWTH_FACTOR_RECEPTOR_BINDING",
    
    "MOSERLE_IFNA_RESPONSE",
    "GO:0023026",
    "GO:0042611",
    "GO:0002501",
    "GO:0002495"
  ),
  display_name = c(
    "Fibrillar collagen trimer",
    "TGFB-SMAD signaling",
    "IGF binding",
    "Breast cancer amplicon",
    "Extracellular matrix assembly",
    
    "Neuroactive ligand-receptor interaction",
    "GPCR ligand binding",
    "Neuropeptide receptor activity",
    "Neurotransmitter receptor activity",
    "FGFR binding",
    
    "IFNA response",
    "MHCII protein complex binding",
    "MHC protein complex",
    "Peptide antigen assembly with MHC protein complex",
    "Antigen processing and presentation"
  )
)

if (nrow(selected_pathways) == 0) {
  stop("Please fill selected_pathways first.")
}
selected_pathways[, input_order := seq_len(.N)]

# ---------- subset ----------
plot_dt <- merge(
  selected_pathways,
  gsea_complete,
  by = "pathway_id",
  all.x = TRUE
)

plot_dt <- plot_dt[!is.na(collection)]
if (nrow(plot_dt) == 0) {
  stop("No selected pathway_id matched in 04_1_GSEA_complete_integrated_table.")
}

dup_name <- plot_dt[, .N, by = display_name][N > 1, display_name]
plot_dt[, row_label := ifelse(
  display_name %in% dup_name,
  paste0(display_name, " [", collection, "]"),
  display_name
)]

# ---------- calculate mean-based subtype scores ----------
plot_dt[, score_ID := rowMeans(cbind(NES_ID_vs_IN, NES_ID_vs_IE), na.rm = TRUE)]
plot_dt[, score_IE := rowMeans(cbind(NES_IE_vs_IN, -NES_ID_vs_IE), na.rm = TRUE)]
plot_dt[, score_IN := rowMeans(cbind(-NES_IE_vs_IN, -NES_ID_vs_IN), na.rm = TRUE)]

plot_dt[, assigned_subtype := c("ID", "IE", "IN")[max.col(.SD, ties.method = "first")],
        .SDcols = c("score_ID", "score_IE", "score_IN")]

# ---------- export selected table ----------
out_selected <- unique(plot_dt[, .(
  pathway_id,
  display_name,
  collection,
  Description,
  row_label,
  input_order,
  assigned_subtype,
  score_ID,
  score_IE,
  score_IN,
  NES_IE_vs_IN,
  p.adjust_IE_vs_IN,
  core_enrichment_IE_vs_IN,
  NES_ID_vs_IN,
  p.adjust_ID_vs_IN,
  core_enrichment_ID_vs_IN,
  NES_ID_vs_IE,
  p.adjust_ID_vs_IE,
  core_enrichment_ID_vs_IE
)])

save_dat_overwrite_local(as.data.frame(out_selected), "04_1_GSEA_custom_selected_pathways_table", ext = "csv")
save_dat_overwrite_local(as.data.frame(out_selected), "04_1_GSEA_custom_selected_pathways_table", ext = "rds")

# ---------- build heatmap long table ----------
heat_dt <- rbindlist(list(
  plot_dt[, .(row_label, input_order, assigned_subtype, subtype = "ID", score = score_ID)],
  plot_dt[, .(row_label, input_order, assigned_subtype, subtype = "IE", score = score_IE)],
  plot_dt[, .(row_label, input_order, assigned_subtype, subtype = "IN", score = score_IN)]
), fill = TRUE)

heat_dt[, subtype := factor(subtype, levels = c("ID", "IE", "IN"))]

row_order_dt <- unique(
  heat_dt[, .(row_label, assigned_subtype, input_order)][
    order(factor(assigned_subtype, levels = c("ID", "IE", "IN")), input_order)
  ]
)
row_order <- rev(row_order_dt$row_label)
heat_dt[, row_label := factor(row_label, levels = row_order)]

# ---------- heatmap ----------
p_custom <- ggplot(heat_dt, aes(x = subtype, y = row_label, fill = score)) +
  geom_tile(color = "white", linewidth = 0.4) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    na.value = "grey92",
    name = "Score"
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "Selected GSEA pathways",
    subtitle = "Mean score derived from pairwise NES"
  ) +
  theme_bw(base_size = 9) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 9, face = "bold"),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(face = "bold", size = 10),
    plot.subtitle = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7)
  )

save_plot_pdf_local(
  p_custom,
  "04_1_GSEA_custom_selected_pathways_heatmap_3col",
  w = 5.6,
  h = max(4.2, 0.26 * length(unique(heat_dt$row_label)))
)

message("Section 13 finished successfully.")

message("04_1 duplicateCorrelation DEG + GO/KEGG/GSEA + GMT-GSEA pipeline finished successfully.")


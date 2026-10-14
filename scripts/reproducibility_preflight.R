# scripts/reproducibility_preflight.R
# Lightweight reproducibility preflight for the TNBC CosMx 6k project.
# Run from the repository root (or any subdirectory inside the repository).
#
# By default this checks renv activation, package availability, frozen input
# filenames + SHA-256 identities, finalized metadata structure, XY/UM cell-ID
# coverage, signature structure, GMT files, and Cairo PDF capability.
# Set check_large_expression <- TRUE to additionally read the ~400 MB normalized
# expression RDS and validate its dimensions/cell overlap.

rm(list = ls())
gc()
options(stringsAsFactors = FALSE)

check_large_expression <- FALSE
check_core_hashes <- TRUE
# Opt in only after reconstructing metadata with 02_3. This does not establish
# equality with the audited metadata; all structural checks remain required.
rebuilt_metadata <- isTRUE(getOption("tnbc.rebuilt_metadata", FALSE))

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
    "Open the R project or run this script from the repository root/subdirectory."
  )
}

proj <- find_project_root()
raw_dir <- file.path(proj, "raw")
gene_set_dir <- file.path(proj, "gene_sets")
out_dir <- file.path(proj, "results", "preflight")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

report <- character()
add_report <- function(...) {
  msg <- paste0(...)
  report <<- c(report, msg)
  message(msg)
}

stop_with_report <- function(...) {
  add_report("ERROR: ", ...)
  writeLines(report, file.path(out_dir, "reproducibility_preflight.txt"))
  stop(paste0(...), call. = FALSE)
}

add_report("Project root: ", proj)
add_report("Preflight time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
add_report("R version: ", R.version.string)
expected_r_version <- "4.5.1"
if (!identical(paste(R.version$major, R.version$minor, sep = "."), expected_r_version)) {
  add_report(
    "WARN: audited environment used R ", expected_r_version,
    "; current R is ", paste(R.version$major, R.version$minor, sep = ".")
  )
}

renv_root <- normalizePath(file.path(proj, "renv", "library"), winslash = "/", mustWork = FALSE)
lib_paths_norm <- normalizePath(.libPaths(), winslash = "/", mustWork = FALSE)
renv_active <- any(startsWith(lib_paths_norm, renv_root))
add_report("Project renv library active: ", renv_active)
if (!renv_active) {
  stop_with_report(
    "Project renv library is not active. Open the .Rproj and restore/activate renv before preflight."
  )
}

cairo_ok <- isTRUE(capabilities("cairo"))
add_report("Cairo graphics capability: ", cairo_ok)
if (!cairo_ok) {
  stop_with_report("R was built without Cairo support, but the workflow writes Cairo PDF figures.")
}
add_report("")

core_packages <- c(
  "data.table", "Matrix", "ggplot2", "Seurat", "scales",
  "GSVA", "BiocParallel", "SummarizedExperiment", "matrixStats", "readxl",
  "emmeans", "lme4", "lmerTest", "limma",
  "clusterProfiler", "org.Hs.eg.db", "enrichplot", "fgsea",
  "cowplot", "deldir", "digest"
)

check_packages <- function(pkgs, label) {
  add_report("== ", label, " ==")
  status <- vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)
  for (pkg in pkgs) {
    if (status[[pkg]]) {
      add_report("OK   ", pkg, " ", as.character(utils::packageVersion(pkg)))
    } else {
      add_report("MISS ", pkg)
    }
  }
  add_report("")
  invisible(status)
}

core_status <- check_packages(core_packages, "Canonical workflow packages")

core_missing <- names(core_status)[!core_status]
if (length(core_missing) > 0) {
  add_report("Canonical workflow packages missing: ", paste(core_missing, collapse = ", "))
  add_report("")
}

core_inputs <- c(
  "TNBC_6k_metadata_Bogang_selected.RDS",
  "TNBC_6k_norm_Bogang.RDS",
  "TNBC_6k_xy_Bogang.RDS",
  "TNBC_6k_um_Bogang.RDS",
  "selected_subtype_specific_genes_tumor.csv",
  "T_cell_signalling.xlsx"
)

add_report("== Core raw inputs ==")
core_paths <- file.path(raw_dir, core_inputs)
xlsx_index <- match("T_cell_signalling.xlsx", core_inputs)
if (!file.exists(core_paths[xlsx_index])) {
  core_paths[xlsx_index] <- file.path(gene_set_dir, "T_cell_signalling.xlsx")
}
for (i in seq_along(core_inputs)) {
  if (file.exists(core_paths[[i]])) {
    size_mb <- file.info(core_paths[[i]])$size / 1024^2
    add_report("OK   ", core_inputs[[i]], " (", format(round(size_mb, 1), nsmall = 1), " MB)")
  } else {
    add_report("MISS ", core_inputs[[i]])
  }
}
add_report("")

if (!all(file.exists(core_paths))) {
  stop_with_report(
    "Missing core raw inputs: ",
    paste(core_inputs[!file.exists(core_paths)], collapse = ", ")
  )
}

add_report("== Core raw SHA-256 identities ==")
hash_manifest <- file.path(proj, "environment", "raw_core_sha256.txt")

if (check_core_hashes) {
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop_with_report("Package 'digest' is required for SHA-256 input verification.")
  }
  if (!file.exists(hash_manifest)) {
    stop_with_report("Missing SHA-256 manifest: environment/raw_core_sha256.txt")
  }
  
  manifest_lines <- trimws(readLines(hash_manifest, warn = FALSE))
  manifest_lines <- manifest_lines[
    nzchar(manifest_lines) & !startsWith(manifest_lines, "#")
  ]
  manifest_parts <- strsplit(manifest_lines, "[[:space:]]+")
  expected_sha256 <- setNames(
    vapply(manifest_parts, function(x) x[[1]], character(1)),
    vapply(manifest_parts, function(x) paste(x[-1], collapse = " "), character(1))
  )
  
  missing_manifest <- setdiff(core_inputs, names(expected_sha256))
  if (length(missing_manifest) > 0) {
    stop_with_report(
      "SHA-256 manifest is missing core inputs: ",
      paste(missing_manifest, collapse = ", ")
    )
  }
  
  observed_sha256 <- vapply(
    core_paths,
    function(fp) digest::digest(fp, algo = "sha256", file = TRUE),
    character(1)
  )
  names(observed_sha256) <- core_inputs
  
  for (nm in core_inputs) {
    ok <- identical(unname(observed_sha256[[nm]]), unname(expected_sha256[[nm]]))
    add_report(if (ok) "OK   " else "MISMATCH ", nm)
  }
  
  bad_hash <- core_inputs[observed_sha256[core_inputs] != expected_sha256[core_inputs]]
  if (rebuilt_metadata) {
    add_report("WARN: rebuilt metadata mode: metadata file identity is not required; biological equality is NOT established by structural checks.")
    add_report("Observed metadata SHA-256: ", observed_sha256[["TNBC_6k_metadata_Bogang_selected.RDS"]])
    bad_hash <- setdiff(bad_hash, "TNBC_6k_metadata_Bogang_selected.RDS")
  }
  if (length(bad_hash) > 0) {
    stop_with_report(
      "Core raw input SHA-256 mismatch: ",
      paste(bad_hash, collapse = ", "),
      ". Confirm that the exact audited raw files were supplied."
    )
  }
} else {
  add_report("SKIP core raw SHA-256 verification (check_core_hashes = FALSE)")
}
add_report("")

gmt_files <- c(
  "h.all.v2025.1.Hs.symbols.gmt",
  "c2.all.v2025.1.Hs.symbols.gmt",
  "c5.all.v2025.1.Hs.symbols.gmt",
  "c7.all.v2025.1.Hs.symbols.gmt"
)

add_report("== Versioned GMT resources ==")
gmt_paths <- file.path(gene_set_dir, gmt_files)
gmt_ok <- file.exists(gmt_paths) & file.info(gmt_paths)$size > 0
for (i in seq_along(gmt_files)) {
  add_report(if (gmt_ok[[i]]) "OK   " else "MISS/EMPTY ", gmt_files[[i]])
}
add_report("")

if (!all(gmt_ok)) {
  stop_with_report("One or more required GMT files are missing or empty.")
}

# Finalized metadata snapshot checks.
meta <- readRDS(file.path(raw_dir, "TNBC_6k_metadata_Bogang_selected.RDS"))
meta_df <- as.data.frame(meta)

required_meta_cols <- c(
  "cell_id", "Run_Tissue_name", "FOV", "BC6k_refined",
  "celltype_new", "disttonearesttumor"
)
missing_meta_cols <- setdiff(required_meta_cols, names(meta_df))
if (length(missing_meta_cols) > 0) {
  stop_with_report(
    "Finalized metadata is missing required columns: ",
    paste(missing_meta_cols, collapse = ", ")
  )
}

n_cells <- nrow(meta_df)
n_fov <- length(unique(meta_df$FOV))
n_slides <- length(unique(meta_df$Run_Tissue_name))
n_dup_cell <- sum(duplicated(meta_df$cell_id))

add_report("== Finalized metadata ==")
add_report("Cells: ", n_cells, " (expected frozen snapshot: 275204)")
add_report("FOVs: ", n_fov, " (expected frozen snapshot: 172)")
add_report("Slides: ", n_slides, " (expected frozen snapshot: 5)")
add_report("Duplicated cell_id: ", n_dup_cell)

if (n_cells != 275204L || n_fov != 172L || n_slides != 5L || n_dup_cell != 0L) {
  stop_with_report(
    "Finalized metadata does not match the audited 275204-cell / 172-FOV / 5-slide snapshot."
  )
}

broad_expected <- c(
  "B.cell", "Dendritic.cell", "Endothelial", "Fibroblast",
  "Macrophage", "Mast.cell", "Monocyte", "Neutrophil", "NK.cell",
  "Plasma", "T.cell.CD4", "T.cell.CD8", "T.cell.regulatory", "tumor"
)
missing_broad <- setdiff(broad_expected, unique(as.character(meta_df$BC6k_refined)))
if (length(missing_broad) > 0) {
  stop_with_report("Missing expected broad cell types: ", paste(missing_broad, collapse = ", "))
}

cd8_expected <- c(
  "cd8_cytotoxic", "cd8_exhausted", "cd8_naive", "cd8_tcm", "cd8_tem"
)
missing_cd8 <- setdiff(cd8_expected, unique(as.character(meta_df$celltype_new)))
if (length(missing_cd8) > 0) {
  stop_with_report("Missing expected fine CD8 subtypes: ", paste(missing_cd8, collapse = ", "))
}

dist <- suppressWarnings(as.numeric(meta_df$disttonearesttumor))
if (all(is.na(dist))) stop_with_report("disttonearesttumor cannot be interpreted as numeric.")
if (any(dist < 0, na.rm = TRUE)) stop_with_report("disttonearesttumor contains negative values.")
add_report("disttonearesttumor median (um): ", format(round(stats::median(dist, na.rm = TRUE), 3), nsmall = 3))
add_report("")

# XY and UMAP objects are small enough to validate by default.
xy <- readRDS(file.path(raw_dir, "TNBC_6k_xy_Bogang.RDS"))
um <- readRDS(file.path(raw_dir, "TNBC_6k_um_Bogang.RDS"))

xy_ids <- rownames(xy)
um_ids <- rownames(um)
if (is.null(xy_ids) || is.null(um_ids)) {
  stop_with_report("XY/UM objects must store cell IDs as row names.")
}

missing_xy <- sum(!meta_df$cell_id %in% xy_ids)
missing_um <- sum(!meta_df$cell_id %in% um_ids)
add_report("== XY / UMAP coverage ==")
add_report("Metadata cells missing from XY: ", missing_xy)
add_report("Metadata cells missing from UMAP: ", missing_um)
if (missing_xy > 0L || missing_um > 0L) {
  stop_with_report("Finalized metadata is not fully covered by XY/UM cell IDs.")
}
add_report("")

# Signature file structure.
sig <- utils::read.csv(
  file.path(raw_dir, "selected_subtype_specific_genes_tumor.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
if (!all(c("Gene", "Subtype") %in% names(sig))) {
  stop_with_report("Signature CSV must contain Gene and Subtype columns.")
}
sig_counts <- table(sig$Subtype)
add_report("== Tumor immunotype signature ==")
for (nm in names(sig_counts)) add_report(nm, ": ", unname(sig_counts[[nm]]), " genes")
expected_sig_counts <- c(tumor_desert = 40L, tumor_excluded = 35L, tumor_inflamed = 40L)
if (!all(names(expected_sig_counts) %in% names(sig_counts)) ||
    any(as.integer(sig_counts[names(expected_sig_counts)]) != expected_sig_counts)) {
  stop_with_report("Signature CSV does not match the audited 40/35/40 subtype definition.")
}
add_report("")

if (check_large_expression) {
  add_report("== Large normalized-expression check ==")
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop_with_report("Matrix package is required for the large-expression check.")
  }
  norm <- readRDS(file.path(raw_dir, "TNBC_6k_norm_Bogang.RDS"))
  d <- dim(norm)
  add_report("Expression dimensions: ", paste(d, collapse = " x "), " (expected 294923 x 6175)")
  if (length(d) != 2L || !identical(as.integer(d), c(294923L, 6175L))) {
    stop_with_report("Normalized-expression dimensions do not match the audited snapshot.")
  }
  expr_ids <- rownames(norm)
  if (is.null(expr_ids)) stop_with_report("Normalized expression must store cell IDs as row names.")
  miss_expr <- sum(!meta_df$cell_id %in% expr_ids)
  add_report("Metadata cells missing from expression: ", miss_expr)
  if (miss_expr > 0L) stop_with_report("Finalized metadata is not fully covered by expression cell IDs.")
  add_report("")
}

writeLines(report, file.path(out_dir, "reproducibility_preflight.txt"))

if (length(core_missing) > 0) {
  stop(
    "Input structure passed, but canonical R packages are missing: ",
    paste(core_missing, collapse = ", "),
    ". See results/preflight/reproducibility_preflight.txt",
    call. = FALSE
  )
}

message("")
message("Preflight PASSED.")
message("Report: ", file.path(out_dir, "reproducibility_preflight.txt"))
if (!check_large_expression) {
  message("For the optional full expression-object check, set check_large_expression <- TRUE and rerun.")
}



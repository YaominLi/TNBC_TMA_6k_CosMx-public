# Reconstruct finalized metadata from the original metadata and Yi update.
# Run once before preflight/03_2; already finalized inputs should skip this step.
rm(list = ls())
gc()
graphics.off()

options(stringsAsFactors = FALSE)

# =========================================================
# 00 paths
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

proj <- find_project_root()

raw_dir <- file.path(proj, "raw")
backup_dir <- file.path(raw_dir, "backup_before_02_2")
out_dir <- file.path(proj, "results", "02_3", "data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# The immutable backup is preferred. For a fresh clone, the target may hold
# the original metadata; it is backed up before any replacement.
dir.create(backup_dir, recursive = TRUE, showWarnings = FALSE)
# clean inputs only
old_meta_backup_file <- file.path(backup_dir, "TNBC_6k_metadata_Bogang_selected_before_02_2.RDS")
new_meta_file <- file.path(raw_dir, "TNBC_6k_metadata_172FOVs_Bogang.RDS")

# overwrite target
target_meta_file <- file.path(raw_dir, "TNBC_6k_metadata_Bogang_selected.RDS")

save_dat_overwrite <- function(x, filename, ext = c("csv", "rds", "txt")) {
  ext <- match.arg(ext)
  out <- file.path(out_dir, paste0(filename, ".", ext))
  if (file.exists(out)) file.remove(out)
  
  if (ext == "csv") write.csv(x, out, row.names = FALSE)
  if (ext == "rds") saveRDS(x, out)
  if (ext == "txt") writeLines(x, out)
  
  message("Saved: ", out)
  invisible(out)
}

# =========================================================
# 01 packages
# =========================================================
suppressPackageStartupMessages({
  library(data.table)
})

# =========================================================
# 02 checks and load
# =========================================================
if (!file.exists(old_meta_backup_file)) {
  if (!file.exists(target_meta_file)) stop("Supply the original metadata or its pre-update backup.")
  candidate <- as.data.frame(readRDS(target_meta_file))
  if (any(c("celltype_new", "BC6k_refined_old", "disttonearesttumor_old") %in% names(candidate))) {
    stop("Metadata is already processed. Skip 02_3, or supply the original pre-update backup. No file changed.")
  }
  if (!file.copy(target_meta_file, old_meta_backup_file, overwrite = FALSE)) {
    stop("Cannot preserve the original metadata backup.")
  }
  rm(candidate)
}
if (!file.exists(new_meta_file)) {
  stop("Missing new metadata file: ", new_meta_file)
}

old_dt <- as.data.table(readRDS(old_meta_backup_file))
new_dt <- as.data.table(readRDS(new_meta_file))

need_old <- c(
  "cell_id", "FOV", "Run_Tissue_name",
  "BC6k_refined", "disttonearesttumor", "immunotype_BW"
)

need_new <- c(
  "cell_id", "FOV", "Run_Tissue_name",
  "celltype", "disttonearestcancer"
)

miss_old <- setdiff(need_old, names(old_dt))
miss_new <- setdiff(need_new, names(new_dt))

if (length(miss_old) > 0) {
  stop("Missing required columns in backup old metadata: ", paste(miss_old, collapse = ", "))
}
if (length(miss_new) > 0) {
  stop("Missing required columns in new metadata: ", paste(miss_new, collapse = ", "))
}

old_dt[, cell_id := as.character(cell_id)]
new_dt[, cell_id := as.character(cell_id)]
old_dt[, FOV := as.character(FOV)]
new_dt[, FOV := as.character(FOV)]
old_dt[, Run_Tissue_name := as.character(Run_Tissue_name)]
new_dt[, Run_Tissue_name := as.character(Run_Tissue_name)]

# Reject ambiguous joins and processed backups before deriving new fields.
if (any(c("celltype_new", "BC6k_refined_old", "disttonearesttumor_old") %in% names(old_dt))) {
  stop("The backup is already processed; use the untouched original metadata.")
}
for (x in list(old_dt, new_dt)) {
  if (anyNA(x$cell_id) || any(!nzchar(x$cell_id)) || anyDuplicated(x$cell_id)) {
    stop("Input cell_id values must be nonmissing, nonempty and unique.")
  }
  if (anyNA(x$FOV) || anyNA(x$Run_Tissue_name)) stop("Missing FOV/slide identifiers.")
}
new_distance <- suppressWarnings(as.numeric(as.character(new_dt$disttonearestcancer)))
if (any(!is.na(new_dt$disttonearestcancer) & is.na(new_distance)) ||
    any(!is.finite(new_distance) & !is.na(new_distance)) ||
    any(new_distance < 0, na.rm = TRUE) || all(is.na(new_distance))) {
  stop("Invalid updated distance values; expected nonnegative millimetres (NA permitted).")
}
new_dt[, disttonearestcancer := new_distance]

# =========================================================
# 03 shared cells only
# =========================================================
shared_ids <- intersect(unique(old_dt$cell_id), unique(new_dt$cell_id))
if (length(shared_ids) == 0) {
  stop("No shared cell_id found between backup old and new metadata.")
}

old_sub <- old_dt[cell_id %in% shared_ids]
new_sub <- unique(
  new_dt[cell_id %in% shared_ids, .(
    cell_id,
    FOV_new = FOV,
    Run_Tissue_name_new = Run_Tissue_name,
    celltype_new = as.character(celltype),
    disttonearestcancer = as.numeric(disttonearestcancer)
  )],
  by = "cell_id"
)

dt <- merge(old_sub, new_sub, by = "cell_id", all = FALSE)

stopifnot("celltype_new" %in% names(dt))
stopifnot("disttonearestcancer" %in% names(dt))

# =========================================================
# 04 consistency check of shared structural fields
# =========================================================
field_check <- data.table(
  field = c("FOV", "Run_Tissue_name"),
  pct_equal = c(
    mean(as.character(dt$FOV) == as.character(dt$FOV_new)),
    mean(as.character(dt$Run_Tissue_name) == as.character(dt$Run_Tissue_name_new))
  )
)
save_dat_overwrite(as.data.frame(field_check), "02_3_shared_field_check", ext = "csv")
if (any(!is.finite(field_check$pct_equal)) || any(field_check$pct_equal != 1)) {
  stop("Original and updated FOV/slide identifiers disagree. Target metadata unchanged.")
}

# =========================================================
# 05 map fine celltype back to old-pipeline-compatible BC6k_refined
# keep fine annotation in celltype_new
# =========================================================
dt[, BC6k_refined_old := BC6k_refined]

map_to_old_like <- function(x) {
  x <- as.character(x)
  
  out <- x
  
  # tumor
  out[out %in% c("cancer")] <- "tumor"
  
  # B / plasma
  out[out %in% c("bcell")] <- "B.cell"
  out[out %in% c("plasma")] <- "Plasma"
  
  # dendritic / myeloid
  out[out %in% c("dendritic")] <- "Dendritic.cell"
  out[out %in% c("macrophage")] <- "Macrophage"
  out[out %in% c("monocyte")] <- "Monocyte"
  out[out %in% c("neutrophil")] <- "Neutrophil"
  out[out %in% c("mast")] <- "Mast.cell"
  
  # NK
  out[out %in% c("nk")] <- "NK.cell"
  
  # CD4
  out[out %in% c("cd4_naive", "cd4_tcm", "cd4_tem", "cd4_th1", "cd4_th17", "cd4_th2")] <- "T.cell.CD4"
  out[out %in% c("cd4_treg")] <- "T.cell.regulatory"
  
  # CD8
  out[out %in% c("cd8_cytotoxic", "cd8_exhausted", "cd8_naive", "cd8_tcm", "cd8_tem")] <- "T.cell.CD8"
  
  # stromal
  out[out %in% c("endothelial")] <- "Endothelial"
  out[out %in% c("fibroblast", "stroma0", "stroma1")] <- "Fibroblast"
  
  out
}

dt[, BC6k_refined := map_to_old_like(celltype_new)]
broad_allowed <- c("tumor", "B.cell", "Plasma", "Dendritic.cell", "Macrophage",
  "Monocyte", "Neutrophil", "Mast.cell", "NK.cell", "T.cell.CD4",
  "T.cell.regulatory", "T.cell.CD8", "Endothelial", "Fibroblast")
if (anyNA(dt$BC6k_refined) || any(!dt$BC6k_refined %in% broad_allowed)) {
  stop("Unrecognized updated cell labels; review mapping before proceeding.")
}
if (nrow(dt) != 275204L || uniqueN(dt$FOV) != 172L || uniqueN(dt$Run_Tissue_name) != 5L) {
  stop("Rebuilt metadata does not match the expected 275204 cells / 172 FOVs / 5 slides. Target unchanged.")
}

# =========================================================
# 06 overwrite distance variable using new definition
# convert to um for old pipeline compatibility
# =========================================================
dt[, disttonearesttumor_old := disttonearesttumor]
dt[, disttonearesttumor := as.numeric(disttonearestcancer) * 1000]

dt[, disttonearesttumor_3group := fifelse(
  is.na(disttonearesttumor), NA_character_,
  fifelse(
    disttonearesttumor <= 50, "<=50um",
    fifelse(disttonearesttumor <= 100, "50-100um", ">100um")
  )
)]

# =========================================================
# 07 drop merge helper columns not needed in final metadata
# keep celltype_new as fine annotation for future use
# =========================================================
drop_cols <- intersect(c("FOV_new", "Run_Tissue_name_new", "disttonearestcancer"), names(dt))
if (length(drop_cols) > 0) {
  dt[, (drop_cols) := NULL]
}

# =========================================================
# 08 reorder columns
# keep original structure first, append added columns after
# =========================================================
old_cols <- names(old_dt)
extra_cols <- setdiff(names(dt), old_cols)
setcolorder(dt, c(intersect(old_cols, names(dt)), extra_cols))

# =========================================================
# 09 summaries
# =========================================================
summary_dt <- data.table(
  n_cells_final = nrow(dt),
  n_FOV_final = uniqueN(dt$FOV),
  n_Run_Tissue_name_final = uniqueN(dt$Run_Tissue_name),
  n_old_cells_backup = nrow(old_dt),
  n_new_cells = nrow(new_dt),
  n_shared_cells = length(shared_ids)
)
save_dat_overwrite(as.data.frame(summary_dt), "02_3_final_summary", ext = "csv")

anno_change_old_vs_fine <- dt[, .N, by = .(BC6k_refined_old, celltype_new)][order(-N)]
save_dat_overwrite(as.data.frame(anno_change_old_vs_fine), "02_3_old_vs_new_fine_annotation", ext = "csv")

anno_change_old_vs_pipeline <- dt[, .N, by = .(BC6k_refined_old, BC6k_refined)][order(-N)]
save_dat_overwrite(as.data.frame(anno_change_old_vs_pipeline), "02_3_old_vs_pipeline_annotation", ext = "csv")

current_bc_levels <- data.table(BC6k_refined = sort(unique(dt$BC6k_refined)))
save_dat_overwrite(as.data.frame(current_bc_levels), "02_3_current_BC6k_refined_levels", ext = "csv")

dist_group_summary <- dt[, .N, by = disttonearesttumor_3group][order(disttonearesttumor_3group)]
save_dat_overwrite(as.data.frame(dist_group_summary), "02_3_distance_group_summary", ext = "csv")

dist_summary <- data.table(
  n = sum(!is.na(dt$disttonearesttumor)),
  min = min(dt$disttonearesttumor, na.rm = TRUE),
  q01 = quantile(dt$disttonearesttumor, 0.01, na.rm = TRUE),
  q05 = quantile(dt$disttonearesttumor, 0.05, na.rm = TRUE),
  q25 = quantile(dt$disttonearesttumor, 0.25, na.rm = TRUE),
  median = median(dt$disttonearesttumor, na.rm = TRUE),
  q75 = quantile(dt$disttonearesttumor, 0.75, na.rm = TRUE),
  q95 = quantile(dt$disttonearesttumor, 0.95, na.rm = TRUE),
  q99 = quantile(dt$disttonearesttumor, 0.99, na.rm = TRUE),
  max = max(dt$disttonearesttumor, na.rm = TRUE)
)
save_dat_overwrite(as.data.frame(dist_summary), "02_3_distance_summary", ext = "csv")

# =========================================================
# 10 backup current overwritten file if exists, then overwrite target
# =========================================================
if (file.exists(target_meta_file)) {
  current_overwritten_backup <- tempfile("metadata_before_prepare_", tmpdir = backup_dir, fileext = ".RDS")
  if (!file.copy(target_meta_file, current_overwritten_backup, overwrite = FALSE)) {
    stop("Cannot back up existing target metadata.")
  }
}

staged_meta <- tempfile("metadata_ready_", tmpdir = raw_dir, fileext = ".RDS")
saveRDS(dt, staged_meta)
stopifnot(identical(as.data.frame(readRDS(staged_meta)), as.data.frame(dt)))
if (!file.copy(staged_meta, target_meta_file, overwrite = TRUE)) stop("Cannot install finalized metadata.")
unlink(staged_meta)
message("Overwritten metadata file: ", target_meta_file)

save_dat_overwrite(dt, "02_3_overwritten_metadata_copy", ext = "rds")
save_dat_overwrite(as.data.frame(dt), "02_3_overwritten_metadata_copy", ext = "csv")

# =========================================================
# 11 report
# =========================================================
report_lines <- c(
  "02_3 rebuild compatible metadata from backup report",
  paste0("Input backup old metadata: ", old_meta_backup_file),
  paste0("Input new metadata: ", new_meta_file),
  paste0("Output overwritten metadata: ", target_meta_file),
  paste0("Original metadata backup: ", old_meta_backup_file),
  "",
  "[strategy]",
  "1. Start from backup old metadata only.",
  "2. Merge with Yi new metadata on shared cell_id only.",
  "3. Keep fine annotation in celltype_new.",
  "4. Collapse fine annotation back to old-pipeline-compatible BC6k_refined.",
  "5. Collapse all cd8_* classes into T.cell.CD8 for downstream 03+ compatibility.",
  "6. Overwrite disttonearesttumor using new disttonearestcancer and convert to um.",
  "7. Recalculate disttonearesttumor_3group using <=50, 50-100, >100.",
  "",
  "[summary]",
  paste(capture.output(print(summary_dt)), collapse = "\n"),
  "",
  "[shared field consistency]",
  paste(capture.output(print(field_check)), collapse = "\n"),
  "",
  "[current BC6k_refined levels]",
  paste(sort(unique(dt$BC6k_refined)), collapse = ", ")
)

save_dat_overwrite(report_lines, "02_3_report", ext = "txt")

message("02_3 finished successfully.")

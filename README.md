# TNBC TMA CosMx 6k analysis

R analysis code for the TNBC CosMx 6k project.

## Data layout

The `raw/` directory is intentionally excluded from Git. Place the five external core inputs under `raw/` at the repository root. The small T-cell gene-set workbook is bundled in `gene_sets/`.

The analysis uses these six core inputs:

| Input | Contents / role |
| --- | --- |
| `TNBC_6k_metadata_Bogang_selected.RDS` | Finalized cell annotations, slide/FOV identifiers, and tumor-distance metadata |
| `TNBC_6k_norm_Bogang.RDS` | Normalized expression used for signatures, pseudobulk, and expression analyses |
| `TNBC_6k_xy_Bogang.RDS` | Spatial cell coordinates |
| `TNBC_6k_um_Bogang.RDS` | UMAP embedding coordinates |
| `selected_subtype_specific_genes_tumor.csv` | ID/IE/IN tumor-signature gene lists |
| `gene_sets/T_cell_signalling.xlsx` | Bundled T-cell signaling gene sets; an existing `raw/T_cell_signalling.xlsx` takes precedence and is checked against the same hash |

The repository supplies code, environment definitions and the T-cell gene-set workbook; the RDS objects and tumor-signature CSV must still be supplied separately. Input identities are recorded in [environment/raw_core_sha256.txt](environment/raw_core_sha256.txt).

Enrichment uses MSigDB **2025.1 human gene-symbol GMT** resources. Obtain the four files listed in [gene_sets/README.md](gene_sets/README.md) before running preflight. GMT resources are not bundled in the public code distribution.

## Input provenance and annotation versions

The initial data package supplied by Yi comprised four preprocessed R objects: `TNBC_6k_metadata_Bogang_selected.RDS` (cell metadata), `TNBC_6k_norm_Bogang.RDS` (normalized expression), `TNBC_6k_um_Bogang.RDS` (UMAP), and `TNBC_6k_xy_Bogang.RDS` (spatial coordinates). Here, `raw/` means project inputs; these are not instrument-level raw reads or images. “Supplied by” identifies the data provider, not authorship of every upstream processing step.

Yi subsequently supplied `TNBC_6k_metadata_172FOVs_Bogang.RDS` (updated cell annotations and cancer-distance values) and `Immunotype_Wu_172FOVs_40genes.csv` (a supplied immunotype reference). The updated metadata is required only when rebuilding the finalized input with 02_3; the immunotype-reference CSV is not required by the current workflow. The immunotype reference CSV is distinct from `selected_subtype_specific_genes_tumor.csv`; the latter is the gene-list input used by 03_2. Provider attribution for that gene list and `T_cell_signalling.xlsx` is not established here.

**The current `TNBC_6k_metadata_Bogang_selected.RDS` is a locally prepared, finalized metadata object, not the unchanged original file with that name.** Preparation matched the original metadata and Yi's updated metadata by shared `cell_id`, retained the updated fine labels in `celltype_new`, and mapped those labels to the broader `BC6k_refined` classes used downstream. Previous annotations and distances were preserved in `BC6k_refined_old` and `disttonearesttumor_old`. Updated `disttonearestcancer` values were converted from mm to µm during preparation to produce `disttonearesttumor`; distance groups were recalculated. Do not convert this finalized distance again.

This metadata preparation did not regenerate the supplied norm/UMAP/XY objects. Align these objects by cell identifiers, not row position: the finalized metadata snapshot contains 275,204 cells across 172 FOVs and 5 slides, while the normalized expression object contains 294,923 cells by 6,175 genes. Use the recorded SHA-256 identities to distinguish the finalized input from an older file bearing the same name.

The RDS written by 03_2 at `results/03_2/data/03_2_meta_with_immunotype_BW_top_40.rds` is a further **derived analysis output**. Its `immunotype_BW_top_40` labels are calculated locally from tumor-only FOV pseudobulk and the supplied signature gene list; they are not simply copied from Yi's immunotype-reference CSV.

## Fresh-clone quick start

For a collaborator starting from a clean machine:

1. Install **R 4.5.1**. RStudio is recommended.
2. Clone this repository.
3. Supply the RDS inputs and tumor-signature CSV listed above. If using Yi's original data, follow the preparation instructions below; if using finalized metadata, skip preparation.
4. Open `TNBC_TMA_6k_CosMx.Rproj` in RStudio. The tracked `.Rprofile` automatically activates the project renv environment. If using R without RStudio, start R with the repository root as the working directory so the same `.Rprofile` is loaded.
5. Restore the frozen R package environment:

```r
renv::restore()
```

6. If needed, prepare metadata using 02_3 as described below. Place the required GMT resources in `gene_sets/`, then run the reproducibility preflight:

```r
source("scripts/reproducibility_preflight.R")
```

7. After `Preflight PASSED`, run the classification entry point:

```r
source("analysis/03_2_redefine_immunotype_top40_tumorOnly_nonTumorComposition.R")
```

Wait for its final success message and confirm the metadata handoff exists before running dependent branches. Run scripts individually using the dependency order below; filenames are not a requirement to execute every script in sequence.

The renv lockfile does not install R itself, so using R 4.5.1 is the closest match to the audited analysis workstation.

## Prepare metadata from collaborator inputs (before 03_2)

Use this route when you have the original data package and Yi's updated annotation RDS, rather than the finalized metadata. Restore renv first. Place these files in `raw/`:

- The **original, untouched** `TNBC_6k_metadata_Bogang_selected.RDS`.
- Yi's `TNBC_6k_metadata_172FOVs_Bogang.RDS`.
- The unchanged norm, xy and um RDS objects and the tumor-signature CSV listed above.

Then run:

```r
source("analysis/02_3_rebuild_compatible_metadata_from_backup.R")
options(tnbc.rebuilt_metadata = TRUE)
source("scripts/reproducibility_preflight.R")
# After Preflight PASSED:
source("analysis/03_2_redefine_immunotype_top40_tumorOnly_nonTumorComposition.R")
```

02_3 preserves the original metadata at `raw/backup_before_02_2/TNBC_6k_metadata_Bogang_selected_before_02_2.RDS` and uses that immutable backup for reruns. It creates the finalized metadata at the original filename. If you already have this backup, you can supply it directly instead of placing an original file at the target path. Do not rename processed metadata as an original backup. The script rejects duplicate cell IDs, inconsistent FOV/slide identifiers, unknown labels and a cohort other than 275,204 cells / 172 FOVs / 5 slides before replacing the target. Supplied expression and coordinate objects are not modified.

The `tnbc.rebuilt_metadata` option exempts **only the regenerated metadata file** from byte-for-byte SHA-256 matching. It records the observed hash and retains the metadata structure/cohort/coverage checks and hashes for all other inputs. Re-serialization may change file bytes, but matching dimensions alone does not establish identical data or biological results. The exact audited metadata is unavailable in this repository for an independent content comparison. Set this option in each new R session using regenerated metadata. With finalized audited inputs, leave it unset (the strict default).

You do not need to run the superseded 02_1/02_2 transition scripts. If your metadata is already finalized, skip 02_3. Reports from preparation are written to `results/02_3/data/`.

## Project root

Open `TNBC_TMA_6k_CosMx.Rproj` before running analyses. Current reproducibility scripts locate the project by searching upward for this file, so they do not depend on a user-specific drive letter.

## Included analysis workflow

For reproduction from the finalized metadata supplied in `raw/`, the main analysis entry point is:

1. `03_2_redefine_immunotype_top40_tumorOnly_nonTumorComposition.R`
   - tumor-only FOV pseudobulk
   - minimum 30 tumor cells per FOV
   - gene-wise z-score across FOVs
   - ID / IE / IN Top40 signature scoring
   - writes `results/03_2/data/03_2_meta_with_immunotype_BW_top_40.rds`
2. After 03_2, the following 03-series analyses are independent branches and can be run in any order:
   - `03_3_build_FOV_summary_table.R`
   - `03_4_distance_curve_by_immunotype.R`
   - `03_4_1_distance_curve_CD8_subtypes.R`
   - `03_4_2_CD8_counts_and_pies.R`
   - `03_4_3_Tcell_signalling_ssGSEA_by_immunotype.R`
   - `03_4_4_CD8_subtypes_inside_outside_distance.R`
   - `03_5_gene_expression_in_selected_celltype_by_immunotype.R`
   - `03_6_DEG_by_celltype_and_immunotype_all_genes.R`
3. Formal enrichment also branches directly from the 03_2 metadata output plus normalized expression:
   - `04_1_DEG_enrichment_tumorOnly_by_03_2_blockSlide.R` — pairwise tumor contrasts
   - `04_2_GSEA_vsRest_tumorOnly_by_03_2_blockSlide.R` — tumor one-vs-rest contrasts
4. Spatial visualization is independent of the 03_2/04 outputs and reads finalized objects from `raw/`:
   - `05_0_1_plot_all_FOV_celltypes.R`
   - `05_1b_plot_selected_FOV_gene_projection_voronoi.R`

`05_0_check_spatial_plot_inputs.R` is a QC utility for the raw spatial objects.

Before running a branch, review its parameter block for target genes, cell types, thresholds, and FOV selection. In particular, 05_1b defaults to three representative FOVs (`s3f1`, `s4f79`, `s2f2`); it is not an all-FOV cohort analysis. Density filtering can omit cells from spatial panels, so input summary counts may differ from the displayed counts.

Generated tables and PDFs are written beneath each script's `results/` directory. Keep the source commit, session information, parameter files, and console completion/error alongside each run. Rerunning a script may replace outputs; copy previous results before comparison.

The analysis starts at **03_2**, with **02_3** as an optional input-preparation step. The 03 branches are independent after that handoff; spatial scripts read finalized raw inputs directly.

## Preflight before analysis

After placing the core inputs under `raw/`, run:

```r
source("scripts/reproducibility_preflight.R")
```

The preflight checks renv activation, required package availability, required input filenames and SHA-256 identities, the finalized 275,204-cell / 172-FOV / 5-slide metadata snapshot, XY/UM cell-ID coverage, tumor-signature structure, GMT resources, and Cairo PDF capability. It writes a report to `results/preflight/reproducibility_preflight.txt`.

The normalized expression RDS is not loaded by default because it is large. Set `check_large_expression <- TRUE` inside the script for the optional full 294,923 × 6,175 expression-object check.

## R packages

The current workflow uses the following non-base packages across the analysis scripts:

- Core/data/plotting: `data.table`, `Matrix`, `ggplot2`, `scales`, `cowplot`, `deldir`, `matrixStats`, `readxl`
- Reproducibility preflight: `digest`
- Expression/object compatibility: `Seurat`
- Statistics: `limma`, `lme4`, `lmerTest`, `emmeans`
- Gene-set/enrichment: `GSVA`, `BiocParallel`, `SummarizedExperiment`, `clusterProfiler`, `org.Hs.eg.db`, `enrichplot`, `fgsea`

Package versions are frozen in `renv.lock`, generated on the analysis workstation under R 4.5.1. The current lockfile records the CRAN/Bioconductor dependency graph needed by the workflow. See `environment/RENV_SETUP.md` for collaborator restore and maintainer refresh instructions.

## Fixed project conventions

- Current immunotype field: `immunotype_BW_top_40`
- Immunotype order: ID → IE → IN
- Immunotype colors: ID `#000000`, IE `#0000FF`, IN `#FF0000`
- `disttonearesttumor` in the finalized metadata is already in µm; downstream scripts must not multiply it by 1000.
- Broad CD8: `BC6k_refined == "T.cell.CD8"`
- Fine CD8 subtypes: `celltype_new`
- FOV is the biological sample/inference unit; `Run_Tissue_name` is treated as a technical slide/batch variable where modeled.

## Reproducibility notes

- `04_1` and `04_2` include online KEGG calls through clusterProfiler. Those results can depend on the annotation state available at run time. Use the fixed H/C2/C5/C7 GMT version listed above for reproducible local enrichment.
- Generated `results/` are excluded from Git and should be regenerated from the supplied inputs.


### Known enrichment limitations

fgsea may return unestimated NES/p values when gene statistics are unbalanced; these remain NA and must not be treated as nonsignificant. In 04_1/04_2, partially missing scores can also produce a missing descriptive best-subtype label.


No raw data or generated analysis results are included. Syntax checks do not independently validate biological results.


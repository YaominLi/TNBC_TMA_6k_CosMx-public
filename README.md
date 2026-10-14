# TNBC TMA CosMx 6k analysis

R analysis code for the TNBC CosMx 6k project.

## Data layout

The `raw/` directory is intentionally excluded from Git. Supply the six core inputs separately and place them under `raw/` at the repository root without renaming files.

Request these six core files from the project owner before attempting a full run:

| File under `raw/` | Contents / role |
| --- | --- |
| `TNBC_6k_metadata_Bogang_selected.RDS` | Finalized cell annotations, slide/FOV identifiers, and tumor-distance metadata |
| `TNBC_6k_norm_Bogang.RDS` | Normalized expression used for signatures, pseudobulk, and expression analyses |
| `TNBC_6k_xy_Bogang.RDS` | Spatial cell coordinates |
| `TNBC_6k_um_Bogang.RDS` | UMAP embedding coordinates |
| `selected_subtype_specific_genes_tumor.csv` | ID/IE/IN tumor-signature gene lists |
| `T_cell_signalling.xlsx` | T-cell signaling gene sets |

The repository supplies code and environment definitions; downloading it alone does not supply these data. Input identities are recorded in [environment/raw_core_sha256.txt](environment/raw_core_sha256.txt).

Enrichment uses MSigDB **2025.1 human gene-symbol GMT** resources. Obtain the four files listed in [gene_sets/README.md](gene_sets/README.md) before running preflight. GMT resources are not bundled in the public code distribution.

## Fresh-clone quick start

For a collaborator starting from a clean machine:

1. Install **R 4.5.1**. RStudio is recommended.
2. Clone this repository.
3. Copy the six separately supplied core inputs listed above into `raw/` without renaming files.
4. Open `TNBC_TMA_6k_CosMx.Rproj` in RStudio. The tracked `.Rprofile` automatically activates the project renv environment. If using R without RStudio, start R with the repository root as the working directory so the same `.Rprofile` is loaded.
5. Restore the frozen R package environment:

```r
renv::restore()
```

6. Place the required GMT resources in `gene_sets/`, then run the reproducibility preflight:

```r
source("scripts/reproducibility_preflight.R")
```

7. After `Preflight PASSED`, run the classification entry point:

```r
source("analysis/03_2_redefine_immunotype_top40_tumorOnly_nonTumorComposition.R")
```

Wait for its final success message and confirm the metadata handoff exists before running dependent branches. Run scripts individually using the dependency order below; filenames are not a requirement to execute every script in sequence.

The renv lockfile does not install R itself, so using R 4.5.1 is the closest match to the audited analysis workstation.

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

The current workflow starts at **03_2**. The 03 branches are independent after that handoff; spatial scripts read finalized raw inputs directly.

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

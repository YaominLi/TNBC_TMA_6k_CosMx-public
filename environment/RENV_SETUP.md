# R environment reproducibility with renv

The project environment is now frozen in `renv.lock`, generated on the analysis workstation with **R 4.5.1**.

## Collaborator: restore on a clean machine

1. Install R 4.5.1 and clone the repository.
2. Place the six separately supplied core inputs listed in [README.md](../README.md) under `raw/` at the repository root.
3. Open `TNBC_TMA_6k_CosMx.Rproj`. The tracked `.Rprofile` sources `renv/activate.R`, which activates the project library and bootstraps the recorded renv version when needed.
4. Restore package versions:

```r
renv::restore()
```

5. Obtain the four fixed-version GMT files listed in [gene_sets/README.md](../gene_sets/README.md) and place them under `gene_sets/`, then run:

```r
source("scripts/reproducibility_preflight.R")
```

6. After `Preflight PASSED`, begin at `analysis/03_2_redefine_immunotype_top40_tumorOnly_nonTumorComposition.R`.

If `.libPaths()` does not show a path under this project's `renv/library/`, renv is not active and the preflight should not be treated as a valid fresh-clone check.

## Maintainer: refresh the lockfile intentionally

Only refresh the lockfile on the actual analysis workstation after intentionally changing the R/package environment and confirming the workflow still behaves as expected.

```r
renv::status()
renv::snapshot()
```

Then review and commit the reproducibility files:

- `.Rprofile`
- `renv.lock`
- `renv/activate.R`
- `renv/settings.json`
- `renv/.gitignore`

Do **not** commit `renv/library/`, `renv/staging/`, `renv/python/`, or other machine-specific caches.

The current core raw-file SHA-256 identities are recorded in `environment/raw_core_sha256.txt` and are also checked by the preflight script. If a core raw input is intentionally replaced, update that manifest only after re-auditing the changed input and downstream consequences.


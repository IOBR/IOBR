## Resubmission

This is a resubmission. In this version I have addressed all CRAN review comments from Benjamin Altmann regarding functions writing to the user's home filespace (including the package directory and getwd()).

### CRAN Policy Fix: No default writing to user filespace

Per CRAN policy, functions must not write by default or in examples/vignettes/tests to the user's home filespace. All functions now only write when an explicit `path`/`save_path`/`output.dir` argument is provided. When these arguments are `NULL` (the default where appropriate), no files are written to disk.

#### Specific changes

**1. `R/creat_folder.R`**
- Relative paths were resolved against `getwd()`, causing directories to be created in the user's working directory by default.
- Changed to resolve relative paths against `tempdir()` instead, so that `creat_folder()` never creates folders in the user's home or project directory unless an absolute path is given.

**2. `R/eset_distribution.R`**
- Default `project = NULL` caused PNG files to be saved into a `"result"` folder under the working directory.
- Changed behavior: when `project = NULL`, plots are generated but **not saved** to disk. Files are only saved when `project` is explicitly provided.
- Updated documentation and example accordingly.

**3. `R/find_outlier_samples.R`**
- Default `save = TRUE` caused PDF/PNG files to be written to the working directory.
- Changed default to `save = FALSE`.
- Updated documentation and example accordingly.

**4. `R/iobr_cor_plot.R`**
- When `path = NULL`, the function created a subdirectory under `getwd()` using `dir.create()` and `file.path(getwd(), ...)`.
- Changed to only create directories and save files when `path` is explicitly provided (`save_results <- !is.null(path)`). All `ggsave()`, `pdf()`, and `save_pdf()` calls are now conditional on `save_results`.

**5. `R/sig_pheatmap.R`**
- Default `path = NULL` triggered creation of a `"Marker-heatmap-average"` folder under the working directory via `creat_folder()`.
- Changed to only save the heatmap when `path` is explicitly non-`NULL`.
- Updated documentation to state that `path = NULL` means the heatmap is not saved.

**6. `R/sig_box_batch.R`**
- Default `path = NULL` triggered creation of a `"1-sig-box-batch"` folder under the working directory.
- Changed to only save plots when `path` is explicitly provided.
- Updated documentation and example accordingly.

**7. `R/get_cor.R`**
- `.save_cor_plot()` helper used `creat_folder(path)` with a hard-coded default folder name when `path = NULL`.
- Changed to issue a warning and return immediately when `path = NULL`, so no plot or `.RData` file is written.

**8. `R/sig_gsea.R`**
- `.get_gene_sets()` wrote a temporary `"sig.csv"` to the current working directory.
- Changed to use `tempfile()` for the temporary file base name, so the file is created in the system temp directory and properly cleaned up.

**9. `R/find_mutations.R`**
- Multiple unconditional `write.csv()` and `ggsave()` calls wrote files even when `save_path = NULL`.
- Added `!is.null(save_path)` guards around all file-writing operations (CSV exports, individual PDF plots, combined PDF plots, and OncoPrint PDFs).

**10. `R/IPS_calculation.R`**
- When `plot = TRUE`, the function created `"IPS-Results"` under the current working directory.
- Changed to create the folder under `tempdir()` instead.

**11. `R/estimate_helper.R` (`plotPurity`)**
- Default `output.dir = "estimated_purity_plots"` caused PNG files to be saved under the working directory.
- Changed default to `output.dir = NULL`.
- All `dir.create()` and `png()/dev.off()` calls are now conditional on `!is.null(output.dir)`.
- Updated documentation accordingly.

### Examples

All examples that write files now either:
- Use `tempdir()` (e.g., `path = tempdir()`), or
- Run with the default `NULL` path so nothing is written, or
- Are wrapped in `\donttest{}` / `\dontrun{}` where appropriate.

### Test environments

- macOS aarch64, R 4.5.2 (local)

### R CMD check results

```
0 errors | 0 warnings | 0 notes
```

The package passes `R CMD check --as-cran` with no errors, warnings, or notes.

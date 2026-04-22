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
- Added validation: `project` must be provided when `save = TRUE`.
- Updated documentation and example accordingly.

**4. `R/format_signatures.R`**
- Changed `output_name` parameter to `output_path` for clarity.
- Added validation: `output_path` must be provided when `save_signature = TRUE`.
- Updated documentation accordingly.

**5. `R/batch_sig_surv_plot.R`**
- Default `save_path = file.path(tempdir(), "Multiple-KM-plot")` was creating directories.
- Changed default to `save_path = NULL` (no files saved by default).
- Updated documentation accordingly.

**6. `R/iobr_cor_plot.R`**
- When `path = NULL`, the function created a subdirectory under `getwd()` using `dir.create()` and `file.file(getwd(), ...)`.
- Changed to only create directories and save files when `path` is explicitly provided (`save_results <- !is.null(path)`).
- Fixed `ggsave()` calls to use correct parameter order (`filename` first, then `plot`).
- All `ggsave()`, `pdf()`, and `save_pdf()` calls are now conditional on `save_results`.

**7. `R/sig_pheatmap.R`**
- Default `path = NULL` triggered creation of a `"Marker-heatmap-average"` folder under the working directory via `creat_folder()`.
- Changed to only save the heatmap when `path` is explicitly non-`NULL`.
- Updated documentation to state that `path = NULL` means the heatmap is not saved.

**8. `R/sig_box_batch.R`**
- Default `path = NULL` triggered creation of a `"1-sig-box-batch"` folder under the working directory.
- Changed to only save plots when `path` is explicitly provided.
- Updated documentation and example accordingly.

**9. `R/get_cor.R`**
- `.save_cor_plot()` helper used `creat_folder(path)` with a hard-coded default folder name when `path = NULL`.
- Changed to issue a warning and return immediately when `path = NULL`, so no plot or `.RData` file is written.

**10. `R/sig_gsea.R`**
- `.get_gene_sets()` wrote a temporary `"sig.csv"` to the current working directory.
- Changed to use `tempfile()` for the temporary file base name, so the file is created in the system temp directory and properly cleaned up.

**11. `R/find_mutations.R`**
- Fixed semantic naming error: renamed `file_name` variable to `output_dir` (it stores a directory path, not a filename).
- Fixed `ggsave()` calls to use correct parameter order (`filename` first, then `plot`).
- Added `!is.null(save_path)` guards around all file-writing operations (CSV exports, individual PDF plots, combined PDF plots, and OncoPrint PDFs).
- Updated documentation: `save_path` description now correctly states "If NULL, no files are saved".

**12. `R/IPS_calculation.R`**
- When `plot = TRUE`, the function created `"IPS-Results"` under the current working directory.
- Changed to create the folder under `tempdir()` instead.

**13. `R/estimate_helper.R` (`plotPurity`)**
- Default `output.dir = "estimated_purity_plots"` caused PNG files to be saved under the working directory.
- Changed default to `output.dir = NULL`.
- All `dir.create()` and `png()/dev.off()` calls are now conditional on `!is.null(output.dir)`.
- Updated documentation accordingly.

**14. `R/surv_group.R`**
- Fixed `ggsave()` calls to use correct parameter order (`filename` first, then `plot`).

**15. `R/roc_time.R`**
- Fixed `ggsave()` calls to use correct parameter order (`filename` first, then `plot`).

### Code Quality Fixes

Fixed `ggsave()` parameter order issues in multiple files. The correct usage is:
```r
ggplot2::ggsave(filename = "...", plot = p, path = "...")
```
Not:
```r
ggplot2::ggsave(p, filename = "...", path = "...")  # Incorrect
```

### Examples

All examples that write files now either:
- Use `tempdir()` (e.g., `path = tempdir()`), or
- Run with the default `NULL` path so nothing is written, or
- Are wrapped in `\donttest{}` / `\dontrun{}` where appropriate.

### Test environments

- macOS aarch64, R 4.5.2 (local)

### R CMD check results

```
0 errors | 0 warnings | 1 note
```

The only note is "unable to verify current time" which is a network/time server issue unrelated to the package code.

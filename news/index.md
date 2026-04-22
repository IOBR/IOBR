# Changelog

## IOBR 2.2.0

CRAN release: 2026-04-22

### CRAN Policy Fix (Resubmission)

- Fixed all functions that wrote to the user’s home filespace (working
  directory) by default. Writing now only occurs when an explicit `path`
  / `save_path` / `output.dir` is provided.
- Affected functions:
  [`eset_distribution()`](https://iobr.github.io/IOBR/reference/eset_distribution.md),
  [`find_outlier_samples()`](https://iobr.github.io/IOBR/reference/find_outlier_samples.md),
  [`iobr_cor_plot()`](https://iobr.github.io/IOBR/reference/iobr_cor_plot.md),
  [`sig_pheatmap()`](https://iobr.github.io/IOBR/reference/sig_pheatmap.md),
  [`sig_box_batch()`](https://iobr.github.io/IOBR/reference/sig_box_batch.md),
  [`plotPurity()`](https://iobr.github.io/IOBR/reference/plotPurity.md),
  [`IPS_calculation()`](https://iobr.github.io/IOBR/reference/IPS_calculation.md),
  [`find_mutations()`](https://iobr.github.io/IOBR/reference/find_mutations.md),
  [`sig_gsea()`](https://iobr.github.io/IOBR/reference/sig_gsea.md),
  [`get_cor()`](https://iobr.github.io/IOBR/reference/get_cor.md),
  [`batch_sig_surv_plot()`](https://iobr.github.io/IOBR/reference/batch_sig_surv_plot.md),
  [`format_signatures()`](https://iobr.github.io/IOBR/reference/format_signatures.md),
  and
  [`creat_folder()`](https://iobr.github.io/IOBR/reference/creat_folder.md).
- Examples updated to use
  [`tempdir()`](https://rdrr.io/r/base/tempfile.html) when file writing
  is demonstrated.

### Bug Fixes

- **[`find_mutations()`](https://iobr.github.io/IOBR/reference/find_mutations.md)**:
  Fixed semantic naming error where `file_name` variable was used to
  store directory paths. Renamed to `output_dir` for clarity. Fixed
  [`ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html)
  parameter order issues.
- **[`iobr_cor_plot()`](https://iobr.github.io/IOBR/reference/iobr_cor_plot.md)**:
  Fixed
  [`ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html)
  parameter order issues. The correct order is `filename` first, then
  `plot`.
- **[`surv_group()`](https://iobr.github.io/IOBR/reference/surv_group.md)**:
  Fixed
  [`ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html)
  parameter order issues.
- **[`roc_time()`](https://iobr.github.io/IOBR/reference/roc_time.md)**:
  Fixed
  [`ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html)
  parameter order issues.
- **[`batch_sig_surv_plot()`](https://iobr.github.io/IOBR/reference/batch_sig_surv_plot.md)**:
  Changed default `save_path` from
  `file.path(tempdir(), "Multiple-KM-plot")` to `NULL` to prevent
  automatic directory creation.
- **[`format_signatures()`](https://iobr.github.io/IOBR/reference/format_signatures.md)**:
  Changed parameter name from `output_name` to `output_path` for
  consistency. Added validation requiring `output_path` when
  `save_signature = TRUE`.
- **[`find_outlier_samples()`](https://iobr.github.io/IOBR/reference/find_outlier_samples.md)**:
  Added validation requiring `project` when `save = TRUE`.
- **[`plotPurity()`](https://iobr.github.io/IOBR/reference/plotPurity.md)**:
  Changed default `output.dir` from `"estimated_purity_plots"` to
  `NULL`.

### Major Changes

- **CRAN Size Compliance**: Moved large datasets (\>5MB total) from
  `R/sysdata.rda` and `data/` to GitHub Releases to meet CRAN package
  size requirements. Data is now downloaded on-demand and cached
  locally.
- **New Data Management System**:
  - Added
    [`load_data()`](https://iobr.github.io/IOBR/reference/load_data.md)
    function for unified data access (supports sysdata, exported data,
    and GitHub-hosted datasets)
  - Added
    [`download_iobr_data()`](https://iobr.github.io/IOBR/reference/download_iobr_data.md)
    with multiple mirror support (GitHub, ghproxy.vip, gh-proxy.org,
    ghfast.top)
  - Added
    [`add_iobr_mirror()`](https://iobr.github.io/IOBR/reference/add_iobr_mirror.md)
    for custom mirror configuration
  - Added
    [`list_github_datasets()`](https://iobr.github.io/IOBR/reference/list_github_datasets.md)
    to show available remote datasets
  - Added
    [`clear_iobr_cache()`](https://iobr.github.io/IOBR/reference/clear_iobr_cache.md)
    to manage downloaded data
- **Automatic Caching**: Downloaded datasets are cached in user’s cache
  directory and persist across package updates.
- Fixed
  [`count2tpm()`](https://iobr.github.io/IOBR/reference/count2tpm.md) to
  handle symbol-based gene IDs more robustly.

### Data Migration

The following datasets are now hosted on GitHub and downloaded on first
use: - Reference matrices: `BRef`, `TRef`, `lm22` - Annotations:
`anno_gc_vm32`, `anno_grch38`, `anno_hug133plus2`, `anno_illumina`,
`anno_rnaseq` - Example datasets: `tcga_stad_sig`, `imvigor210_sig`,
`eset_stad`, `sig_stad`, `eset_gse62254`, etc. - Gene sets: `hallmark`,
`kegg`, `go_bp`, `go_cc`, `go_mf`, `reactome`, `msig_immune`,
`msig_sc` - Cell markers: `cancer_type_genes`, `cellmarkers`,
`common_genes`, `immuneCuratedData`, `ips_gene_set`, `SI_geneset`,
`mRNA_cell_default`, `mus_human_gene_symbol`, `onco_sig`,
`PurityDataAffy` - Reference data: `xCell.data`, `quantiseq_data` -
Signatures: `signature_collection_citation`, `signature_metabolism`,
`signature_sc`, `signature_tumor`

### Internal Updates

- Updated all internal code to use
  [`load_data()`](https://iobr.github.io/IOBR/reference/load_data.md)
  instead of direct object references for migrated datasets
- Updated `globalVariables.R` to reflect data migration
- Enhanced error messages with manual download instructions when all
  mirrors fail

## IOBR 2.1.1

### New Features

- Added Docker support with automated image builds and publishing to
  GitHub Container Registry (`ghcr.io/iobr/iobr`). The image is based on
  `rocker/tidyverse` and includes IOBR with all dependencies
  pre-installed (#Docker).
- Added Conda/Mamba installation instructions to README for easier
  environment setup.
- Expanded
  [`iobr_cor_plot()`](https://iobr.github.io/IOBR/reference/iobr_cor_plot.md)
  to support multi-group (3+) comparisons using Kruskal-Wallis test, in
  addition to the existing two-group Wilcoxon test
  ([\#28](https://github.com/IOBR/IOBR/issues/28)).

### Bug Fixes

- Fixed
  [`generateRef_DEseq2()`](https://iobr.github.io/IOBR/reference/generateRef_DEseq2.md)
  to handle sparse single-cell RNA-seq data by using
  `type = "poscounts"` for size factor estimation, preventing errors
  when all genes contain at least one zero.
- Fixed
  [`batch_kruskal()`](https://iobr.github.io/IOBR/reference/batch_kruskal.md)
  output format to correctly display mean-centered values for all
  groups.
- Fixed
  [`count2tpm()`](https://iobr.github.io/IOBR/reference/count2tpm.md)
  computation results to align with previous correct versions.
- Updated
  [`count2tpm()`](https://iobr.github.io/IOBR/reference/count2tpm.md)
  documentation to clarify that gene identifiers are converted to gene
  symbols in the output regardless of input ID type.
- Wrapped some long-running examples in `\dontrun{}` and `\donttest{}`
  to comply with CRAN check time limits.

### Improvements

- Refactored internal helper functions in `count2tpm.R` for better code
  organization and maintainability.
- Updated
  [`load_data()`](https://iobr.github.io/IOBR/reference/load_data.md)
  documentation to include additional available datasets.

## IOBR 2.1.0

- Initial CRAN submission.

## Test environments
* local macOS Tahoe (26.5), R 4.5.2
* ubuntu-latest (on GitHub Actions), R-release, R-devel
* windows-latest (on GitHub Actions), R-release
* macos-latest (on GitHub Actions), R-release

## R CMD check results

0 errors | 0 warnings | 0 notes

## Resubmission
This is a resubmission to address CRAN policy violations regarding file system usage and internet access during checks reported by Prof Brian Ripley.

* Changed the default cache directory from `tools::R_user_dir("IOBR", which = "cache")` to a session-specific temporary directory (`file.path(tempdir(), "IOBR_cache")`). This ensures that the package does not write to the user's home directory (`~/.cache/R/IOBR`) by default, complying with CRAN policies.
* Wrapped all examples in documentation that trigger internet downloads (via `download_iobr_data()` or `load_data()` for GitHub-hosted datasets) in `\donttest{}` blocks. This prevents automated checks from downloading large datasets and ensures they fail gracefully or are skipped during the check process.
* Updated `get_iobr_cache_dir()`, `set_iobr_cache_dir()`, and `reset_iobr_cache_dir()` documentation to clarify the new default behavior and provide instructions for users to opt-in to persistent caching.
* Added `Additional_repositories` to `DESCRIPTION` to resolve the `Package required but not available: ‘GSVA’` error on the `r-oldrel-macos-arm64` flavor, ensuring Bioconductor dependencies are correctly discovered by CRAN's automated check tools.
* Regenerated all `.Rd` files to reflect these changes.

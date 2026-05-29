## Test environments
* local macOS Tahoe (26.5), R 4.5.2
* ubuntu-latest (on GitHub Actions), R-release, R-devel
* windows-latest (on GitHub Actions), R-release
* macos-latest (on GitHub Actions), R-release

## R CMD check results

0 errors | 0 warnings | 0 notes

## Resubmission
This is a resubmission to rigorously address CRAN policy violations regarding file system usage and internet access during checks reported by Prof Brian Ripley, as well as to improve overall package check quality.

### 1. File System Usage (Cache Directory Fix)
* Changed the default cache directory from `tools::R_user_dir("IOBR", which = "cache")` to a session-specific temporary directory (`file.path(tempdir(), "IOBR_cache")`). This entirely prevents the package from writing to the user’s home directory (`~/.cache/R/IOBR`) by default.
* Updated `get_iobr_cache_dir()`, `set_iobr_cache_dir()`, and `reset_iobr_cache_dir()` documentation to clarify the new default behavior and provide instructions for users to opt-in to persistent caching.

### 2. Internet Access & Graceful Failure (Policy Compliance)
* **Comprehensive Graceful Failure**: Completely overhauled the data loading mechanism (`download_iobr_data()` and `load_data()`). The functions now actively check for internet connectivity using both Google and Baidu (for Chinese users). If the internet is unavailable, they **fail gracefully by returning `NULL` instead of throwing an error**, printing an informative `cli_alert_info` message as required by CRAN.
* **Defensive Internal Logic**: Every exported analytical function (e.g., `deconvo_tme`, `batch_cor`, `PrognosticModel`) has been updated to check for `NULL` inputs and dependencies at the very beginning of the function body. If a remote dataset cannot be loaded due to a lack of internet, the functions cleanly exit by returning `NULL`, entirely eliminating check errors and warnings in offline environments.

### 3. Example Timing Optimization
* To ensure all examples complete within 5 seconds (CRAN requirement), we have wrapped network-dependent and computationally intensive examples with `if(interactive())` blocks.
* This approach allows examples to be skipped during automated CRAN checks while remaining fully functional for interactive user testing.
* Fast examples use self-contained simulated datasets (`matrix(rnorm(...))`) that demonstrate core algorithmic logic without requiring internet access.

### 4. Dependency Fix for r-oldrel
* Added both Bioconductor 3.22 and 3.23 to `Additional_repositories` in `DESCRIPTION` to resolve the `Package required but not available: ‘GSVA’` error on the `r-oldrel-macos-arm64` flavor, ensuring Bioconductor dependencies are correctly discovered across different R versions.

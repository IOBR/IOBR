## Test environments
* local macOS Tahoe (26.5), R 4.5.2
* ubuntu-latest (on GitHub Actions), R-release, R-devel
* windows-latest (on GitHub Actions), R-release
* macos-latest (on GitHub Actions), R-release

## R CMD check results

0 errors | 0 warnings | 0 notes

## Resubmission
This is a resubmission to address dependency issues reported on `r-oldrel-macos-arm64`.

* Added `BiocManager` to `Suggests` to ensure Bioconductor repositories are correctly identified by automated check tools.
* Fixed GitHub Actions configuration for `r-devel` which was causing false positive dependency resolution errors.

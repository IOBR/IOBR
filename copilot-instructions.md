# IOBR Package Development Instructions

## Development Environment Setup

### R Installation (Using rig + pak - RECOMMENDED)

Use the modern R installation manager rig and package installer pak for optimal development experience:

```bash
# Install rig (R Installation Manager)
curl -Ls https://github.com/r-lib/rig/releases/download/latest/rig-linux-latest.tar.gz | sudo tar xz -C /usr/local

# Install latest R release
sudo rig add release

# Verify installation
rig list
R --version
```

**Timing:** rig installation ~1 minute, R installation ~5-8 minutes

### Package Dependencies Installation

```bash
# Navigate to package directory
cd /path/to/IOBR

# Install all package dependencies using pak (much faster than install.packages)
R --slave -e "pak::local_install_deps()"

# Install additional development tools
R --slave -e "pak::pak(c('devtools', 'knitr', 'rmarkdown', 'testthat', 'pkgdown', 'roxygen2', 'styler', 'lintr'))"
```

**Timing:** Dependencies installation ~4-6 minutes, dev tools ~2-3 minutes

### System Dependencies

Essential system packages (automatically handled by pak):

- pandoc (for vignettes and documentation)
- libxml2-dev, libcurl4-openssl-dev, libssl-dev (for various R packages)
- Graphics libraries: libfreetype6-dev, libjpeg-dev, libpng-dev, libtiff-dev

## Build and Testing Commands

### Package Building

```bash
# Quick build (without vignettes) - recommended for development
R CMD build . --no-build-vignettes

# Full build (with vignettes) - for release
R CMD build .
```

**Timing:** Quick build ~10-20 seconds, Full build ~5-10 minutes

### Package Checking

```bash
# Quick check (skip suggested packages)
_R_CHECK_FORCE_SUGGESTS_=false R CMD check IOBR_*.tar.gz --no-manual

# Full check (requires all suggested packages)
R CMD check IOBR_*.tar.gz --no-manual

# Using devtools (recommended)
R -e "devtools::check()"
```

**Timing:** Quick check ~2-3 minutes, Full check ~10-15 minutes

### Running Tests

```bash
# Install package locally first
R --slave -e "pak::local_install('.')"

# Run tests
R --slave -e "testthat::test_dir('tests/testthat')"

# Test basic functionality
R --slave -e "library(IOBR); data('eset_stad'); print(head(eset_stad))"
```

### Code Formatting and Quality

```bash
# Format all R code
R -e "styler::style_pkg()"

# Check code quality
R -e "lintr::lint_package()"

# Update documentation
R -e "roxygen2::roxygenise()"
```

## Development Workflow

### 1. Initial Setup

```bash
cd /path/to/IOBR
R --slave -e "pak::local_install_deps()"
R --slave -e "pak::pak(c('devtools', 'roxygen2', 'styler', 'lintr'))"
```

### 2. Make Changes

- Edit R files in `R/` directory
- Update documentation in roxygen2 comments
- Add/update tests in `tests/testthat/`

### 3. Format and Check

```bash
# Format code
R -e "styler::style_pkg()"

# Update documentation
R -e "devtools::document()"

# Run checks
R -e "devtools::check()"
```

### 4. Test

```bash
# Install locally
R -e "devtools::install()"

# Run tests
R -e "devtools::test()"
```

### 5. Build and Verify

```bash
# Build package
R CMD build .

# Final check
R CMD check IOBR_*.tar.gz
```

## Common Issues and Solutions

### Issue: Missing dependencies

**Solution:** Install with pak: `R -e "pak::local_install_deps()"`

### Issue: Roxygen version mismatch

**Solution:** Update DESCRIPTION: `RoxygenNote: 7.3.3`

### Issue: Vignette build fails

**Solution:** Install knitr and rmarkdown: `R -e "pak::pak(c('knitr', 'rmarkdown'))"`

### Issue: Check fails with "cannot open connection"

**Solution:** Ensure all system dependencies are installed

## Package Standards

### Code Style

- Use tidyverse style guide
- Run `styler::style_pkg()` before committing
- Maximum line length: 80 characters (documentation can exceed)

### Dependencies

- Minimize dependencies
- Use specific packages instead of meta-packages (e.g., dplyr instead of tidyverse)
- Declare all imports in DESCRIPTION
- Use `@importFrom` for specific functions

### Documentation

- Document all exported functions with roxygen2
- Include examples for all exported functions
- Keep vignettes up-to-date
- Run `devtools::document()` after changes

### Testing

- Write tests for all new functions
- Aim for >80% code coverage
- Run `devtools::test()` before committing

## Quick Reference

```bash
# Complete check workflow
cd /path/to/IOBR
R -e "styler::style_pkg()"
R -e "devtools::document()"
R -e "devtools::check()"
R -e "devtools::test()"

# Build for distribution
R CMD build .
R CMD check IOBR_*.tar.gz --as-cran
```

## Resources

- [R Packages Book](https://r-pkgs.org/)
- [Writing R Extensions](https://cran.r-project.org/doc/manuals/r-release/R-exts.html)
- [rig Documentation](https://github.com/r-lib/rig)
- [pak Documentation](https://pak.r-lib.org/)

#!/bin/bash
# Complete IOBR Package Check Script
# This script performs a full devtools::check() on the IOBR package

set -e  # Exit on error

echo "=== IOBR Package Complete Check Workflow ==="
echo ""

# Step 1: Verify R installation
echo "Step 1: Verifying R installation..."
R --version
echo ""

# Step 2: Install package dependencies
echo "Step 2: Installing package dependencies (this may take 10-15 minutes)..."
cd /home/runner/work/IOBR/IOBR
R --slave -e "pak::local_install_deps()"
echo ""

# Step 3: Install development tools
echo "Step 3: Installing development tools..."
R --slave -e "pak::pak(c('devtools', 'roxygen2', 'styler', 'lintr', 'knitr', 'rmarkdown', 'testthat'))"
echo ""

# Step 4: Document package
echo "Step 4: Updating package documentation..."
R --slave -e "devtools::document()"
echo ""

# Step 5: Install package locally
echo "Step 5: Installing IOBR package locally..."
R --slave -e "devtools::install()"
echo ""

# Step 6: Run devtools::check()
echo "Step 6: Running devtools::check() (this may take 10-15 minutes)..."
R -e "devtools::check(args = c('--no-manual'), error_on = 'warning')"
echo ""

# Step 7: Run tests
echo "Step 7: Running package tests..."
R --slave -e "devtools::test()"
echo ""

# Step 8: Check code quality
echo "Step 8: Checking code quality with lintr..."
R --slave -e "lintr::lint_package()" | head -50
echo ""

# Step 9: Verify package can be loaded
echo "Step 9: Verifying package can be loaded..."
R --slave -e "library(IOBR); cat('IOBR loaded successfully\n')"
echo ""

#  Step 10: Build package
echo "Step 10: Building package tarball..."
R CMD build .
echo ""

# Step 11: Final R CMD check on built package
echo "Step 11: Running R CMD check on built package..."
PACKAGE_FILE=$(ls -t IOBR_*.tar.gz | head -1)
R CMD check "$PACKAGE_FILE" --no-manual
echo ""

echo "=== Package Check Complete ==="
echo ""
echo "Check the results above for any errors, warnings, or notes."
echo "For CRAN submission, aim for 0 errors, 0 warnings, 0 notes."

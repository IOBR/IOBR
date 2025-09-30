# IOBR Roxygen2 Documentation Standards

**Version:** 1.0  
**Date:** September 30, 2025  
**Status:** Applied to all R files

## Overview

This document defines the standardized roxygen2 documentation style for the IOBR package. All R files should follow these conventions to ensure consistent, professional documentation.

## General Structure

### Order of Roxygen Tags

Standard order for roxygen2 documentation blocks:

```r
#' @title (optional - only if different from function name)
#' @description
#' @param parameter1 Description
#' @param parameter2 Description
#' @return Description
#' @author Name
#' @export
#' @import package1
#' @importFrom package2 function1
#' @examples
#' # Example code
```

## Spacing and Formatting

### 1. Space After `#'`
✅ **Correct:**
```r
#' This is a description
```

❌ **Incorrect:**
```r
#'This is a description
```

### 2. Empty Lines
- One empty roxygen line (`#'`) between description and @param
- One empty roxygen line before @examples
- No empty lines between consecutive @param tags

✅ **Correct:**
```r
#' Function description here.
#'
#' @param x Input data
#' @param y Another parameter
#' @return Output description
#'
#' @examples
#' result <- my_function(x, y)
```

### 3. Tag Capitalization
- Use lowercase for all roxygen tags: `@param`, `@return`, `@export`
- Never `@Param`, `@Return`, or `@Export`

## Description Standards

### 1. Function Descriptions
- Start with a capital letter
- End with a period
- Be concise but complete (1-3 sentences)
- Avoid redundant phrases like "This function..."

✅ **Good:**
```r
#' Calculate tumor microenvironment scores using multiple methods.
```

❌ **Avoid:**
```r
#' This function is used to calculate tumor microenvironment scores.
```

### 2. Long Descriptions
- For complex functions, use `@description` tag explicitly
- Break into multiple paragraphs if needed

```r
#' @description
#' Calculate tumor microenvironment scores using multiple deconvolution methods.
#' Supports CIBERSORT, TIMER, xCell, and other algorithms.
#'
#' The function performs normalization, validation, and score calculation
#' with appropriate error handling.
```

## Parameter Documentation

### 1. Parameter Names
- Match exactly with function signature
- Document all parameters (no exceptions)

### 2. Parameter Descriptions
- Start with capital letter
- End with period
- Include type information when helpful
- Specify default values in description

✅ **Good:**
```r
#' @param eset Expression matrix with genes as rows and samples as columns.
#' @param method Deconvolution method to use. Options: "cibersort", "timer", "xcell". Default: "cibersort".
#' @param normalize Logical indicating whether to normalize data. Default: TRUE.
```

### 3. Special Parameter Types
- **Logical parameters:** Specify TRUE/FALSE and meaning
- **Character parameters:** List valid options
- **Numeric parameters:** Mention valid ranges
- **Data frames/matrices:** Describe expected structure

## Return Value Documentation

### 1. Always Document Returns
Every exported function must have `@return` tag

### 2. Return Descriptions
- Describe the structure of return value
- Specify class/type when relevant
- Mention special cases (NULL, errors)

✅ **Good:**
```r
#' @return A data frame with deconvolution scores. Columns represent cell types,
#'   rows represent samples. Returns NULL if deconvolution fails.
```

❌ **Insufficient:**
```r
#' @return Results
```

## Examples

### 1. Always Provide Examples
- Every exported function should have at least one working example
- Use realistic data when possible

### 2. Example Structure
```r
#' @examples
#' # Load example data
#' data(eset_stad, package = "IOBR")
#'
#' # Run analysis
#' result <- my_function(eset_stad, method = "cibersort")
#'
#' # View results
#' head(result)
```

### 3. Long-Running Examples
- Wrap time-consuming examples in `\donttest{}`
- Wrap examples requiring external data in `\dontrun{}`

```r
#' @examples
#' \donttest{
#' # This example requires external data
#' result <- complex_analysis(large_dataset)
#' }
```

## Import Documentation

### 1. Package Imports
- Use `@import` for full package imports (sparingly)
- Prefer `@importFrom package function` for specific functions

✅ **Preferred:**
```r
#' @importFrom dplyr filter mutate select
#' @importFrom ggplot2 ggplot aes geom_point
```

⚠️ **Use sparingly:**
```r
#' @import dplyr
```

### 2. Import Placement
- Place import tags after `@export`
- Group imports by package
- Alphabetize within groups

## Author Information

### 1. Standard Format
```r
#' @author Dongqiang Zeng
```

### 2. Multiple Authors
```r
#' @author Dongqiang Zeng, Author Two, Author Three
```

## Special Cases

### 1. Internal Functions
- Functions not exported should still be documented
- Use descriptive comments even without full roxygen

```r
#' Internal helper function for data validation
#'
#' @param data Input data to validate
#' @return TRUE if valid, FALSE otherwise
#' @keywords internal
internal_validator <- function(data) {
  # Function body
}
```

### 2. S3 Methods
- Document the generic, not each method
- Methods inherit documentation from generic

### 3. Data Objects
```r
#' @title Dataset Name
#' @description Description of the dataset
#' @format A data frame with X rows and Y columns:
#' \describe{
#'   \item{column1}{Description}
#'   \item{column2}{Description}
#' }
#' @source Source of data
"dataset_name"
```

## Common Issues to Avoid

### 1. Missing Documentation
❌ Never leave parameters undocumented
❌ Never omit @return for exported functions
❌ Never skip @examples for exported functions

### 2. Inconsistent Formatting
❌ Mixed spacing after #'
❌ Inconsistent capitalization
❌ Missing periods at end of sentences

### 3. Unhelpful Documentation
❌ "See function" (be specific!)
❌ "Parameter x" (describe what x is!)
❌ Empty @examples blocks

## Validation

### Automated Checks
Run these commands to validate documentation:

```r
# Check for documentation completeness
devtools::document()

# Check for style issues
styler::style_pkg()

# Full package check
devtools::check()
```

### Manual Review Checklist
- [ ] All parameters documented
- [ ] Return value described
- [ ] Examples provided and working
- [ ] Consistent formatting throughout
- [ ] No spelling errors in documentation
- [ ] Appropriate cross-references

## References

- [Writing R Extensions](https://cran.r-project.org/doc/manuals/r-release/R-exts.html)
- [roxygen2 Documentation](https://roxygen2.r-lib.org/)
- [R Packages Book](https://r-pkgs.org/man.html)

---

**Note:** This document serves as the official style guide for IOBR package documentation. All contributors should follow these standards when adding or modifying function documentation.

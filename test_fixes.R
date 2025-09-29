# Test script to validate our fixes
cat("Testing R package fixes...\n")

# Check if DESCRIPTION file has Author field
desc_content <- readLines("DESCRIPTION")
has_author <- any(grepl("^Author:", desc_content))
cat("Has Author field:", has_author, "\n")

# Check if library() calls are removed from key files
r_files <- c("R/zzz.R", "R/LR_cal.R", "R/find_marker_in_bulk.R", 
             "R/find_outlier_samples.R", "R/generateRef_seurat.R", "R/iobr_deg.R")

for (file in r_files) {
  if (file.exists(file)) {
    content <- readLines(file)
    lib_calls <- grep("^[[:space:]]*library\\(", content)
    cat("File", file, "- library() calls:", length(lib_calls), "\n")
  }
}

# Check for global assignments
global_assigns <- system("grep -r '<<-' R/", intern = TRUE)
cat("Global assignments (<<-) found:", length(global_assigns), "\n")

cat("Fixes validation complete.\n")

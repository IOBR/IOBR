#' Stratified Random Sampling of Cells
#'
#' @description
#' Performs stratified random sampling of cells from single-cell data,
#' ensuring proportional representation of each cell type while respecting
#' minimum and maximum count constraints.
#'
#' @param input A data frame or Seurat object containing cell annotations.
#' @param group Character string specifying the column name for cell type
#'   grouping.
#' @param proportion Numeric value between 0 and 1 specifying the sampling
#'   proportion. Default is 0.1.
#' @param minimum_count_include Integer specifying the minimum count threshold
#'   for a cell type to be included in sampling. Default is 300.
#' @param minimum_count Integer specifying the minimum number of cells to
#'   sample per cell type. Default is 200.
#' @param maximum_count Integer specifying the maximum number of cells to
#'   sample per cell type. Default is 1000.
#' @param sub_cluster Optional character string specifying a sub-cluster
#'   column for filtering. Default is NULL.
#' @param cell_type Optional character string specifying the cell type value
#'   to filter when `sub_cluster` is provided. Default is NULL.
#'
#' @return A data frame containing the sampled cells with preserved cell type
#'   proportions.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Sample cells from a data frame
#' sampled_cells <- random_strata_cells(
#'   input = cell_annotations,
#'   group = "cell_type",
#'   proportion = 0.1,
#'   minimum_count_include = 300,
#'   minimum_count = 200,
#'   maximum_count = 1000
#' )
#'
#' # Sample cells from a Seurat object
#' sampled_cells <- random_strata_cells(
#'   input = seurat_object,
#'   group = "seurat_clusters",
#'   proportion = 0.2
#' )
#' }
random_strata_cells <- function(input,
                                group,
                                proportion = 0.1,
                                minimum_count_include = 300,
                                minimum_count = 200,
                                maximum_count = 1000,
                                sub_cluster = NULL,
                                cell_type = NULL) {
  # Validate inputs
  stopifnot(
    "input must be a data frame or Seurat object" =
      is.data.frame(input) || inherits(input, "Seurat"),
    "group must be a character string" = is.character(group) && length(group) == 1,
    "proportion must be between 0 and 1" =
      is.numeric(proportion) && proportion > 0 && proportion < 1,
    "minimum_count_include must be positive" =
      is.numeric(minimum_count_include) && minimum_count_include > 0,
    "minimum_count must be positive" =
      is.numeric(minimum_count) && minimum_count > 0,
    "maximum_count must be greater than minimum_count" =
      is.numeric(maximum_count) && maximum_count > minimum_count
  )

  # Extract metadata from Seurat object
  if (inherits(input, "Seurat")) {
    input <- input@meta.data
  }

  # Check if group column exists
  if (!group %in% colnames(input)) {
    cli::cli_abort("Column '{group}' not found in input data.")
  }

  # Filter invalid entries
  input <- input[!is.na(input[[group]]), ]
  input <- input[!input[[group]] %in% c("Undetermined", "NA", " ", ""), ]

  # Apply sub-cluster filter if specified
  if (!is.null(sub_cluster) && !is.null(cell_type)) {
    if (!sub_cluster %in% colnames(input)) {
      cli::cli_abort("Column '{sub_cluster}' not found in input data.")
    }
    input <- input[input[[sub_cluster]] == cell_type, ]
  }

  # Calculate cell type frequencies
  cell_freq <- as.data.frame(table(input[[group]]))
  colnames(cell_freq) <- c("CellType", "Count")

  cli::cli_alert_info("Cell type counts before sampling:")
  if (interactive()) print(cell_freq)

  # Select cell types meeting minimum count threshold
  valid_types <- as.character(
    cell_freq[cell_freq$Count > minimum_count_include, "CellType"]
  )
  valid_types <- valid_types[!valid_types %in% c("Undetermined", "NA", " ", "")]

  if (length(valid_types) == 0) {
    cli::cli_abort("No cell types meet the minimum count threshold of {minimum_count_include}.")
  }

  cli::cli_alert_info("Cell types included in sampling:")
  if (interactive()) print(valid_types)

  # Filter to valid cell types
  input <- input[input[[group]] %in% valid_types, ]
  input <- input[order(input[[group]]), ]

  # Check minimum count constraint
  filtered_freq <- as.data.frame(table(input[[group]]))
  if (minimum_count > min(filtered_freq$Freq)) {
    cli::cli_abort(c(
      "minimum_count ({minimum_count}) exceeds the smallest cell type count.",
      "i" = "Smallest count: {min(filtered_freq$Freq)}",
      "*" = "Reduce minimum_count or adjust minimum_count_include."
    ))
  }

  # Calculate sample sizes
  sample_sizes <- round(table(input[[group]]) * proportion)

  cli::cli_alert_info("Initial sample sizes (proportion = {proportion}):")
  if (interactive()) print(sample_sizes)

  # Perform stratified sampling
  rlang::check_installed("sampling")
  strata_result <- sampling::strata(
    input,
    group,
    size = as.numeric(sample_sizes),
    method = "srswor"
  )
  sampled_data <- sampling::getdata(input, strata_result)
  sampled_data <- input[rownames(input) %in% rownames(sampled_data), ]

  # Adjust cell types below minimum count
  current_freq <- as.data.frame(table(sampled_data[[group]]))
  colnames(current_freq) <- c("CellType", "Count")

  below_min <- current_freq[current_freq$Count < minimum_count, "CellType"]
  above_max <- current_freq[current_freq$Count > maximum_count, "CellType"]

  # Upsample cell types below minimum
  if (length(below_min) > 0) {
    cli::cli_alert_info("Upsampling cell types below minimum_count ({minimum_count}):")
    print(as.character(below_min))

    for (cell in below_min) {
      cell_data <- input[input[[group]] == cell, ]
      sampled_data <- sampled_data[sampled_data[[group]] != cell, ]
      upsampled <- cell_data[sample(rownames(cell_data), minimum_count), ]
      sampled_data <- rbind(sampled_data, upsampled)
    }
  }

  # Downsample cell types above maximum
  if (length(above_max) > 0) {
    cli::cli_alert_info("Downsampling cell types above maximum_count ({maximum_count}):")
    print(as.character(above_max))

    for (cell in above_max) {
      cell_data <- input[input[[group]] == cell, ]
      sampled_data <- sampled_data[sampled_data[[group]] != cell, ]
      downsampled <- cell_data[sample(rownames(cell_data), maximum_count), ]
      sampled_data <- rbind(sampled_data, downsampled)
    }
  }

  # Report final counts
  final_freq <- as.data.frame(table(sampled_data[[group]]))
  colnames(final_freq) <- c("CellType", "Count")

  cli::cli_alert_info("Cell type counts after sampling:")
  if (interactive()) print(final_freq)

  return(sampled_data)
}

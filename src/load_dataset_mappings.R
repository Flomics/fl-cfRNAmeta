#' Load dataset_mappings.json for plotting and analysis.
#'
#' datasetVisualOrder may include parent placeholder ids (e.g. "giraldez",
#' "block") that group child batches but never appear as sample-level
#' dataset_batch values. Only leaf batches listed in datasetAnalysisBatch
#' are returned as core_order.
#'
#' @param path Path to dataset_mappings.json.
#' @param quiet If FALSE, report filtered parent placeholders.
#' @return A list with mappings, core_order, parent_levels, datasetsLabels,
#'   datasetsPalette, and datasetsMarkers (when present).
load_dataset_mappings <- function(path = "src/dataset_mappings.json", quiet = TRUE) {
  if (!file.exists(path)) {
    stop("Dataset mappings file not found: ", path)
  }

  mappings <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  visual_order <- unlist(mappings$datasetVisualOrder)
  leaf_batches <- names(mappings$datasetAnalysisBatch)

  if (is.null(leaf_batches) || length(leaf_batches) == 0) {
    stop("datasetAnalysisBatch is missing or empty in ", path)
  }

  core_order <- visual_order[visual_order %in% leaf_batches]
  parent_levels <- setdiff(visual_order, leaf_batches)

  unknown_leaves <- setdiff(leaf_batches, visual_order)
  if (length(unknown_leaves) > 0 && !quiet) {
    message(
      "Leaf batches missing from datasetVisualOrder (appended at end): ",
      paste(unknown_leaves, collapse = ", ")
    )
  }
  core_order <- c(core_order, unknown_leaves)

  if (length(parent_levels) > 0 && !quiet) {
    message(
      "Excluded ", length(parent_levels), " parent placeholder level(s) from plot order: ",
      paste(parent_levels, collapse = ", ")
    )
  }

  datasetsLabels <- unlist(mappings$datasetsLabels)
  datasetsPalette <- unlist(mappings$datasetsPalette)
  datasetsMarkers <- if (!is.null(mappings$datasetsMarkers)) {
    unlist(mappings$datasetsMarkers)
  } else {
    NULL
  }

  list(
    mappings = mappings,
    core_order = core_order,
    parent_levels = parent_levels,
    leaf_batches = leaf_batches,
    datasetsLabels = datasetsLabels,
    datasetsPalette = datasetsPalette,
    datasetsMarkers = datasetsMarkers
  )
}

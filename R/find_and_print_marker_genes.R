#' Find Marker Genes for Seurat Clusters (with Optional Parallel Execution)
#'
#' @description
#' Identifies marker genes for each cluster specified in a Seurat object using
#' the `Seurat::FindMarkers` function. It iterates through each cluster, finds
#' markers comparing that cluster against all others, and compiles the results
#' into a single data frame. Supports parallel execution across clusters for
#' improved performance on multi-core machines (works on Windows, macOS, Linux).
#'
#' @details
#' This function sets the default assay and identity (`Idents`) of the Seurat
#' object based on the `assay_name` and `cluster_column` parameters, respectively.
#' It then iterates through each cluster level defined in `Idents(seurat_object)`.
#' For each cluster, it calls `Seurat::FindMarkers` with the specified parameters
#' (`logfc_threshold`, `min_pct`, `only_pos`).
#'
#' Parallel execution, if enabled (`parallel = TRUE`), uses the `parallel` package
#' to create and manage a socket cluster, distributing tasks across cores.
#'
#' The results for all clusters are combined into a single data frame. If an
#' `output_file` path is provided, this data frame is saved as a tab-separated file.
#' If `assign_to_env = TRUE`, the resulting markers data frame is also assigned
#' to the variable `markers_table` in the global environment.
#'
#' @param seurat_object A Seurat object containing single-cell data.
#' @param assay_name Character string. The name of the assay to use for finding
#'   markers (e.g., "RNA", "SCT"). Default is `"RNA"`.
#' @param cluster_column Character string. The name of the metadata column in
#'   `seurat_object@meta.data` that contains the cluster assignments to use for
#'   identity. Default is `"seurat_clusters"`.
#' @param logfc_threshold Numeric. Limit testing to genes which show, on average,
#'   at least X-fold difference (log-scale) between the two groups of cells.
#'   Default is `0.1`.
#' @param min_pct Numeric. Only test genes that are detected in a minimum fraction
#'   of cells in either of the two populations. Default is `0.25`.
#' @param only_pos Logical. If `TRUE`, only return positive markers (genes more
#'   highly expressed in the cluster of interest). Default is `TRUE`.
#' @param parallel Logical. If `TRUE`, attempts to run the marker finding process
#'   in parallel across clusters using a socket cluster from the `parallel` package.
#'   Default is `FALSE`.
#' @param num_cores Optional. Numeric. The number of cores to use if `parallel = TRUE`.
#'   If `NULL` (default) when `parallel = TRUE`, it attempts to use
#'   `parallel::detectCores() - 1` cores, ensuring at least 1 core is used.
#'   Ignored if `parallel = FALSE`.
#' @param output_file Optional. Character string. If provided, the path to save
#'   the resulting marker table as a tab-separated file. If `NULL` (default),
#'   the table is not saved to a file by the function.
#' @param assign_to_env Logical. If `TRUE`, the resulting markers data frame is
#'   assigned to the variable name `markers_table` in the global environment
#'   (`.GlobalEnv`). Default is `FALSE` to avoid cluttering the global namespace.
#'
#' @return A data frame containing marker genes for all clusters, typically ordered
#'   by significance within each cluster group by `Seurat::FindMarkers`. Includes
#'   columns like `p_val`, `avg_log2FC`, `pct.1`, `pct.2`, `p_val_adj`,
#'   `gene` (gene symbol), and `cluster` (the cluster ID). Returns an empty
#'   data frame if no markers are found.
#'
#' @importFrom Seurat DefaultAssay Idents<- Idents levels FindMarkers CreateSeuratObject NormalizeData FindVariableFeatures ScaleData
#' @importFrom utils write.table head
#' @importFrom methods inherits
#' @importFrom dplyr bind_rows
#' @import parallel
#'
#' @seealso \code{\link[Seurat]{FindMarkers}}, \code{\link[parallel]{makeCluster}}, \code{\link[parallel]{parLapply}}, \code{\link{print_markers}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # This example requires the Seurat package.
#' # install.packages("Seurat")
#' # library(Seurat)
#' # library(parallel) # Needed for parallel execution functionality
#' # library(dplyr) # Needed by print_markers helper
#'
#' # Create a dummy Seurat object for demonstration
#' counts <- matrix(rpois(10000, lambda = 1), nrow = 100, ncol = 100)
#' rownames(counts) <- paste0("Gene", 1:100)
#' colnames(counts) <- paste0("Cell", 1:100)
#' seu_obj <- CreateSeuratObject(counts = counts)
#' seu_obj@meta.data$seurat_clusters <- factor(sample(0:2, 100, replace = TRUE))
#' seu_obj <- NormalizeData(seu_obj, verbose = FALSE)
#'
#' # --- Find markers (Sequential, no global assignment) ---
#' marker_results <- find_cluster_markers(
#'   seurat_object = seu_obj,
#'   parallel = FALSE,
#'   assign_to_env = FALSE # Default, explicit here
#' )
#' head(marker_results)
#'
#' # --- Find markers (Parallel, assign to global env) ---
#' # available_cores <- parallel::detectCores()
#' # cores_to_use <- if(available_cores > 1) 2 else 1
#' # find_cluster_markers(
#' #  seurat_object = seu_obj,
#' #  parallel = TRUE,
#' #  num_cores = cores_to_use,
#' #  assign_to_env = TRUE, # Assigns to markers_table globally
#' #  output_file = file.path(tempdir(), "parallel_markers.tsv")
#' # )
#' # Does markers_table exist globally now?
#' # if(exists("markers_table")) head(markers_table)
#' }
#' @examples
#' \dontrun{
#' # Assume 'marker_results' is a data frame from find_cluster_markers()
#' # (See example in find_cluster_markers documentation)
#' # library(dplyr) # Ensure dplyr is loaded if running example standalone
#'
#' # Check if marker_results exists and has data
#' # if (exists("marker_results") && nrow(marker_results) > 0) {
#' #
#' #   # --- Recommended Usage: Pass cluster_id as character string ---
#' #   # This handles simple numeric names ("0", "1") and complex names ("0_1", "T_cells") robustly.
#' #
#' #   # Print top 10 markers for cluster "0"
#' #   print("Top 10 markers for cluster '0':")
#' #   print_markers(marker_table = marker_results, cluster_id = "0", n_genes = 10)
#' #
#' #   # Print top 5 markers for cluster "1"
#' #   print("Top 5 markers for cluster '1':")
#' #   print_markers(marker_table = marker_results, cluster_id = "1", n_genes = 5)
#' #
#' #   # Example if a cluster was named "0_1"
#' #   # print_markers(marker_table = marker_results, cluster_id = "0_1", n_genes = 10)
#' #
#' # } else {
#' #  print("Marker results not available or empty.")
#' # }
#' #
#' # # If find_cluster_markers was run with assign_to_env = TRUE:
#' # if (exists("markers_table") && nrow(markers_table) > 0) {
#' #    print("Top 15 markers for cluster '0' from global markers_table:")
#' #    print_markers(markers_table, cluster_id = "0", n_genes = 15)
#' # }
#' }
#'
find_cluster_markers <- function(seurat_object,
                                 assay_name = "RNA",
                                 cluster_column = "seurat_clusters",
                                 logfc_threshold = 0.1,
                                 min_pct = 0.25,
                                 only_pos = TRUE,
                                 parallel = FALSE,
                                 num_cores = NULL,
                                 output_file = NULL,
                                 assign_to_env = FALSE) { # Added assign_to_env parameter
  
  # --- Input Validation ---
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required but not installed.")
  }
  if (!inherits(seurat_object, "Seurat")) {
    stop("'seurat_object' must be a Seurat object.")
  }
  if (!assay_name %in% names(seurat_object@assays)) {
    stop("Assay '", assay_name, "' not found in the Seurat object.")
  }
  if (!cluster_column %in% names(seurat_object@meta.data)) {
    stop("Cluster column '", cluster_column, "' not found in the Seurat object metadata.")
  }
  if (parallel && !requireNamespace("parallel", quietly = TRUE)) {
    stop("Package 'parallel' is required for parallel execution (parallel = TRUE).")
  }
  # --- End Validation ---
  
  # Set default assay and ident
  Seurat::DefaultAssay(seurat_object) <- assay_name
  Seurat::Idents(seurat_object) <- cluster_column
  
  # Get cluster levels
  if (!is.factor(Seurat::Idents(seurat_object))) {
    warning("Cluster column '", cluster_column, "' is not a factor. Converting to factor.")
    Seurat::Idents(seurat_object) <- factor(Seurat::Idents(seurat_object))
  }
  clusters <- levels(Seurat::Idents(seurat_object))
  if (length(clusters) == 0) {
    stop("No cluster levels found for column '", cluster_column, "'. Cannot find markers.")
  }
  
  message(sprintf("Finding markers for %d clusters: %s", length(clusters), paste(clusters, collapse=", ")))
  message(sprintf("Using Assay: %s, Cluster Column: %s", assay_name, cluster_column))
  message(sprintf("Parameters: logfc_threshold=%.2f, min_pct=%.2f, only_pos=%s",
                  logfc_threshold, min_pct, as.character(only_pos)))
  
  
  # Define the worker function (used by both lapply and parLapply)
  find_markers_for_cluster_worker <- function(cluster_id,
                                              seurat_obj, # Pass object explicitly
                                              logfc_thresh,
                                              min_pct_val,
                                              only_pos_flag) {
    temp_markers <- tryCatch({
      Seurat::FindMarkers(
        seurat_obj, # Use the passed object
        ident.1 = cluster_id,
        logfc.threshold = logfc_thresh,
        min.pct = min_pct_val,
        only.pos = only_pos_flag,
        verbose = FALSE # Keep verbose FALSE
      )
    }, error = function(e) {
      warning("Could not find markers for cluster '", cluster_id, "'. Error: ", e$message, call. = FALSE)
      return(NULL)
    })
    
    if (!is.null(temp_markers) && nrow(temp_markers) > 0) {
      temp_markers$gene <- row.names(temp_markers)
      temp_markers$cluster <- as.character(cluster_id)
      return(temp_markers)
    } else {
      return(NULL)
    }
  }
  
  # --- Execute Marker Finding (Sequential or Parallel) ---
  cluster_markers_list <- NULL
  
  if (parallel) {
    if (!requireNamespace("parallel", quietly = TRUE)) {
      stop("Package 'parallel' needed for parallel = TRUE")
    }
    # --- Parallel Execution ---
    n_cores <- num_cores
    if (is.null(n_cores)) {
      n_cores <- parallel::detectCores() - 1
      n_cores <- max(1, n_cores) # Ensure at least 1 core
    } else {
      n_cores <- max(1, as.integer(num_cores)) # Ensure at least 1 core
    }
    total_cores <- parallel::detectCores()
    if (n_cores > total_cores) {
      warning(sprintf("Requested %d cores, but only %d are available. Using %d cores.", n_cores, total_cores, total_cores))
      n_cores <- total_cores
    }
    
    message(sprintf("Starting parallel execution using %d cores.", n_cores))
    cl <- NULL # Initialize cluster variable
    tryCatch({
      cl <- parallel::makeCluster(n_cores)
      parallel::clusterExport(cl, varlist = c("seurat_object", "logfc_threshold", "min_pct", "only_pos"), envir = environment())
      parallel::clusterEvalQ(cl, { library(Seurat) })
      
      cluster_markers_list <- parallel::parLapply(
        cl = cl, X = clusters, fun = find_markers_for_cluster_worker,
        seurat_obj = seurat_object, logfc_thresh = logfc_threshold,
        min_pct_val = min_pct, only_pos_flag = only_pos
      )
    }, error = function(e) {
      if (!is.null(cl)) parallel::stopCluster(cl)
      stop("Error during parallel execution: ", e$message) # Re-throw error
    }, finally = {
      if (!is.null(cl)) { message("Stopping parallel cluster..."); parallel::stopCluster(cl) }
    })
    message("Parallel execution finished.")
    
  } else {
    # --- Sequential Execution ---
    message("Starting sequential execution...")
    cluster_markers_list <- lapply(
      X = clusters, FUN = find_markers_for_cluster_worker,
      seurat_obj = seurat_object, logfc_thresh = logfc_threshold,
      min_pct_val = min_pct, only_pos_flag = only_pos
    )
    message("Sequential execution finished.")
  }
  
  # --- Process Results ---
  cluster_markers_list <- cluster_markers_list[!sapply(cluster_markers_list, is.null)]
  
  if (length(cluster_markers_list) > 0) {
    if (requireNamespace("dplyr", quietly = TRUE)) {
      markers_table <- dplyr::bind_rows(cluster_markers_list)
    } else {
      message("dplyr not found, using base::do.call(rbind, ...).")
      markers_table <- do.call(rbind, cluster_markers_list)
    }
    message(sprintf("Successfully compiled markers for %d clusters.", length(unique(markers_table$cluster))))
  } else {
    warning("No markers were found for any cluster after processing.")
    markers_table <- data.frame() # Return empty data frame
  }
  
  # --- Save Output File ---
  if (!is.null(output_file)) {
    if (nrow(markers_table) > 0) {
      output_dir <- dirname(output_file)
      if (!dir.exists(output_dir)) { dir.create(output_dir, recursive = TRUE) }
      tryCatch({
        write.table(markers_table, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)
        message("Markers table saved to: ", output_file)
      }, error = function(e) {
        warning("Failed to save markers table to '", output_file, "'. Error: ", e$message)
      })
    } else {
      message("No markers found, output file '", output_file, "' was not created.")
    }
  }
  
  # --- Optional Global Assignment (Side Effect) ---
  if (assign_to_env) {
    assign("markers_table", markers_table, envir = .GlobalEnv)
    message("Note: Markers table assigned to 'markers_table' in the global environment.")
  } else {
    message("Note: Markers table was returned but not assigned to the global environment (assign_to_env=FALSE).")
  }
  
  
  # Return the combined table
  return(markers_table)
}


# ======================================================
# Standalone Helper Function to Print Top Markers
# ======================================================

#' Print Top Marker Genes for a Specific Cluster
#'
#' @description
#' A helper function to quickly print the top marker genes for a specified cluster
#' from a marker table generated by functions like `find_cluster_markers`. It
#' assumes the input table is already reasonably ordered by significance for
#' the relevant cluster (as is typical for `Seurat::FindMarkers` output).
#'
#' @param marker_table A data frame containing marker gene results, typically
#'   output from `find_cluster_markers` or `Seurat::FindMarkers`. Must contain
#'   at least 'gene' and 'cluster' columns.
#' @param cluster_id The specific cluster identifier (matching values in the
#'   'cluster' column of `marker_table`) for which to print markers.
#' @param n_genes Integer. The maximum number of top genes to print for the
#'   specified cluster. Default is 200.
#'
#' @return This function does not return a value (`NULL` implicitly). It prints
#'   the selected gene names to the console, one gene per line.
#'
#' @importFrom dplyr %>% filter select head
#'
#' @export
#'
#' @seealso \code{\link{find_cluster_markers}}
#'
#' @examples
#' \dontrun{
#' # Assume 'marker_results' is a data frame from find_cluster_markers()
#' # (See example in find_cluster_markers documentation)
#'
#' # Check if marker_results exists and has data
#' # if (exists("marker_results") && nrow(marker_results) > 0) {
#' #   # Print top 10 markers for cluster '0' (if cluster '0' exists)
#' #   print_markers(marker_table = marker_results, cluster_id = 0, n_genes = 10)
#' #
#' #   # Print top 5 markers for cluster '1'
#' #   print_markers(marker_table = marker_results, cluster_id = 1, n_genes = 5)
#' # } else {
#' #  print("Marker results not available or empty.")
#' # }
#'
#' # If find_cluster_markers was run with assign_to_env = TRUE:
#' # if (exists("markers_table") && nrow(markers_table) > 0) {
#' #    print_markers(markers_table, cluster_id = 0, n_genes = 15)
#' # }
#' }
#'
print_markers <- function(marker_table, cluster_id, n_genes = 200) {
  # --- Input Validation ---
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required for this function.")
  }
  if (!is.data.frame(marker_table)) {
    stop("'marker_table' must be a data frame.")
  }
  if (!all(c("gene", "cluster") %in% names(marker_table))) {
    stop("'marker_table' must contain 'gene' and 'cluster' columns.")
  }
  # Convert cluster_id to character for safe comparison, as cluster column is character
  cluster_id_char <- as.character(cluster_id)
  available_clusters <- unique(as.character(marker_table$cluster))
  if (!cluster_id_char %in% available_clusters) {
    warning("Cluster '", cluster_id_char, "' not found in the 'cluster' column of the marker table. Available clusters: ", paste(available_clusters, collapse=", "))
    return(invisible(NULL)) # Exit gracefully
  }
  # --- End Validation ---
  
  # Filter, select top N, get gene names
  result_genes <- marker_table %>%
    dplyr::filter(cluster == cluster_id_char) %>% # Use character version for filtering
    # dplyr::arrange(p_val_adj) %>% # Optionally uncomment to explicitly sort by adjusted p-value
    dplyr::head(n = n_genes) %>%
    dplyr::select(gene) %>%
    dplyr::pull(gene) # Extract gene names as a vector
  
  if (length(result_genes) > 0) {
    # Print genes, one per line
    cat(result_genes, sep = "\n")
  } else {
    message("No genes found for cluster '", cluster_id_char, "' after filtering.")
  }
  
  # Return NULL implicitly
  invisible(NULL)
}
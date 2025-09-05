#' Find Cluster Markers with Diagnostic Odds Ratio (DOR)
#'
#' @description
#' Identifies marker genes for each cluster in a Seurat object using
#' `Seurat::FindMarkers`. Optionally calculates the Diagnostic Odds Ratio (DOR)
#' and its logarithm (logDOR) for each marker as an additional measure of
#' specificity. Supports parallel execution and optional saving to text or Excel files.
#'
#' @details
#' This function iterates through each cluster defined by `cluster_column`. For each
#' cluster, it runs `FindMarkers` against all other cells.
#'
#' If `calculate_dor = TRUE`:
#' \itemize{
#'   \item It calculates `DOR = (pct.1 * (1 - pct.2)) / ((1 - pct.1) * pct.2)`.
#'   \item To avoid division by zero or log(0), a small epsilon (1e-9) is added
#'     to `pct.1` and `pct.2` before calculating DOR and `logDOR = log2(DOR)`.
#' }
#'
#' Parallel execution (`parallel = TRUE`) uses the `parallel` package (socket
#' clusters) for cross-platform compatibility.
#'
#' Results can be saved as a single TSV file (`output_file`) and/or an Excel
#' file (`save_excel = TRUE`, requires `WriteXLS` package) with markers for each
#' cluster on a separate sheet. Sheet names are sanitized and made unique.
#'
#' Optionally, the final combined marker table can be assigned to the variable
#' `markers_table` in the global environment (`assign_to_env = TRUE`).
#'
#' @param seurat_object A Seurat object.
#' @param assay_name Character string. Assay to use. Default: `"RNA"`.
#' @param cluster_column Character string. Metadata column with cluster IDs.
#'   Default: `"seurat_clusters"`.
#' @param logfc_threshold Numeric. Log fold-change threshold for `FindMarkers`.
#'   Default: `0.1`.
#' @param min_pct Numeric. Minimum fraction of cells expressing gene in either
#'   population for `FindMarkers`. Default: `0.25`.
#' @param only_pos Logical. Only return positive markers (higher in cluster)?
#'   Default: `TRUE`.
#' @param calculate_dor Logical. Calculate DOR and logDOR? Default: `TRUE`.
#' @param parallel Logical. Use parallel processing? Default: `FALSE`. Requires
#'   the `parallel` package.
#' @param num_cores Optional. Numeric. Number of cores for parallel processing.
#'   Default (if `parallel=TRUE` and `num_cores=NULL`): `detectCores() - 1`.
#' @param output_file Optional. Character string. Path to save the combined marker
#'   table as a tab-separated file. Default: `NULL` (don't save TSV).
#' @param save_excel Logical. Save results to an Excel file with one sheet per
#'   cluster? Default: `FALSE`. Requires the `WriteXLS` package.
#' @param excel_file Optional. Character string. Path for the output Excel file.
#'   Default (if `save_excel=TRUE` and `excel_file=NULL`): `"cluster_markers_dor.xlsx"`
#'   in the current directory.
#' @param assign_to_env Logical. Assign the final results table to `markers_table`
#'   in the global environment? Default: `FALSE`.
#'
#' @return A data frame containing marker genes for all clusters, combined. Includes
#'   standard `FindMarkers` output columns, plus `gene`, `cluster`, and optionally
#'   `DOR` and `logDOR`. Returns an empty data frame if no markers are found.
#'
#' @importFrom Seurat DefaultAssay Idents<- Idents levels FindMarkers CreateSeuratObject NormalizeData
#' @importFrom utils write.table head
#' @importFrom methods requireNamespace inherits
#' @import parallel
#'
#' @seealso \code{\link[Seurat]{FindMarkers}}, \code{\link[WriteXLS]{WriteXLS}}, \code{\link[parallel]{makeCluster}}, \code{\link[parallel]{parLapply}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Requires Seurat, parallel, dplyr (optional), WriteXLS (optional)
#' # install.packages(c("Seurat", "WriteXLS"))
#' # library(Seurat)
#' # library(parallel)
#' # library(dplyr) # Optional, for potential viewing/manipulation
#'
#' # --- Create Dummy Data ---
#' counts <- matrix(rpois(10000, lambda = 1), nrow = 100, ncol = 100)
#' rownames(counts) <- paste0("Gene", 1:100)
#' colnames(counts) <- paste0("Cell", 1:100)
#' seu_obj <- CreateSeuratObject(counts = counts)
#' seu_obj$seurat_clusters <- factor(sample(0:2, 100, replace = TRUE)) # Factor is important
#' seu_obj <- NormalizeData(seu_obj, verbose = FALSE)
#' # --- End Dummy Data ---
#'
#' # --- Run Sequentially, Calculate DOR, Save TSV, Return Table ---
#' markers_dor <- find_cluster_markers_with_dor(
#'   seurat_object = seu_obj,
#'   cluster_column = "seurat_clusters",
#'   calculate_dor = TRUE,
#'   parallel = FALSE,
#'   output_file = file.path(tempdir(), "markers_with_dor.tsv"),
#'   save_excel = FALSE,
#'   assign_to_env = FALSE
#' )
#' head(markers_dor)
#'
#' # --- Run in Parallel, Save Excel, Assign Globally ---
#' # available_cores <- parallel::detectCores()
#' # cores_to_use <- if(available_cores > 1) 2 else 1
#' # find_cluster_markers_with_dor(
#' #   seurat_object = seu_obj,
#' #   cluster_column = "seurat_clusters",
#' #   calculate_dor = TRUE,
#' #   parallel = TRUE,
#' #   num_cores = cores_to_use,
#' #   save_excel = TRUE,
#' #   excel_file = file.path(tempdir(), "markers_dor_parallel.xlsx"),
#' #   assign_to_env = TRUE # Will create 'markers_table' globally
#' # )
#' # Check if 'markers_table' exists globally
#' # if(exists("markers_table")) head(markers_table)
#' # Check if Excel file exists
#' # list.files(tempdir(), pattern = "\\.xlsx$")
#' }
#'
find_cluster_markers_with_dor <- function(seurat_object,
                                          assay_name = "RNA",
                                          cluster_column = "seurat_clusters",
                                          logfc_threshold = 0.1,
                                          min_pct = 0.25,
                                          only_pos = TRUE,
                                          calculate_dor = TRUE,
                                          parallel = FALSE,
                                          num_cores = NULL,
                                          output_file = NULL,
                                          save_excel = FALSE,
                                          excel_file = NULL,
                                          assign_to_env = FALSE) {
  
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
  if (save_excel && !requireNamespace("WriteXLS", quietly = TRUE)) {
    stop("Package 'WriteXLS' is required for Excel output (save_excel = TRUE). Please install it.")
  }
  # --- End Validation ---
  
  # Set default assay and ident
  Seurat::DefaultAssay(seurat_object) <- assay_name
  Seurat::Idents(seurat_object) <- cluster_column
  
  # Get cluster levels (ensure it's a factor)
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
  message(sprintf("Parameters: logfc_threshold=%.2f, min_pct=%.2f, only_pos=%s, calculate_dor=%s",
                  logfc_threshold, min_pct, as.character(only_pos), as.character(calculate_dor)))
  
  # --- Define Worker Function ---
  find_markers_worker <- function(cluster_id, seurat_obj, logfc_thresh, min_pct_val, only_pos_flag, calc_dor_flag) {
    # message(paste0("Starting worker for cluster: ", cluster_id)) # Optional debug message
    temp_markers <- tryCatch({
      Seurat::FindMarkers(
        seurat_obj,
        ident.1 = cluster_id,
        logfc.threshold = logfc_thresh,
        min.pct = min_pct_val,
        only.pos = only_pos_flag,
        verbose = FALSE
      )
    }, error = function(e) {
      warning("Could not find markers for cluster '", cluster_id, "'. Error: ", e$message, call. = FALSE)
      return(NULL)
    })
    
    if (is.null(temp_markers) || nrow(temp_markers) == 0) {
      # message(paste0("No markers/error for cluster: ", cluster_id)) # Optional debug
      return(NULL)
    }
    
    # Calculate DOR if requested and results are valid
    if (calc_dor_flag) {
      # Add epsilon for numerical stability (avoids 0 and 1)
      epsilon <- 1e-9
      pct1_stable <- pmin(pmax(temp_markers$pct.1, epsilon), 1 - epsilon)
      pct2_stable <- pmin(pmax(temp_markers$pct.2, epsilon), 1 - epsilon)
      
      # Calculate DOR using stable percentages
      temp_markers$DOR <- (pct1_stable * (1 - pct2_stable)) / ((1 - pct1_stable) * pct2_stable)
      
      # Calculate log2(DOR) - should now be finite
      temp_markers$logDOR <- log2(temp_markers$DOR)
      
      # Optional: check for any remaining non-finite values (shouldn't happen with epsilon)
      if(any(!is.finite(temp_markers$logDOR))) {
        warning("Non-finite logDOR values encountered for cluster '", cluster_id, "' despite epsilon. Check input pct values.")
        # Replace with NA or handle as needed
        temp_markers$logDOR[!is.finite(temp_markers$logDOR)] <- NA
      }
    }
    
    temp_markers$gene <- rownames(temp_markers)
    temp_markers$cluster <- as.character(cluster_id)
    # message(paste0("Finished worker for cluster: ", cluster_id)) # Optional debug
    return(temp_markers)
  }
  
  # --- Execute Marker Finding (Sequential or Parallel) ---
  cluster_markers_list <- NULL
  
  if (parallel) {
    # --- Parallel Execution ---
    if (!requireNamespace("parallel", quietly = TRUE)) {
      stop("Package 'parallel' needed for parallel = TRUE")
    }
    n_cores <- num_cores
    if (is.null(n_cores)) {
      n_cores <- parallel::detectCores() - 1
      n_cores <- max(1, n_cores)
    } else {
      n_cores <- max(1, as.integer(num_cores))
    }
    total_cores <- parallel::detectCores()
    if (n_cores > total_cores) {
      warning(sprintf("Requested %d cores, but only %d are available. Using %d cores.", n_cores, total_cores, total_cores))
      n_cores <- total_cores
    }
    
    message(sprintf("Starting parallel execution using %d cores.", n_cores))
    cl <- NULL
    tryCatch({
      cl <- parallel::makeCluster(n_cores)
      parallel::clusterExport(cl, varlist = c("seurat_object", "logfc_threshold", "min_pct", "only_pos", "calculate_dor"), envir = environment())
      parallel::clusterEvalQ(cl, { library(Seurat) })
      
      cluster_markers_list <- parallel::parLapply(
        cl = cl, X = clusters, fun = find_markers_worker,
        seurat_obj = seurat_object, logfc_thresh = logfc_threshold,
        min_pct_val = min_pct, only_pos_flag = only_pos, calc_dor_flag = calculate_dor
      )
    }, error = function(e) {
      if (!is.null(cl)) parallel::stopCluster(cl)
      stop("Error during parallel execution: ", e$message)
    }, finally = {
      if (!is.null(cl)) { message("Stopping parallel cluster..."); parallel::stopCluster(cl) }
    })
    message("Parallel execution finished.")
    
  } else {
    # --- Sequential Execution ---
    message("Starting sequential execution...")
    cluster_markers_list <- lapply(
      X = clusters, FUN = find_markers_worker,
      seurat_obj = seurat_object, logfc_thresh = logfc_threshold,
      min_pct_val = min_pct, only_pos_flag = only_pos, calc_dor_flag = calculate_dor
    )
    message("Sequential execution finished.")
  }
  
  # --- Process Results ---
  cluster_markers_list <- cluster_markers_list[!sapply(cluster_markers_list, is.null)] # Remove NULLs
  
  if (length(cluster_markers_list) > 0) {
    # Combine results
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
  
  # --- Save TSV Output File ---
  if (!is.null(output_file)) {
    if (nrow(markers_table) > 0) {
      output_dir <- dirname(output_file)
      if (!dir.exists(output_dir)) { dir.create(output_dir, recursive = TRUE) }
      tryCatch({
        write.table(markers_table, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)
        message("Combined markers table saved to TSV: ", output_file)
      }, error = function(e) {
        warning("Failed to save markers table to TSV '", output_file, "'. Error: ", e$message)
      })
    } else {
      message("No markers found, TSV file '", output_file, "' was not created.")
    }
  }
  
  # --- Save Excel Output File ---
  if (save_excel) {
    if (nrow(markers_table) > 0) {
      # Define default excel file path if not provided
      if (is.null(excel_file)) {
        excel_file <- "cluster_markers_dor.xlsx"
        message("Excel file path not specified, defaulting to: ", excel_file)
      }
      # Ensure directory exists
      excel_dir <- dirname(excel_file)
      if (!dir.exists(excel_dir)) { dir.create(excel_dir, recursive = TRUE) }
      
      message("Preparing Excel file...")
      tryCatch({
        # Split table by cluster for separate sheets
        marker_list_split <- split(markers_table, markers_table$cluster)
        
        # Create safe and unique sheet names
        sheet_names_base <- sapply(names(marker_list_split), function(x) {
          truncated <- substr(x, 1, 31)
          gsub("[^a-zA-Z0-9_.-]", "_", truncated) # Sanitize
        })
        
        sheet_names_final <- make.unique(sheet_names_base, sep = "_")
        # Double check length after make.unique adds separators/numbers
        sheet_names_final <- sapply(sheet_names_final, function(nm) substr(nm, 1, 31))
        # Check again for uniqueness after potential second truncation
        if(any(duplicated(sheet_names_final))) {
          warning("Could not generate unique Excel sheet names under 31 characters after sanitization. Some sheets might be overwritten or cause errors.")
          # Use base names and hope WriteXLS handles it or risk error
          sheet_names_final <- sheet_names_base
        }
        
        
        names(marker_list_split) <- sheet_names_final
        
        # Write to Excel using WriteXLS
        WriteXLS::WriteXLS(marker_list_split, ExcelFileName = excel_file, SheetNames = sheet_names_final)
        message("Excel file saved to: ", excel_file)
        
      }, error = function(e) {
        warning("Failed to save markers to Excel file '", excel_file, "'. Error: ", e$message)
      })
    } else {
      message("No markers found, Excel file was not created.")
    }
  }
  
  
  # --- Optional Global Assignment ---
  if (assign_to_env) {
    assign("markers_table", markers_table, envir = .GlobalEnv)
    message("Note: Markers table assigned to 'markers_table' in the global environment.")
  }
  
  # --- Return Results ---
  return(markers_table)
}

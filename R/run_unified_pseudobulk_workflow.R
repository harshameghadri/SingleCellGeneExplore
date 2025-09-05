#' Unified Workflow for Pseudo-bulk Analysis from Seurat Object
#'
#' This function orchestrates a multi-step workflow:
#' 1. Calculates average gene expression per celltype-subject combination.
#'    (Uses rowMeans by default, with Seurat::AverageExpression as a fallback).
#' 2. Splits this data into cell-type specific long-format tibbles.
#' 3. Applies a user-defined statistical test to these tibbles. If no test
#'    function is provided, it defaults to `Test.wilcox.onAverage_internal`.
#' It provides options to save intermediate and final results in various formats
#' (.Rds, CSV, Excel).
#'
#' @param seurat_obj A Seurat object.
#' @param cell_type_column Character string. Metadata column for cell type identity.
#'        Default is "cell.type.ident".
#' @param subject_column Character string. Metadata column for subject/sample identity.
#'        This column should ideally combine condition and subject ID if they are separate
#'        in the original metadata, e.g., "PAH_Subject1", "CTRL_Subject2", as the
#'        `split_data_by_cell_type2_internal` function expects to parse condition from this.
#'        Default is "orig.ident".
#' @param test_wilcox_fn A user-provided function (or NULL) that takes a single tibble
#'        (as structured by `split_data_by_cell_type2_internal`) and returns a data frame
#'        of statistical results. This function is applied per cell type.
#'        If NULL (default), the internal `Test.wilcox.onAverage_internal` function will be used,
#'        which compares "PAH" vs "CTRL" conditions defined in the 'disease' column
#'        parsed by `split_data_by_cell_type2_internal`.
#'        Default is NULL.
#' @param output_path Character string. Base directory path where all output files
#'        will be saved. If NULL, files are saved in the current working directory.
#'        Default is "./PseudoBulk_Workflow_Outputs".
#'
#' @param save_avg_exp_matrix Logical. Whether to save the wide-format average
#'        expression matrix (output of step 1) as an .Rds file. Default is TRUE.
#' @param avg_exp_matrix_filename Character string. Filename for the saved average
#'        expression matrix .Rds file. Default is "average_expression_matrix.Rds".
#'
#' @param save_cell_type_tibbles Logical. Whether to save the list of cell-type
#'        specific long-format tibbles (output of step 2) as an .Rds file.
#'        Default is TRUE.
#' @param cell_type_tibbles_filename Character string. Filename for the saved
#'        list of tibbles .Rds file. Default is "cell_type_tibbles_list.Rds".
#'
#' @param save_wilcox_rds Logical. Whether to save the final list of Wilcox/statistical
#'        results (output of step 3) as an .Rds file. Default is TRUE.
#' @param wilcox_rds_filename Character string. Filename for the saved Wilcox results
#'        .Rds file. Default is "wilcox_results_list.Rds".
#'
#' @param save_wilcox_csv Logical. Whether to save the combined Wilcox/statistical
#'        results as a single CSV file. Default is TRUE.
#' @param wilcox_csv_filename Character string. Filename for the saved combined
#'        Wilcox results CSV file. Default is "combined_wilcox_results.csv".
#'
#' @param save_wilcox_excel Logical. Whether to save the Wilcox/statistical results
#'        to an Excel file, with each cell type's results on a separate sheet.
#'        Default is TRUE.
#' @param wilcox_excel_filename Character string. Filename for the saved Wilcox
#'        results Excel file. Default is "wilcox_results_by_celltype.xlsx".
#'
#' @param params_for_process_tibbles List. A list of additional parameters to pass
#'        to the `Process_tibbles_Calculate_average_wilcox_internal` function.
#'        Example: `list(user_suffix = "MySuffix", verbose = TRUE, show_progress = TRUE)`.
#'        The `tibble_list` and `test_wilcox_fn` are passed directly by the wrapper.
#'        Default is `list(user_suffix = "APS", verbose = TRUE, show_progress = TRUE)`.
#'
#' @return A list containing the final statistical results (output of step 3),
#'         invisibly. Paths to saved files are printed to the console.
#'
#' @export
#' @examples
#' \dontrun{
#' # Assume 'my_seurat_object' is loaded and preprocessed
#'
#' # Create dummy Seurat object for example
#' if (requireNamespace("Seurat", quietly = TRUE) &&
#'     requireNamespace("SeuratObject", quietly = TRUE) &&
#'     requireNamespace("dplyr", quietly = TRUE) &&
#'     requireNamespace("stats", quietly = TRUE) &&
#'     requireNamespace("tibble", quietly = TRUE)) { 
#'
#'   counts <- matrix(rpois(10000, lambda = 1), nrow = 100, ncol = 100)
#'   rownames(counts) <- paste0("Gene", 1:100)
#'   colnames(counts) <- paste0("Cell", 1:100)
#'   my_seurat_object <- Seurat::CreateSeuratObject(counts = counts)
#'   my_seurat_object$cell.type.ident <- sample(
#'     c("B_cell_long_name_example", "T_cell", "Macrophage_another_very_long_name"),
#'     100, replace = TRUE
#'   )
#'   my_seurat_object$orig.ident <- paste0(
#'     sample(c("PAH", "CTRL"), 100, replace = TRUE), "_",
#'     sample(paste0("S", 1:5), 100, replace = TRUE)
#'   )
#'   my_seurat_object <- Seurat::NormalizeData(my_seurat_object, verbose = FALSE)
#'
#'   # Example 1: Using the default internal Wilcoxon test
#'   default_test_results <- run_unified_pseudobulk_workflow(
#'     seurat_obj = my_seurat_object,
#'     cell_type_column = "cell.type.ident",
#'     subject_column = "orig.ident",
#'     output_path = "./MyWorkflowOutput_DefaultTest_Fallback"
#'   )
#'   if(length(default_test_results) > 0 && !is.null(default_test_results[[1]])) {
#'     print(head(default_test_results[[1]]))
#'   }
#' }
#' }
run_unified_pseudobulk_workflow <- function(
    seurat_obj,
    cell_type_column = "cell.type.ident",
    subject_column = "orig.ident",
    test_wilcox_fn = NULL, 
    output_path = "./PseudoBulk_Workflow_Outputs",
    save_avg_exp_matrix = TRUE,
    avg_exp_matrix_filename = "average_expression_matrix.Rds",
    save_cell_type_tibbles = TRUE,
    cell_type_tibbles_filename = "cell_type_tibbles_list.Rds",
    save_wilcox_rds = TRUE,
    wilcox_rds_filename = "wilcox_results_list.Rds",
    save_wilcox_csv = TRUE,
    wilcox_csv_filename = "combined_wilcox_results.csv",
    save_wilcox_excel = TRUE,
    wilcox_excel_filename = "wilcox_results_by_celltype.xlsx",
    params_for_process_tibbles = list(user_suffix = "APS", verbose = TRUE, show_progress = TRUE)
) {
  
  # --- 0. Preliminaries ---
  if (!requireNamespace("Seurat", quietly = TRUE) || !requireNamespace("SeuratObject", quietly = TRUE)) {
    stop("Packages 'Seurat' and 'SeuratObject' are required.")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' is required.")
  if (!requireNamespace("tidyr", quietly = TRUE)) stop("Package 'tidyr' is required.")
  if (!requireNamespace("purrr", quietly = TRUE)) stop("Package 'purrr' is required (for `map`).")
  if (!requireNamespace("stats", quietly = TRUE)) stop("Package 'stats' is required.")
  if (!requireNamespace("utils", quietly = TRUE)) stop("Package 'utils' is required.")
  if (!requireNamespace("tibble", quietly = TRUE)) stop("Package 'tibble' is required (for rownames_to_column).")
  
  
  verbose_param_val <- function(params_list) {
    val <- params_list$verbose
    if (is.null(val) || !is.logical(val) || length(val) != 1) return(TRUE)
    return(val)
  }
  show_progress_param_val <- function(params_list) {
    val <- params_list$show_progress
    if (is.null(val) || !is.logical(val) || length(val) != 1) return(TRUE)
    return(val)
  }
  
  if (show_progress_param_val(params_for_process_tibbles) &&
      !requireNamespace("progress", quietly = TRUE)) {
    warning("Package 'progress' not found. Progress bars in step 3 will not be shown.")
  }
  if (save_wilcox_excel && !requireNamespace("writexl", quietly = TRUE)) {
    stop("Package 'writexl' is required to save as Excel. Set save_wilcox_excel = FALSE or install 'writexl'.")
  }
  
  actual_test_fn <- NULL
  if (is.null(test_wilcox_fn)) {
    if(verbose_param_val(params_for_process_tibbles)) message("`test_wilcox_fn` is NULL. Using default internal Wilcoxon test: Test.wilcox.onAverage_internal.")
    actual_test_fn <- Test.wilcox.onAverage_internal
  } else if (is.function(test_wilcox_fn)) {
    if(verbose_param_val(params_for_process_tibbles)) message("Using user-provided `test_wilcox_fn`.")
    actual_test_fn <- test_wilcox_fn
  } else {
    stop("`test_wilcox_fn` must be a function or NULL.")
  }
  
  if (!is.null(output_path)) {
    if (!dir.exists(output_path)) {
      if(verbose_param_val(params_for_process_tibbles)) message("Output path '", output_path, "' does not exist. Creating it.")
      dir.create(output_path, recursive = TRUE)
    }
  } else {
    output_path <- "." 
  }
  
  message_verbose_main <- function(msg) {
    if(verbose_param_val(params_for_process_tibbles)) message(msg)
  }
  
  # --- 1. Calculate Average Expression (with fallback) ---
  message_verbose_main("Step 1: Calculating average expression matrix...")
  avg_exp_matrix <- calculate_average_expression_internal(
    seurat_obj = seurat_obj,
    cell_type_column = cell_type_column,
    subject_column = subject_column,
    verbose_main = verbose_param_val(params_for_process_tibbles) # Pass verbose flag
  )
  
  if (save_avg_exp_matrix) {
    file_path <- file.path(output_path, avg_exp_matrix_filename)
    saveRDS(avg_exp_matrix, file = file_path)
    message_verbose_main(paste("  - Average expression matrix saved to:", file_path))
  }
  message_verbose_main("Step 1: Completed.\n")
  
  # --- 2. Split Data by Cell Type into Tibbles ---
  message_verbose_main("Step 2: Splitting data by cell type into tibbles...")
  cell_type_tibbles_list <- split_data_by_cell_type2_internal(df = avg_exp_matrix)
  
  if (save_cell_type_tibbles) {
    file_path <- file.path(output_path, cell_type_tibbles_filename)
    saveRDS(cell_type_tibbles_list, file = file_path)
    message_verbose_main(paste("  - List of cell type tibbles saved to:", file_path))
  }
  message_verbose_main("Step 2: Completed.\n")
  
  # --- 3. Process Tibbles & Calculate Wilcox (or other stats) ---
  message_verbose_main("Step 3: Processing tibbles and applying statistical test...")
  process_args <- list(
    tibble_list = cell_type_tibbles_list,
    test_wilcox_fn = actual_test_fn
  )
  user_params <- params_for_process_tibbles
  user_params$save_path <- NULL
  user_params$save_name <- NULL
  default_process_params <- list(user_suffix = "APS", verbose = TRUE, show_progress = TRUE)
  for (p_name in names(user_params)) {
    if (!is.null(user_params[[p_name]])) {
      default_process_params[[p_name]] <- user_params[[p_name]]
    }
  }
  for (p_name in names(default_process_params)) {
    process_args[[p_name]] <- default_process_params[[p_name]]
  }
  wilcox_results_list <- do.call(Process_tibbles_Calculate_average_wilcox_internal, process_args)
  
  if (save_wilcox_rds) {
    file_path <- file.path(output_path, wilcox_rds_filename)
    saveRDS(wilcox_results_list, file = file_path)
    message_verbose_main(paste("  - List of Wilcox/statistical results saved to (.Rds):", file_path))
  }
  
  if (save_wilcox_csv) {
    if (length(wilcox_results_list) > 0 && sum(sapply(wilcox_results_list, function(x) !is.null(x) && nrow(x) > 0)) > 0) {
      valid_results <- Filter(function(x) !is.null(x) && is.data.frame(x) && nrow(x) > 0, wilcox_results_list)
      if(length(valid_results) > 0) {
        combined_df <- dplyr::bind_rows(valid_results, .id = "cell_type_source_list_name")
        if("cell_type_source_list_name" %in% colnames(combined_df) && !is.null(process_args$user_suffix)){
          suffix_to_remove <- paste0(".", process_args$user_suffix, ".Average")
          combined_df$cell_type <- gsub(suffix_to_remove, "", combined_df$cell_type_source_list_name, fixed = TRUE)
          combined_df <- dplyr::select(combined_df, -"cell_type_source_list_name")
          combined_df <- dplyr::relocate(combined_df, cell_type, .before = 1)
        }
        file_path <- file.path(output_path, wilcox_csv_filename)
        utils::write.csv(combined_df, file = file_path, row.names = FALSE)
        message_verbose_main(paste("  - Combined Wilcox/statistical results saved to (.csv):", file_path))
      } else { message_verbose_main("  - No valid data frames in results list to combine. Skipping CSV save.") }
    } else { message_verbose_main("  - Results list is empty or contains no data. Skipping CSV save.") }
  }
  
  if (save_wilcox_excel) {
    if (length(wilcox_results_list) > 0 && sum(sapply(wilcox_results_list, function(x) !is.null(x) && nrow(x) > 0)) > 0) {
      valid_results_for_excel <- Filter(function(x) !is.null(x) && is.data.frame(x) && nrow(x) > 0, wilcox_results_list)
      if(length(valid_results_for_excel) > 0) {
        original_sheet_names <- names(valid_results_for_excel)
        sanitized_sheet_names <- sapply(original_sheet_names, function(name) {
          base_name <- name
          if(!is.null(process_args$user_suffix)){
            suffix_to_remove <- paste0(".", process_args$user_suffix, ".Average")
            base_name <- gsub(suffix_to_remove, "", name, fixed = TRUE)
          }
          clean_name <- gsub("[*:/\\\\?\\[\\]]", "", base_name)
          truncated_name <- substr(clean_name, 1, 30)
          return(truncated_name)
        })
        if (any(duplicated(sanitized_sheet_names))) {
          sanitized_sheet_names <- make.unique(sanitized_sheet_names, sep = "_")
          sanitized_sheet_names <- sapply(sanitized_sheet_names, function(name) substr(name, 1, 31))
        }
        if (any(duplicated(sanitized_sheet_names))) {
          sanitized_sheet_names <- make.unique(sapply(original_sheet_names, function(name) {
            base_name <- name; if(!is.null(process_args$user_suffix)){ suffix_to_remove <- paste0(".", process_args$user_suffix, ".Average"); base_name <- gsub(suffix_to_remove, "", name, fixed = TRUE) }
            substr(gsub("[*:/\\\\?\\[\\]]", "", base_name),1,28)
          }), sep="_")
          message_verbose_main("  - Complex duplicate sheet names encountered, used simpler unique naming for Excel.")
        }
        names(valid_results_for_excel) <- sanitized_sheet_names
        file_path <- file.path(output_path, wilcox_excel_filename)
        writexl::write_xlsx(valid_results_for_excel, path = file_path)
        message_verbose_main(paste("  - Wilcox/statistical results saved to (Excel):", file_path))
      } else { message_verbose_main("  - No valid data frames in results list for Excel. Skipping Excel save.") }
    } else { message_verbose_main("  - Results list is empty or contains no data. Skipping Excel save.") }
  }
  message_verbose_main("Step 3: Completed.\n")
  message_verbose_main("Unified workflow finished successfully!")
  return(invisible(wilcox_results_list))
}


# --- Internal Helper Functions ---

calculate_average_expression_internal <- function(seurat_obj,
                                                  cell_type_column = "cell.type.ident",
                                                  subject_column = "orig.ident",
                                                  verbose_main = TRUE) { # Added verbose_main for internal messages
  
  # Helper for messaging within this function
  message_verbose_calc <- function(msg) { if(verbose_main) message(msg) }
  
  results_df <- NULL
  primary_method_failed <- FALSE
  
  # --- Attempt 1: Primary method (current rowMeans approach) ---
  tryCatch({
    message_verbose_calc("  Attempting average expression with primary (rowMeans) method...")
    
    if (any(duplicated(colnames(seurat_obj)))) {
      warning("Duplicate cell barcodes (colnames) found in Seurat object. This may lead to incorrect aggregation.", call. = FALSE)
    }
    seurat_obj$cellBarcode <- colnames(seurat_obj)
    
    if(!cell_type_column %in% colnames(seurat_obj@meta.data)) stop(paste("`cell_type_column`:", cell_type_column, "not found."))
    if(!subject_column %in% colnames(seurat_obj@meta.data)) stop(paste("`subject_column`:", subject_column, "not found."))
    
    metadata <- seurat_obj@meta.data[, c(cell_type_column, subject_column, "cellBarcode"), drop = FALSE]
    celltype_subject_combinations <- metadata %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(c(cell_type_column, subject_column)))) %>%
      dplyr::summarise(cellBarcodes = list(unique(as.character(cellBarcode))), .groups = "drop")
    
    if (nrow(celltype_subject_combinations) == 0) {
      stop("Primary method: No cell type/subject combinations found.")
    }
    
    default_assay <- Seurat::DefaultAssay(seurat_obj)
    if (is.null(default_assay) || !default_assay %in% SeuratObject::Assays(seurat_obj)) stop("Default assay not found.")
    
    assay_obj <- seurat_obj[[default_assay]]
    if ("data" %in% SeuratObject::Layers(assay_obj)) {
      assay_data <- SeuratObject::GetAssayData(assay_obj, layer = "data")
    } else if ("data" %in% methods::slotNames(assay_obj) && !is.null(methods::slot(assay_obj, "data")) && nrow(methods::slot(assay_obj, "data")) > 0) {
      assay_data <- SeuratObject::GetAssayData(assay_obj, slot = "data")
    } else { stop(paste0("Primary method: 'data' slot/layer in '",default_assay,"' not available/empty.")) }
    
    if (nrow(assay_data) == 0 || ncol(assay_data) == 0) stop(paste0("Primary method: 'data' slot/layer in '",default_assay,"' is empty."))
    gene_names <- rownames(assay_data)
    
    average_expression_list_df <- celltype_subject_combinations %>%
      dplyr::mutate(average_expression = purrr::map(cellBarcodes, ~ {
        cells_to_avg <- .x
        valid_cells <- intersect(cells_to_avg, colnames(assay_data))
        current_group_df <- celltype_subject_combinations[purrr::map_lgl(cellBarcodes, function(cb_list) identical(sort(cb_list), sort(.x))), ]
        ct_info <- if(nrow(current_group_df) > 0) current_group_df[[1, cell_type_column]] else "UnknownCT"
        subj_info <- if(nrow(current_group_df) > 0) current_group_df[[1, subject_column]] else "UnknownSubj"
        if (length(valid_cells) == 0) {
          warning(paste("Primary method: No valid cells for cell_type '", ct_info, "', subject '", subj_info, "'. NAs returned."), call. = FALSE)
          return(rep(NA_real_, length(gene_names)))
        }
        subset_assay_data <- assay_data[, valid_cells, drop = FALSE]
        if (!is.matrix(subset_assay_data) && !is(subset_assay_data, "dgCMatrix")) subset_assay_data <- as.matrix(subset_assay_data)
        if (ncol(subset_assay_data) > 1) rowMeans(subset_assay_data, na.rm = TRUE)
        else if (ncol(subset_assay_data) == 1) as.numeric(subset_assay_data[,1])
        else rep(NA_real_, length(gene_names))
      }))
    
    results_primary <- data.frame(gene = gene_names)
    if (nrow(results_primary) == 0) stop("Primary method: No genes found.")
    
    for (i in 1:nrow(celltype_subject_combinations)) {
      ct_val <- celltype_subject_combinations[[i, cell_type_column]]; subj_val <- celltype_subject_combinations[[i, subject_column]]
      ct_val_clean <- ifelse(is.na(ct_val), "NA_celltype", as.character(ct_val)); subj_val_clean <- ifelse(is.na(subj_val), "NA_subject", as.character(subj_val))
      col_name <- make.names(paste(ct_val_clean, subj_val_clean, sep = "_"), unique = FALSE)
      current_avg_expr <- average_expression_list_df$average_expression[[i]]
      if (length(current_avg_expr) == nrow(results_primary)) results_primary[[col_name]] <- current_avg_expr
      else { results_primary[[col_name]] <- rep(NA_real_, nrow(results_primary)); warning(paste("Primary method: Length mismatch for", col_name), call. = FALSE)}
    }
    
    if (is.null(results_primary) || ncol(results_primary) <= 1) stop("Primary method yielded empty or invalid results (no group columns).") # ncol<=1 means only 'gene' column
    results_df <- results_primary
    message_verbose_calc("  Primary (rowMeans) method for average expression succeeded.")
    
  }, error = function(e_primary) {
    warning(paste("Primary method (rowMeans) for average expression failed: ", e_primary$message, 
                  "\n  Attempting fallback with Seurat::AverageExpression()."), call. = FALSE)
    primary_method_failed <<- TRUE # Assign to variable in parent env of tryCatch
  })
  
  # --- Attempt 2: Fallback method (Seurat::AverageExpression) if primary failed ---
  if (primary_method_failed || is.null(results_df)) {
    tryCatch({
      message_verbose_calc("  Attempting average expression with fallback (Seurat::AverageExpression) method...")
      default_assay_name <- Seurat::DefaultAssay(seurat_obj)
      assay_obj_fb <- seurat_obj[[default_assay_name]]
      
      features_to_avg <- NULL
      if ("data" %in% SeuratObject::Layers(assay_obj_fb)) {
        features_to_avg <- rownames(SeuratObject::GetAssayData(assay_obj_fb, layer = "data"))
      } else if ("data" %in% methods::slotNames(assay_obj_fb) && !is.null(methods::slot(assay_obj_fb, "data")) && nrow(methods::slot(assay_obj_fb, "data")) > 0) {
        features_to_avg <- rownames(SeuratObject::GetAssayData(assay_obj_fb, slot = "data"))
      } else { stop("Fallback: 'data' slot/layer not available/empty.") }
      
      if (length(features_to_avg) == 0) stop("Fallback: No features to average.")
      
      avg_exp_data_list <- Seurat::AverageExpression(
        seurat_obj,
        features = features_to_avg,
        group.by = c(cell_type_column, subject_column), # Groups by these two factors
        slot = "data", 
        assay = default_assay_name 
      )
      
      avg_exp_matrix_from_seurat <- avg_exp_data_list[[default_assay_name]]
      if (is.null(avg_exp_matrix_from_seurat) || nrow(avg_exp_matrix_from_seurat) == 0 || ncol(avg_exp_matrix_from_seurat) == 0) {
        stop("Fallback method (Seurat::AverageExpression) yielded empty or invalid matrix.")
      }
      
      results_fallback <- as.data.frame(avg_exp_matrix_from_seurat)
      results_fallback <- tibble::rownames_to_column(results_fallback, var = "gene")
      
      # Adjust column names: Seurat::AverageExpression with group.by = c("A", "B")
      # creates colnames like "ValFromA.ValFromB".
      # Our original method creates "ValFromA_ValFromB".
      # We need to ensure consistency for split_data_by_cell_type2_internal.
      # The subject_column already contains "Disease_SubjectID".
      # So AverageExpression columns will be CellType.Disease_SubjectID.
      # We need CellType_Disease_SubjectID.
      
      # Reconstruct column names to match primary method's expected output for split_data_by_cell_type2
      # Original col names from AverageExpression are like: "B_cell.PAH_S1"
      # We need "B_cell_PAH_S1"
      # This assumes subject_column has the Disease_SubjectID format.
      # And cell_type_column is the pure cell type.
      
      # Get the original combinations that formed the columns in AverageExpression
      # This is tricky as AverageExpression doesn't directly return this mapping easily.
      # We know the groups are combinations of unique values from cell_type_column and subject_column.
      # The order of columns from AverageExpression should correspond to unique(interaction(meta[[c1]], meta[[c2]]))
      # For simplicity and robustness, if column names from AverageExpression use '.', replace with '_'.
      # This should generally work if the structure is CellType.Disease_SubjectID.
      
      new_colnames <- colnames(results_fallback)
      # Replace the first dot (typically separating cell_type from subject_column content) with an underscore.
      # This is a common pattern from AverageExpression when group.by has multiple elements.
      # Example: Bcell.PAH_S1 -> Bcell_PAH_S1
      new_colnames <- sub(pattern = "\\.", replacement = "_", x = new_colnames)
      colnames(results_fallback) <- new_colnames
      # Ensure 'gene' column is still named 'gene'
      if (colnames(results_fallback)[1] != "gene" && "gene" %in% colnames(results_fallback)) {
        # if 'gene' is not first, try to find it and reorder, or rename if first col is 'gene' but with different case
        if(tolower(colnames(results_fallback)[1]) == "gene") colnames(results_fallback)[1] <- "gene"
      }
      
      
      results_df <- results_fallback
      message_verbose_calc("  Fallback (Seurat::AverageExpression) method for average expression succeeded.")
    }, error = function(e_fallback) {
      stop(paste0("Both primary (rowMeans) and fallback (Seurat::AverageExpression) methods failed to calculate average expression. ",
                  "Please check your Seurat object, `cell_type_column`, and `subject_column`.\n",
                  "Primary error was: (see warning above if one occurred)\n",
                  "Fallback error: ", e_fallback$message), call. = FALSE)
    })
  }
  
  if (is.null(results_df) || ncol(results_df) <= 1) { # ncol <=1 implies only 'gene' col or empty
    stop("Failed to calculate and format average expression using any method, or result is empty.", call. = FALSE)
  }
  return(results_df)
}

Test.wilcox.onAverage_internal <- function(data) {
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame or tibble for Test.wilcox.onAverage_internal.")
  }
  required_cols <- c("gene", "average_expression", "disease")
  if (!all(required_cols %in% names(data))) {
    stop(sprintf("'data' must contain columns: %s for Test.wilcox.onAverage_internal. Missing: %s",
                 paste(required_cols, collapse = ", "),
                 paste(setdiff(required_cols, names(data)), collapse = ", ")))
  }
  current_cell_type_info <- if("cell_type" %in% names(data) && length(unique(data$cell_type)) == 1) {
    paste0("for cell type '", unique(data$cell_type), "'")
  } else {""}
  
  if (!any("PAH" %in% data$disease) || !any("CTRL" %in% data$disease)) {
    warning(paste0("Input data ", current_cell_type_info,
                   " does not contain both 'PAH' and 'CTRL' values in the 'disease' column for comparison. ",
                   "Results for relevant genes will likely be NA or Inf."), call. = FALSE)
  }
  
  result <- data %>%
    dplyr::group_by(.data$gene) %>%
    dplyr::summarise(
      Average_PAH = mean(.data$average_expression[.data$disease == "PAH"], na.rm = TRUE),
      Average_CTRL = mean(.data$average_expression[.data$disease == "CTRL"], na.rm = TRUE),
      log2FoldChange = dplyr::case_when(
        !is.nan(Average_PAH) & Average_PAH > 0 & !is.nan(Average_CTRL) & Average_CTRL > 0 ~ log2(Average_PAH / Average_CTRL),
        !is.nan(Average_PAH) & Average_PAH > 0 & (is.nan(Average_CTRL) | Average_CTRL == 0) ~ Inf,
        (is.nan(Average_PAH) | Average_PAH == 0) & !is.nan(Average_CTRL) & Average_CTRL > 0 ~ -Inf,
        TRUE ~ NaN
      ),
      p_value = {
        pah_values <- .data$average_expression[.data$disease == "PAH" & !is.na(.data$average_expression)]
        ctrl_values <- .data$average_expression[.data$disease == "CTRL" & !is.na(.data$average_expression)]
        if (length(pah_values) >= 1 && length(ctrl_values) >= 1) {
          tryCatch({
            if (length(unique(c(pah_values, ctrl_values))) <= 1) NA_real_
            else stats::wilcox.test(pah_values, ctrl_values, exact = FALSE)$p.value
          }, error = function(e) NA_real_)
        } else NA_real_
      },
      .groups = 'drop'
    ) %>%
    dplyr::mutate(
      adjusted_p_value.BH = stats::p.adjust(.data$p_value, method = "BH"),
      adjusted_p_value.FDR = stats::p.adjust(.data$p_value, method = "fdr"),
      adjusted_p_value.Bonferroni = stats::p.adjust(.data$p_value, method = "bonferroni")
    ) %>%
    dplyr::filter(!is.na(.data$log2FoldChange) & !is.nan(.data$log2FoldChange) & is.finite(.data$log2FoldChange))
  return(result)
}

split_data_by_cell_type2_internal <- function(df) {
  if (!"gene" %in% colnames(df)) stop("Input to split_data_by_cell_type2_internal must contain 'gene' column.", call. = FALSE)
  all_cols <- colnames(df)[colnames(df) != "gene"]
  if (length(all_cols) == 0) {
    warning("No data columns (besides 'gene') in input to split_data_by_cell_type2_internal. Empty list returned.", call. = FALSE)
    return(list())
  }
  col_info <- data.frame(column = all_cols, stringsAsFactors = FALSE)
  parsed_cols <- strsplit(col_info$column, "_", fixed = TRUE)
  col_info$cell_type <- sapply(parsed_cols, function(x) if(length(x) >= 1) x[1] else NA_character_)
  col_info$disease   <- sapply(parsed_cols, function(x) if(length(x) >= 2) x[2] else NA_character_)
  col_info$subject   <- sapply(parsed_cols, function(x) {
    if (length(x) >= 3) paste(x[3:length(x)], collapse = "_")
    else if (length(x) == 2) "NoSubjectID_Parsed" 
    else NA_character_
  })
  if(all(is.na(col_info$cell_type))) warning("Could not parse 'cell_type' from any column names in split_data. Example col: ", all_cols[1], call. = FALSE)
  if(all(is.na(col_info$disease))) warning("Could not parse 'disease' from any column names in split_data. Example col: ", all_cols[1], call. = FALSE)
  
  cell_type_data <- list()
  unique_cell_types <- unique(stats::na.omit(col_info$cell_type))
  if (length(unique_cell_types) == 0) {
    warning("No unique cell types determined after parsing column names in split_data. Empty list returned.", call. = FALSE)
    return(list())
  }
  for (ct in unique_cell_types) {
    ct_cols_to_select <- col_info$column[which(col_info$cell_type == ct & !is.na(col_info$cell_type))]
    ct_info_subset <- col_info[col_info$column %in% ct_cols_to_select, ]
    if (length(ct_cols_to_select) == 0) next
    ct_df_wide <- df %>% dplyr::select(dplyr::all_of(c("gene", ct_cols_to_select)))
    ct_long <- ct_df_wide %>%
      tidyr::pivot_longer(cols = -dplyr::all_of("gene"), names_to = "column", values_to = "average_expression") %>%
      dplyr::left_join(ct_info_subset, by = "column") %>%
      dplyr::select(gene, cell_type, disease, subject, average_expression)
    cell_type_data[[ct]] <- ct_long
  }
  return(cell_type_data)
}

Process_tibbles_Calculate_average_wilcox_internal <- function(
    tibble_list, test_wilcox_fn, user_suffix = "APS", verbose = TRUE, show_progress = TRUE
) {
  if (!is.list(tibble_list)) stop("`tibble_list` must be a list.")
  if (is.null(test_wilcox_fn) || !is.function(test_wilcox_fn)) stop("`test_wilcox_fn` must be a function.")
  verbose_msg_internal <- function(...) if (verbose) message(...)
  pb_overall <- NULL; pb_test_application <- NULL
  if (show_progress && requireNamespace("progress", quietly = TRUE)) {
    pb_overall <- progress::progress_bar$new(format = "Overall Step 3 Progress [:bar] :percent eta: :eta", total = 1, clear = FALSE, width = 80)
    pb_test_application <- progress::progress_bar$new(format = "Applying Stat Test [:bar] :percent eta: :eta", total = length(tibble_list), clear = FALSE, width = 80)
  } else if (show_progress) { verbose_msg_internal("Package 'progress' not available. Disabling progress bars."); show_progress <- FALSE }
  verbose_msg_internal("\nApplying statistical test to each cell type tibble...")
  if (show_progress && !is.null(pb_overall)) pb_overall$tick(0)
  statistical_results_list <- list()
  if (length(tibble_list) == 0) { verbose_msg_internal("Input `tibble_list` empty. No tests run."); return(statistical_results_list) }
  for (name_of_cell_type_tibble in names(tibble_list)) {
    current_tibble <- tibble_list[[name_of_cell_type_tibble]]
    result_key_name <- make.names(paste0(name_of_cell_type_tibble, ".", user_suffix, ".Average"))
    verbose_msg_internal(paste("  - Processing:", name_of_cell_type_tibble, "-> applying test_wilcox_fn"))
    single_result_df <- tryCatch({
      if(is.null(current_tibble) || nrow(current_tibble) == 0) {warning(paste("Tibble '", name_of_cell_type_tibble, "' NULL/empty. Skipping."), call.=F); NULL}
      else test_wilcox_fn(current_tibble)
    }, error = function(e) {warning(paste("Error in `test_wilcox_fn` for '", name_of_cell_type_tibble, "': ", e$message,". NULL returned."), call.=F); NULL})
    if (!is.null(single_result_df) && !is.data.frame(single_result_df)){warning(paste0("`test_wilcox_fn` for '", name_of_cell_type_tibble,"' non-DF. NULL used."), call.=F); single_result_df <- NULL}
    statistical_results_list[[result_key_name]] <- single_result_df
    if (show_progress && !is.null(pb_test_application)) pb_test_application$tick()
  }
  if (show_progress && !is.null(pb_overall)) pb_overall$tick()
  if (verbose) {
    verbose_msg_internal("\nSummary of statistical processing:"); verbose_msg_internal(sprintf("  Tibbles processed: %d", length(tibble_list)))
    successful_results <- Filter(Negate(is.null), statistical_results_list); verbose_msg_internal(sprintf("  Non-NULL DFs generated: %d", length(successful_results)))
  }
  return(statistical_results_list)
}

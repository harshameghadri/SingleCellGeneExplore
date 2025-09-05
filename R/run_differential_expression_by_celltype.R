#' Run Differential Expression Analysis Per Cell Type Between Two Conditions
#'
#' @description
#' Performs differential gene expression (DE) analysis using `Seurat::FindMarkers`
#' for each unique cell type present in a Seurat object. The comparison is made
#' between two specified conditions (e.g., treatment vs. control, disease vs. healthy)
#' found within a designated metadata column.
#'
#' @details
#' The function iterates through each unique identifier found in the `celltype_col`
#' of the Seurat object's metadata. For each cell type:
#' \enumerate{
#'   \item It subsets the Seurat object to include only cells of that specific type.
#'   \item It checks if there are cells belonging to both `condition1` and `condition2`
#'         (as defined in the `condition_col`) within the subset. If not, it skips
#'         that cell type with a warning.
#'   \item It sets the identity (`Idents`) of the subset to the `condition_col`.
#'   \item It runs `Seurat::FindMarkers` comparing `ident.1 = condition1` against
#'         `ident.2 = condition2` using the specified statistical test and thresholds.
#'         Positive `avg_log2FC` values indicate higher expression in `condition1`.
#'   \item The results (including gene names) for the cell type are saved as a
#'         tab-delimited text file in the `output_dir`.
#'   \item Optionally, if `save_as_excel = TRUE` and the `openxlsx` package is
#'         installed, all results are compiled into a single Excel workbook. This
#'         workbook includes a summary sheet and separate sheets for each cell type's
#'         detailed DE results, with conditional formatting applied to the logFC column.
#'         If `openxlsx` is not installed, a message is printed, and only individual
#'         `.txt` files are saved.
#' }
#' The function can optionally return the results as a named list.
#'
#' @param seurat_obj A Seurat object containing single-cell data and metadata.
#' @param celltype_col Character string. The name of the metadata column in
#'   `seurat_obj@meta.data` that contains the cell type annotations.
#'   Default is `"cell.type.ident"`.
#' @param condition_col Character string. The name of the metadata column in
#'   `seurat_obj@meta.data` that contains the condition labels (e.g., "disease",
#'   "treatment"). Default is `"disease.ident"`.
#' @param condition1 Character string. The label for the first condition (ident.1
#'   in `FindMarkers`, the numerator for fold change). Default is `"PAH"`.
#' @param condition2 Character string. The label for the second condition (ident.2
#'   in `FindMarkers`, the denominator for fold change). Default is `"CTRL"`.
#' @param output_dir Character string. The path to the directory where output files
#'   (`.txt` files per cell type and optionally a `.xlsx` file) will be saved.
#'   The directory will be created if it doesn't exist. Default is `"./DE_results"`.
#' @param output_filename Character string. The base name for the combined Excel
#'   output file (without the .xlsx extension). Default is `"DE_results"`. Ignored
#'   if `save_as_excel = FALSE`.
#' @param test_use Character string. The statistical test to use in `FindMarkers`.
#'   See `?Seurat::FindMarkers` for options (e.g., "wilcox", "MAST", "DESeq2").
#'   Default is `"wilcox"` (Wilcoxon rank sum test).
#' @param min_pct Numeric. Only test genes that are detected in a minimum fraction
#'   of cells in either of the two populations being compared. Passed to `FindMarkers`.
#'   Default is `0.1`.
#' @param logfc_threshold Numeric. Limit testing to genes which show, on average,
#'   at least X-fold difference (log-scale) between the two conditions. Passed to
#'   `FindMarkers`. Default is `0.25`. Setting to 0 tests all genes meeting `min_pct`.
#' @param only_pos Logical. If `TRUE`, only return genes with positive `avg_log2FC`
#'   (i.e., more highly expressed in `condition1`). Passed to `FindMarkers`.
#'   Default is `FALSE`.
#' @param return_results Logical. If `TRUE`, the function returns a named list
#'   where each element contains the differential expression results (as a data frame)
#'   for a specific cell type. If `FALSE` (default), the function returns `NULL`
#'   invisibly, primarily relying on side effects (saving files).
#' @param save_as_excel Logical. If `TRUE` (default), attempts to save a combined
#'   Excel file containing all results using the `openxlsx` package (which must be
#'   installed). If `FALSE`, only individual `.txt` files are saved.
#' @param verbose Logical. If `TRUE` (default), print progress messages to the
#'   console during execution.
#'
#' @return If `return_results` is `TRUE`, a named list where names are cell types
#'   and values are data frames containing the DE results from `FindMarkers` for that
#'   cell type. Each data frame includes a 'gene' column.
#'   If `return_results` is `FALSE`, returns `invisible(NULL)`.
#'   The function primarily acts via side effects: creating text files for each
#'   cell type's DE results in `output_dir`, and optionally creating a combined
#'   Excel file.
#'
#' @importFrom Seurat Idents<- FindMarkers subset
#' @importFrom utils write.table setdiff head
#' @importFrom methods exists # Not explicitly used, but good practice for pkg checks
#'
#' @seealso \code{\link[Seurat]{FindMarkers}}, \code{\link[openxlsx]{createWorkbook}} (if saving Excel)
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # This example requires the Seurat package.
#' # install.packages("Seurat")
#' # install.packages("openxlsx") # Optional, for Excel output
#' # library(Seurat)
#' # library(dplyr) # For potential data manipulation before/after
#'
#' # --- Create a dummy Seurat object ---
#' counts <- matrix(rpois(20000, lambda = 2), nrow = 200, ncol = 100)
#' rownames(counts) <- paste0("Gene", 1:200)
#' colnames(counts) <- paste0("Cell", 1:100)
#' seu <- CreateSeuratObject(counts = counts)
#'
#' # Add dummy metadata
#' set.seed(123)
#' seu$cell.type.ident <- sample(c("Tcell", "Bcell", "Myeloid"), ncol(seu), replace = TRUE)
#' seu$disease.ident <- sample(c("PAH", "CTRL"), ncol(seu), replace = TRUE)
#' seu <- NormalizeData(seu, verbose = FALSE)
#' # --- End Dummy Data ---
#'
#' # Define output directory (use tempdir for example)
#' out_dir <- file.path(tempdir(), "DE_Example")
#'
#' # Run the DE analysis
#' de_results_list <- run_differential_expression_by_celltype(
#'  seurat_obj = seu,
#'  celltype_col = "cell.type.ident",
#'  condition_col = "disease.ident",
#'  condition1 = "PAH",
#'  condition2 = "CTRL",
#'  output_dir = out_dir,
#'  output_filename = "PAH_vs_CTRL_DE",
#'  test_use = "wilcox",
#'  min_pct = 0.1,
#'  logfc_threshold = 0.1, # Lower threshold for more results in dummy data
#'  return_results = TRUE, # Get the results back as a list
#'  save_as_excel = TRUE, # Try to save Excel file
#'  verbose = TRUE
#' )
#'
#' # Explore results
#' print(paste("Results saved in:", out_dir))
#' list.files(out_dir)
#'
#' # If return_results = TRUE, inspect the list
#' if (!is.null(de_results_list)) {
#'  names(de_results_list) # See results for which cell types
#'  # head(de_results_list$Tcell) # View head of Tcell results
#' }
#'
#' # Clean up temporary directory
#' # unlink(out_dir, recursive = TRUE)
#' }
#'
run_differential_expression_by_celltype <- function(
    seurat_obj,
    celltype_col = "cell.type.ident",  # Column containing cell type info
    condition_col = "disease.ident",   # Column containing condition info
    condition1 = "PAH",                # First condition (numerator)
    condition2 = "CTRL",               # Second condition (denominator)
    output_dir = "./DE_results",       # Directory for output files
    output_filename = "DE_results",    # Base name for output files
    test_use = "wilcox",               # Statistical test to use
    min_pct = 0.1,                     # Min percentage of cells expressing gene
    logfc_threshold = 0.25,            # Log fold-change threshold
    only_pos = FALSE,                  # Whether to only return positive markers
    return_results = FALSE,            # Whether to return results as a list
    save_as_excel = TRUE,              # Whether to save results as Excel file
    verbose = TRUE                     # Print progress messages
) {
  
  # --- Input Validation ---
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat package is required but not installed")
  }
  if (!inherits(seurat_obj, "Seurat")) {
    stop("'seurat_obj' must be a Seurat object.")
  }
  if (!celltype_col %in% names(seurat_obj@meta.data)) {
    stop("Cell type column '", celltype_col, "' not found in Seurat object metadata.")
  }
  if (!condition_col %in% names(seurat_obj@meta.data)) {
    stop("Condition column '", condition_col, "' not found in Seurat object metadata.")
  }
  # Check if conditions exist in the specified column
  available_conditions <- unique(seurat_obj@meta.data[[condition_col]])
  if (!condition1 %in% available_conditions) {
    stop("condition1 '", condition1, "' not found in metadata column '", condition_col, "'. Available: ", paste(available_conditions, collapse=", "))
  }
  if (!condition2 %in% available_conditions) {
    stop("condition2 '", condition2, "' not found in metadata column '", condition_col, "'. Available: ", paste(available_conditions, collapse=", "))
  }
  
  
  # Check required packages conditionally
  excel_pkg_available <- FALSE
  if (save_as_excel) {
    if (!requireNamespace("openxlsx", quietly = TRUE)) {
      if(verbose) message("NOTE: 'openxlsx' package needed for Excel output (save_as_excel=TRUE) but not installed. Will save only .txt files.")
      save_as_excel <- FALSE # Fallback to saving only txt
    } else {
      excel_pkg_available <- TRUE
    }
  }
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    if(verbose) message("Creating output directory: ", output_dir)
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Get unique cell types
  cell_types <- unique(as.character(seurat_obj@meta.data[[celltype_col]]))
  if(length(cell_types) == 0) {
    stop("No unique cell types found in column '", celltype_col, "'.")
  }
  
  
  if(verbose) {
    message(sprintf("Found %d unique cell types in '%s': %s",
                    length(cell_types), celltype_col, paste(cell_types, collapse=", ")))
    message(sprintf("Starting DE analysis comparing '%s' vs '%s' using column '%s'.",
                    condition1, condition2, condition_col))
  }
  
  # Initialize results list
  all_results <- list()
  
  # --- Process each cell type ---
  for(cell_type in cell_types) {
    current_cell_type_label <- gsub("[^a-zA-Z0-9_.-]", "_", cell_type) # Sanitize for messages/files
    
    if(verbose) {
      message(paste("Processing cell type:", cell_type))
    }
    
    # Subset to specific cell type using cell identifiers for efficiency
    cells_to_subset <- rownames(seurat_obj@meta.data[seurat_obj@meta.data[[celltype_col]] == cell_type, , drop = FALSE])
    
    if(length(cells_to_subset) == 0) {
      warning(paste("No cells found for cell type", cell_type, ". Skipping."))
      next
    }
    cell_subset <- subset(seurat_obj, cells = cells_to_subset)
    
    # Check if both conditions exist in this cell type subset
    subset_meta <- cell_subset@meta.data
    cond1_cells_idx <- subset_meta[[condition_col]] == condition1
    cond2_cells_idx <- subset_meta[[condition_col]] == condition2
    cond1_count <- sum(cond1_cells_idx)
    cond2_count <- sum(cond2_cells_idx)
    
    # Check for sufficient cells in both conditions
    min_cells_per_group <- 3 # Seurat default, maybe make this a parameter?
    if(cond1_count < min_cells_per_group || cond2_count < min_cells_per_group) {
      warning(paste("Cell type", cell_type, "has fewer than", min_cells_per_group,
                    "cells in one or both conditions (", condition1, ":", cond1_count, ", ",
                    condition2, ":", cond2_count, "). Skipping DE analysis for this type."))
      next
    }
    
    if(verbose) {
      message(paste("  Found", cond1_count, condition1, "cells and", cond2_count, condition2, "cells"))
    }
    
    # Set identity to condition (ensure it's a factor for FindMarkers)
    Seurat::Idents(cell_subset) <- factor(subset_meta[[condition_col]])
    
    # Run differential expression
    markers <- NULL # Initialize markers variable
    tryCatch({
      markers <- Seurat::FindMarkers(
        object = cell_subset,
        ident.1 = condition1,
        ident.2 = condition2,
        test.use = test_use,
        min.pct = min_pct,
        only.pos = only_pos,
        logfc.threshold = logfc_threshold,
        verbose = FALSE # Reduce Seurat's own verbosity within the loop
      )
    }, error = function(e) {
      warning(paste("Error running FindMarkers for cell type", cell_type, ":", e$message))
      # markers remains NULL
    })
    
    
    # Check if markers were successfully computed and not empty
    if (!is.null(markers) && nrow(markers) > 0) {
      # Add gene names as column and move it to the front
      markers$gene <- rownames(markers)
      markers <- markers[, c("gene", setdiff(colnames(markers), "gene")), drop = FALSE]
      
      # Store in results list
      all_results[[cell_type]] <- markers
      
      # --- Save individual result as .txt ---
      # Sanitize cell type name for filename
      safe_cell_type_name <- gsub("[^a-zA-Z0-9_.-]", "_", cell_type)
      output_txt_name <- paste0(safe_cell_type_name, "_", condition1, "_vs_", condition2, ".txt")
      file_path_txt <- file.path(output_dir, output_txt_name)
      
      tryCatch({
        write.table(
          markers,
          file = file_path_txt,
          sep = "\t",
          row.names = FALSE,
          quote = FALSE
        )
        if(verbose) {
          message(paste("  Saved results to", file_path_txt))
        }
      }, error = function(e) {
        warning(paste("Failed to save .txt results for", cell_type, "to", file_path_txt, ":", e$message))
      })
      
    } else {
      if(verbose) {
        message(paste("  No significant markers found or error occurred for cell type:", cell_type))
      }
      # Ensure an entry exists if returning list, maybe NULL or empty df? Let's use NULL.
      all_results[[cell_type]] <- NULL
    }
  } # End loop through cell types
  
  # --- Save combined results as Excel file if requested ---
  valid_results <- all_results[!sapply(all_results, is.null)] # Filter out NULL results
  if(save_as_excel && excel_pkg_available && length(valid_results) > 0) {
    excel_file <- file.path(output_dir, paste0(output_filename, ".xlsx"))
    if(verbose) message("Preparing combined Excel file...")
    
    tryCatch({
      wb <- openxlsx::createWorkbook()
      
      # Add a summary sheet first
      openxlsx::addWorksheet(wb, "Summary")
      summary_list <- lapply(names(valid_results), function(ct) {
        res_df <- valid_results[[ct]]
        up_genes <- res_df$gene[res_df$avg_log2FC > logfc_threshold] # Use threshold consistent? Or just >0? Let's use >0.
        down_genes <- res_df$gene[res_df$avg_log2FC < -logfc_threshold] # Use threshold consistent? Or just <0? Let's use <0.
        # Order by p-value or FC before taking top 5? Assume FindMarkers order is okay.
        data.frame(
          Cell_Type = ct,
          Num_DEGs = nrow(res_df),
          Num_Upregulated = sum(res_df$avg_log2FC > 0),
          Num_Downregulated = sum(res_df$avg_log2FC < 0),
          Top5_Upregulated = paste(utils::head(res_df$gene[res_df$avg_log2FC > 0], 5), collapse = ", "),
          Top5_Downregulated = paste(utils::head(res_df$gene[res_df$avg_log2FC < 0], 5), collapse = ", ")
        )
      })
      summary_df <- do.call(rbind, summary_list)
      openxlsx::writeData(wb, "Summary", summary_df)
      openxlsx::setColWidths(wb, "Summary", cols = 1:ncol(summary_df), widths = "auto")
      
      # Add detailed sheets for each cell type
      for(cell_type in names(valid_results)) {
        res_df <- valid_results[[cell_type]]
        # Sanitize sheet name
        sheet_name <- gsub("[^a-zA-Z0-9]", "_", cell_type)
        if(nchar(sheet_name) > 31) { sheet_name <- substr(sheet_name, 1, 31) }
        # Handle potential duplicate sheet names after sanitization (append number)
        original_sheet_name <- sheet_name
        sheet_counter <- 1
        while (sheet_name %in% names(wb)) {
          sheet_name <- paste0(substr(original_sheet_name, 1, 31 - nchar(as.character(sheet_counter)) - 1), "_", sheet_counter)
          sheet_counter <- sheet_counter + 1
          if(nchar(sheet_name) > 31) sheet_name <- substr(sheet_name, 1, 31) # Failsafe
        }
        
        
        openxlsx::addWorksheet(wb, sheet_name)
        openxlsx::writeData(wb, sheet_name, res_df)
        openxlsx::freezePane(wb, sheet_name, firstRow = TRUE) # Freeze header row
        openxlsx::setColWidths(wb, sheet_name, cols = 1:ncol(res_df), widths = "auto")
        
        
        # Add conditional formatting (check if columns exist first)
        if("avg_log2FC" %in% colnames(res_df) && nrow(res_df)>0) {
          style_positive <- openxlsx::createStyle(fontColour = "#006100", bgFill = "#C6EFCE") # Dark Green font, Light Green fill
          style_negative <- openxlsx::createStyle(fontColour = "#9C0006", bgFill = "#FFC7CE") # Dark Red font, Light Red fill
          logfc_col_idx <- which(colnames(res_df) == "avg_log2FC")
          
          openxlsx::conditionalFormatting(wb, sheet_name, cols = logfc_col_idx,
                                          rows = 2:(nrow(res_df) + 1), rule = ">0", style = style_positive)
          openxlsx::conditionalFormatting(wb, sheet_name, cols = logfc_col_idx,
                                          rows = 2:(nrow(res_df) + 1), rule = "<0", style = style_negative)
        }
      }
      
      openxlsx::saveWorkbook(wb, excel_file, overwrite = TRUE)
      
      if(verbose) {
        message(paste("Saved combined results to Excel file:", excel_file))
      }
    }, error = function(e) {
      warning(paste("Failed to create or save Excel workbook:", e$message))
    })
    
  } else if (save_as_excel && !excel_pkg_available && verbose) {
    # Message already printed earlier about package missing
  } else if (save_as_excel && length(valid_results) == 0 && verbose) {
    message("No valid DE results found for any cell type, skipping Excel file creation.")
  }
  
  
  # --- Optional return of results ---
  if(return_results) {
    # Return only the valid results
    return(valid_results)
  } else {
    return(invisible(NULL))
  }
}

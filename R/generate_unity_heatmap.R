#' Create Unity or Z-score Normalized Heatmap by Subject, Disease, and Cell Type
#'
#' @description
#' This function generates a heatmap visualizing gene expression patterns across
#' specified cell types, conditions (diseases), and subjects. It calculates
#' pseudo-bulk expression profiles by averaging gene expression within each
#' unique combination of cell type, disease status, and subject ID, applying
#' a minimum cell count filter per combination. The resulting matrix can be
#' normalized using either Unity (0-1 scaling) or Z-score normalization.
#'
#' @param seurat_obj A Seurat object containing the single-cell RNA-seq data.
#'   The specified `assay` should contain the expression data (e.g., counts,
#'   normalized data).
#' @param celltypes_to_plot A character vector specifying the cell types
#'   (matching values in `celltype_column`) to include in the heatmap.
#' @param genes_for_heatmap A character vector of gene names to display
#'   in the heatmap rows. Genes not found in the Seurat object will be skipped.
#' @param Cell_Colors A named vector or list providing color codes for each
#'   cell type specified in `celltypes_to_plot`. Names should match the cell types.
#' @param disease_colors A named vector or list providing color codes for each
#'   disease status level (e.g., `c(CTRL = "blue", PAH = "red")`). Names should
#'   match the levels in `disease_column`.
#' @param subject_colors A named vector or list providing color codes for each
#'   subject ID (values in `subject_column`). Names should match the subject IDs.
#' @param celltype_column A character string specifying the column name in
#'   `seurat_obj@meta.data` that contains cell type annotations.
#' @param disease_column A character string specifying the column name in
#'   `seurat_obj@meta.data` that contains disease status annotations (e.g., "CTRL", "PAH").
#' @param subject_column A character string specifying the column name in
#'   `seurat_obj@meta.data` that contains subject IDs (e.g., "orig.ident", "patient_id").
#'   Defaults to "orig.ident".
#' @param assay A character string specifying the name of the assay in the Seurat
#'   object containing the expression data to use. Defaults to "RNA".
#' @param min.cells.per.subject An integer specifying the minimum number of cells
#'   required for a specific cell type within a specific subject under a specific
#'   disease condition to be included in the pseudo-bulk calculation. Defaults to 3.
#'   Combinations with fewer cells are excluded.
#' @param cores An integer specifying the number of cores to use for parallel
#'   processing (`mclapply`). Defaults to 10. Note: `mclapply` only works on
#'   Linux/macOS, not Windows. Consider using `BiocParallel` for cross-platform compatibility.
#' @param pdf_filename A character string specifying the filename for the output
#'   PDF heatmap. Defaults to "heatmap.pdf".
#' @param png_filename A character string specifying the filename for the output
#'   PNG heatmap. Defaults to "heatmap.png".
#' @param normalization_method A character string specifying the normalization
#'   method to apply to the heatmap data. Can be "unity" (0-1 scaling) or "zscore".
#'   Defaults to "unity".
#' @param zscore_colors A color ramp function or a vector of colors to use for
#'   the heatmap when `normalization_method` is "zscore". If a vector, a color
#'   interpolation function will be created. Defaults to a divergent color scale
#'   (`circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))`).
#' @param verbose A logical value indicating whether to print progress messages
#'   and warnings (e.g., skipped genes, filtered combinations). Defaults to TRUE.
#'
#' @return Invisibly returns the `Heatmap` object created by the
#'   `ComplexHeatmap` package. This allows further modification or drawing
#'   of the heatmap object if needed.
#'
#' @details
#' The function performs the following steps:
#' 1. Validates inputs, including checking for the existence of specified columns
#'    and the assay, and filters genes present in the Seurat object.
#' 2. Identifies unique combinations of cell type, disease status, and subject ID
#'    based on the specified metadata columns.
#' 3. Filters these combinations based on `min.cells.per.subject`. Combinations
#'    with fewer cells are excluded from further analysis and the heatmap.
#' 4. For each valid combination, calculates the average expression (pseudo-bulk)
#'    of the specified `genes_for_heatmap` using cells belonging to that combination.
#'    Uses `Matrix::rowMeans` for efficiency with sparse matrices.
#' 5. Collates the pseudo-bulk profiles into a matrix.
#' 6. Normalizes the pseudo-bulk matrix row-wise (per gene) using the specified
#'    `normalization_method`. "unity" scales values between 0 and 1. "zscore"
#'    calculates the standard score for each value relative to the gene's mean
#'    and standard deviation across all samples. Genes with constant expression
#'    (zero range or standard deviation) are normalized to 0.
#' 7. Creates column annotations for the heatmap based on cell type, disease,
#'    and subject, using the provided color mappings.
#' 8. Orders the heatmap columns based on cell type, then disease, then subject ID
#'    according to the levels found in the filtered data and the order of
#'    `celltypes_to_plot`.
#' 9. Generates the heatmap using `ComplexHeatmap::Heatmap`, disabling column clustering
#'    and splitting columns by cell type. The color scale used depends on the
#'    `normalization_method` and the `zscore_colors` parameter if applicable.
#' 10. Saves the heatmap to both PDF and PNG files.
#'
#' **Notes:**
#' * The function name `generate_Unity_heatmap` uses dots and reflects only one
#'   normalization method. Consider renaming to `generatePseudoBulkHeatmap` or
#'   `generate_pseudo_bulk_heatmap` for better clarity and style consistency.
#' * Parallel processing uses `parallel::mclapply`, which is not available on Windows.
#'   For cross-platform compatibility in a package, consider alternatives like
#'   `BiocParallel::bplapply` or adding an option for sequential processing.
#'
#' @importFrom dplyr filter pull arrange sym
#' @importFrom parallel mclapply
#' @importFrom Matrix rowMeans
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation draw
#' @importFrom viridis inferno
#' @importFrom methods as
#' @importFrom grDevices pdf dev.off png colorRampPalette
#' @importFrom stats sd
#' @importFrom circlize colorRamp2
#' @import SeuratObject
#' @import Seurat
#'
#' @export
generatePseudoBulkHeatmap <- function(seurat_obj,
                                   celltypes_to_plot,
                                   genes_for_heatmap,
                                   Cell_Colors,
                                   disease_colors,
                                   subject_colors,
                                   celltype_column,
                                   disease_column,
                                   subject_column = "orig.ident",
                                   assay = "RNA",
                                   min.cells.per.subject = 3,
                                   cores = 10,
                                   pdf_filename = "heatmap.pdf",
                                   png_filename = "heatmap.png",
                                   normalization_method = "unity",
                                   zscore_colors = circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
                                   verbose = TRUE) {

  # --- Input validation ---
  if(!all(celltypes_to_plot %in% unique(seurat_obj@meta.data[[celltype_column]]))) {
    missing_types <- setdiff(celltypes_to_plot, unique(seurat_obj@meta.data[[celltype_column]]))
    stop("Some cell types not found in data: ", paste(missing_types, collapse = ", "))
  }

  if (!subject_column %in% colnames(seurat_obj@meta.data)) {
    stop("Subject ID column '", subject_column, "' not found in metadata.")
  }

  if (!disease_column %in% colnames(seurat_obj@meta.data)) {
    stop("Disease column '", disease_column, "' not found in metadata.")
  }

  if (!celltype_column %in% colnames(seurat_obj@meta.data)) {
    stop("Celltype column '", celltype_column, "' not found in metadata.")
  }

  if (!assay %in% Assays(seurat_obj)) {
    stop("Assay '", assay, "' not found in Seurat object.")
  }

  # --- Validate and filter genes ---
  # Use the specified assay for gene availability check
  available_genes <- rownames(GetAssayData(seurat_obj, assay = assay))
  valid_genes <- intersect(genes_for_heatmap, available_genes)
  missing_genes <- setdiff(genes_for_heatmap, available_genes)

  if(length(valid_genes) == 0) {
    stop("None of the specified genes were found in the dataset in assay '", assay, "'!")
  }

  if(length(missing_genes) > 0 && verbose) {
    message("The following genes were not found in the dataset in assay '", assay, "' and will be skipped: \n",
            paste(missing_genes, collapse = ", "), "\n")
  }

  # Update genes_for_heatmap to only include valid genes
  genes_for_heatmap <- valid_genes

  if(verbose) {
    message("Proceeding with ", length(genes_for_heatmap), " genes using assay '", assay, "'")
    message("Using min.cells.per.subject = ", min.cells.per.subject)
    message("Normalization method: ", normalization_method)
  }

  # --- Setup cell barcodes and metadata subset ---
  seurat_obj$cellBarcode <- colnames(seurat_obj)
  # Ensure cell types are factor for consistent ordering if needed later
  seurat_obj@meta.data[[celltype_column]] <- factor(seurat_obj@meta.data[[celltype_column]])
  cellTypes <- levels(seurat_obj@meta.data[[celltype_column]])
  cellTypes <- cellTypes[cellTypes %in% celltypes_to_plot] # Filter levels based on input

  # Select necessary metadata columns, using the specified subject_column
  meta.data.sub <- seurat_obj@meta.data[, c(disease_column, celltype_column, subject_column, "cellBarcode")]

  # --- Helper functions for normalization ---
  # Helper function to normalize data (0-1 scale) - Unity Normalization
  myUnityNormalize <- function(x) {
    # Handle cases with all NA, all zero, or constant values
    if(all(is.na(x))) return(x)
    min_x <- min(x, na.rm = TRUE)
    max_x <- max(x, na.rm = TRUE)
    range_x = max_x - min_x
    if (is.na(range_x) || range_x == 0) {
      # Return 0 if range is zero or NA (constant value)
      return(rep(0, length(x)))
    }
    (x - min_x) / range_x
  }

  # Helper function for Z-score Normalization
  myZscoreNormalize <- function(x) {
    # Handle cases with all NA or constant values (sd = 0)
    if(all(is.na(x))) return(x)
    mean_x <- mean(x, na.rm = TRUE)
    sd_x <- sd(x, na.rm = TRUE)
    if (is.na(sd_x) || sd_x == 0) {
      # Return 0 if standard deviation is zero or NA (constant value)
      return(rep(0, length(x)))
    }
    (x - mean_x) / sd_x
  }

  # --- Generate cell type, disease, and subject combinations with filtering ---
  get.CT.DS.subj.vector <- function(currentCellType){
    # Filter metadata for the current cell type being processed
    tmp.meta.data <- meta.data.sub %>% dplyr::filter(!!dplyr::sym(celltype_column) == currentCellType)
    # Find unique diseases and subjects *within this cell type*
    diseases <- unique(tmp.meta.data[[disease_column]])
    subjects <- unique(tmp.meta.data[[subject_column]]) # Use subject_column
    tmp.CT.DS.subj.vector <- vector() # Initialize empty vector for results
    filtered_count <- 0 # Counter for filtered combinations

    # Iterate through diseases found for this cell type
    for(j in 1:length(diseases)){
      # Iterate through subjects found for this cell type
      for(k in 1:length(subjects)){
        # Get cell barcodes for the specific combination
        temp.cells <- tmp.meta.data %>%
          dplyr::filter(!!dplyr::sym(disease_column) == diseases[j] & !!dplyr::sym(subject_column) == subjects[k]) %>% # Use subject_column
          dplyr::pull(cellBarcode)

        # Only include if cell count meets minimum threshold
        if (length(temp.cells) >= min.cells.per.subject) {
          # Add the valid combination string to the result vector
          tmp.CT.DS.subj.vector <- c(tmp.CT.DS.subj.vector,
                                     paste(currentCellType, diseases[j], subjects[k], sep = "__"))
        } else if (length(temp.cells) > 0) {
          # If cells exist but are below threshold, count and optionally report
          filtered_count <- filtered_count + 1
          if(verbose) {
            # Corrected message to use min.cells.per.subject
            cat("Filtered out: ", currentCellType, "/", diseases[j], "/", subjects[k],
                " - only has ", length(temp.cells), " cells (minimum ",
                min.cells.per.subject, " required)\n", sep = "")
          }
        }
        # If length(temp.cells) == 0, the combination doesn't exist, do nothing
      }
    }

    # Optional reporting after processing a cell type
    if(verbose) {
      if(filtered_count > 0) {
        # Corrected message to use min.cells.per.subject
        cat("Completed for ", currentCellType, ". Filtered out ", filtered_count,
            " combinations with fewer than ", min.cells.per.subject, " cells.\n", sep = "")
      } else {
        cat("Completed for ", currentCellType, ".\n", sep = "")
      }
    }
    return(tmp.CT.DS.subj.vector) # Return vector of valid combinations for this cell type
  }

  # Apply the function in parallel across the specified cell types
  # Note: mclapply does not work on Windows
  if (verbose) message("Generating valid combinations using ", cores, " cores...")
  celltype_disease_subject.list <- parallel::mclapply(cellTypes, get.CT.DS.subj.vector,
                                                      mc.cores = cores)
  # Combine results from all cell types into a single vector
  celltype_disease_subject <- unlist(celltype_disease_subject.list)

  # Stop if no combinations remain after filtering
  if(length(celltype_disease_subject) == 0) {
    stop("No cell type/disease/subject combinations remain after filtering. Try reducing min.cells.per.subject.")
  }

  # --- Calculate average expression for valid combinations ---
  get.SubjectdiseaseCellTypeAvg <- function(combo_string){
    # Parse the combination string
    split_result <- strsplit(as.character(combo_string), "__")[[1]]
    temp.cell.type <- split_result[1]
    temp.disease <- split_result[2]
    temp.subject <- split_result[3]

    # Re-filter metadata to get cell barcodes for this specific combination
    temp.meta.data <- seurat_obj@meta.data[, c(disease_column, celltype_column,
                                               subject_column, "cellBarcode")] # Use subject_column
    temp.cells <- temp.meta.data %>%
      dplyr::filter(!!dplyr::sym(celltype_column) == temp.cell.type &
                      !!dplyr::sym(disease_column) == temp.disease &
                      !!dplyr::sym(subject_column) == temp.subject) %>% # Use subject_column
      dplyr::pull(cellBarcode)

    # Check cell count again (should match filter, but acts as safety check)
    if (length(temp.cells) < min.cells.per.subject && length(temp.cells) > 0) {
      if(verbose) {
        # Corrected message to use min.cells.per.subject
        cat("Warning: Unexpectedly found fewer than ", min.cells.per.subject,
            " cells for ", temp.cell.type, " in subject ", temp.subject,
            " for ", temp.disease, " - using available ", length(temp.cells), " cells\n", sep="")
      }
      # Current logic proceeds if length(temp.cells) > 0
    } else if (length(temp.cells) == 0) {
      # If somehow no cells are found (shouldn't happen due to filter), return NAs
      tmp.df <- as.data.frame(matrix(NA, nrow = length(genes_for_heatmap), ncol = 1))
      rownames(tmp.df) <- genes_for_heatmap
      colnames(tmp.df) <- combo_string
      if(verbose) {
        cat("Warning: No cells found for combination ", combo_string, " during average calculation.\n", sep = "")
      }
      return(tmp.df)
    }


    # Calculate pseudo-bulk profile using the specified assay
    if (length(temp.cells) > 1) {
      # Use Matrix::rowMeans for sparse matrix efficiency
      tmp.df <- as.data.frame(Matrix::rowMeans(GetAssayData(seurat_obj, assay = assay)[genes_for_heatmap, temp.cells, drop = FALSE])) # Specify assay and subset genes
    } else { # length(temp.cells) == 1
      # If only one cell, use its expression directly
      tmp.df <- as.data.frame(GetAssayData(seurat_obj, assay = assay)[genes_for_heatmap, temp.cells, drop = FALSE]) # Specify assay and subset genes
      if(verbose) {
        cat("Note: Subject ", temp.subject, " only has 1 ", temp.cell.type,
            " cell for ", temp.disease, ", using singlet expression.\n", sep = "")
      }
    }

    # Set column name of the resulting pseudo-bulk vector
    colnames(tmp.df) <- combo_string
    return(tmp.df)
  }

  # Apply pseudo-bulk calculation in parallel
  # Note: mclapply does not work on Windows
  if (verbose) message("Calculating pseudo-bulk profiles using ", cores, " cores...")
  collapsed.mtx.list <- parallel::mclapply(celltype_disease_subject,
                                           get.SubjectdiseaseCellTypeAvg,
                                           mc.cores = cores)

  # Combine the list of single-column data frames into one matrix
  collapsed.SubjectdiseaseCellTypeAvg.mtx <- do.call(cbind, collapsed.mtx.list)

  if(verbose) {
    total_combinations <- length(celltype_disease_subject)
    message("Using ", total_combinations, " cell type/disease/subject combinations for heatmap")
    message("Dimensions of pseudo-bulk matrix: ",
            paste(dim(collapsed.SubjectdiseaseCellTypeAvg.mtx), collapse = " x "))
  }

  # --- Prepare for Heatmap ---

  # Create metadata frame for heatmap columns
  heatmap_metadata <- data.frame(
    cell.ident = colnames(collapsed.SubjectdiseaseCellTypeAvg.mtx),
    cell_type = sapply(strsplit(colnames(collapsed.SubjectdiseaseCellTypeAvg.mtx), "__"), `[`, 1),
    disease = sapply(strsplit(colnames(collapsed.SubjectdiseaseCellTypeAvg.mtx), "__"), `[`, 2),
    subject = sapply(strsplit(colnames(collapsed.SubjectdiseaseCellTypeAvg.mtx), "__"), `[`, 3),
    stringsAsFactors = FALSE
  )
  # Assign proper column names based on input parameters
  colnames(heatmap_metadata) <- c("cell.ident", celltype_column, disease_column, subject_column)


  # Convert relevant columns to factors with specific levels for ordering and coloring
  # Use the order of celltypes_to_plot for cell type factor levels
  heatmap_metadata[[celltype_column]] <- factor(heatmap_metadata[[celltype_column]],
                                                levels = celltypes_to_plot)

  # Infer disease levels from the data present in the filtered metadata
  disease_levels_present <- unique(heatmap_metadata[[disease_column]])
  heatmap_metadata[[disease_column]] <- factor(heatmap_metadata[[disease_column]],
                                               levels = disease_levels_present) # Use levels found in data

  # Infer subject levels from the data present in the filtered metadata
  subject_levels_present <- unique(heatmap_metadata[[subject_column]])
  heatmap_metadata[[subject_column]] <- factor(heatmap_metadata[[subject_column]],
                                               levels = subject_levels_present) # Use levels found in data


  # Order the metadata frame to determine column order in heatmap
  heatmap_metadata_ordered <- heatmap_metadata %>%
    dplyr::arrange(!!dplyr::sym(celltype_column), !!dplyr::sym(disease_column), !!dplyr::sym(subject_column)) # Arrange by subject_column

  # Extract ordered vectors for annotation and data ordering
  cell_order <- heatmap_metadata_ordered$cell.ident
  disease_order <- heatmap_metadata_ordered[[disease_column]]
  celltype_order <- heatmap_metadata_ordered[[celltype_column]]
  subject_order <- heatmap_metadata_ordered[[subject_column]] # Use subject_column

  # Create heatmap data matrix, ordered by columns
  # Ensure only valid genes are selected and columns are ordered
  heatmap_df <- as.matrix(collapsed.SubjectdiseaseCellTypeAvg.mtx[genes_for_heatmap, cell_order, drop = FALSE])

  # Apply selected normalization row-wise (per gene)
  if (verbose) message("Applying ", normalization_method, " normalization...")
  if (normalization_method == "unity") {
    heatmap_df_normalized <- t(apply(heatmap_df, MARGIN = 1, FUN = myUnityNormalize))
    norm_legend_title <- "Unity Norm. Expr."
    heatmap_color_scale <- viridis::inferno(256) # Default Unity color scale
  } else if (normalization_method == "zscore") {
    heatmap_df_normalized <- t(apply(heatmap_df, MARGIN = 1, FUN = myZscoreNormalize))
    norm_legend_title <- "Z-score Norm. Expr."
    # Use provided zscore_colors. If it's a vector, create a colorRamp2.
    if (is.vector(zscore_colors)) {
      # Create a colorRamp2 function from the vector.
      # A simple approach for a 3-color vector is to map to -1, 0, 1.
      if (length(zscore_colors) == 3) {
        heatmap_color_scale <- circlize::colorRamp2(c(-1, 0, 1), zscore_colors)
      } else {
        # For other vector lengths, just use colorRampPalette
        # This might not give the desired mapping if the vector isn't symmetric around a midpoint
        if(verbose) message("Note: zscore_colors provided as a vector of length != 3. Using colorRampPalette.")
        heatmap_color_scale <- grDevices::colorRampPalette(zscore_colors)(256)
      }
    } else if (inherits(zscore_colors, "function")) {
      # If it's already a color ramp function (like colorRamp2)
      heatmap_color_scale <- zscore_colors
    } else {
      # Default Z-score color scale if zscore_colors is not valid
      if(verbose) message("Invalid zscore_colors provided, using default divergent scale (-1, 0, 1).")
      heatmap_color_scale <- circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
    }

  } else {
    stop("Invalid normalization_method specified. Choose 'unity' or 'zscore'.")
  }

  # Ensure row/column names are preserved after apply
  rownames(heatmap_df_normalized) <- rownames(heatmap_df)
  colnames(heatmap_df_normalized) <- colnames(heatmap_df)

  # Handle potential NAs introduced by normalization (e.g., for genes with constant expression resulting in 0/0 or NA/0)
  # Replace NA with 0, as constant genes are normalized to 0 by the helper functions
  na_count <- sum(is.na(heatmap_df_normalized))
  if (na_count > 0 && verbose) {
    message("Replacing ", na_count, " NA values (likely from constant genes) with 0 after normalization.")
  }
  heatmap_df_normalized[is.na(heatmap_df_normalized)] <- 0

  # --- Create heatmap annotation ---
  # Filter color lists to only include levels present in the data
  valid_cell_colors <- Cell_Colors[names(Cell_Colors) %in% levels(celltype_order)]
  valid_disease_colors <- disease_colors[names(disease_colors) %in% levels(disease_order)]
  valid_subject_colors <- subject_colors[names(subject_colors) %in% levels(subject_order)]

  # Check if all levels in data have colors provided and issue warnings if not
  missing_cell_colors <- setdiff(levels(celltype_order), names(valid_cell_colors))
  missing_disease_colors <- setdiff(levels(disease_order), names(valid_disease_colors))
  missing_subject_colors <- setdiff(levels(subject_order), names(valid_subject_colors))

  if (length(missing_cell_colors) > 0) {
    warning("Missing colors in Cell_Colors for: ", paste(missing_cell_colors, collapse=", "))
  }
  if (length(missing_disease_colors) > 0) {
    warning("Missing colors in disease_colors for: ", paste(missing_disease_colors, collapse=", "))
  }
  if (length(missing_subject_colors) > 0) {
    warning("Missing colors in subject_colors for: ", paste(missing_subject_colors, collapse=", "))
  }

  # Create the annotation list, ensuring names match the annotation data frame columns
  annotation_colors <- list()
  if (length(valid_cell_colors) > 0) annotation_colors[[celltype_column]] <- valid_cell_colors
  if (length(valid_disease_colors) > 0) annotation_colors[[disease_column]] <- valid_disease_colors
  if (length(valid_subject_colors) > 0) annotation_colors[[subject_column]] <- valid_subject_colors

  # Create annotation data frame using the ordered vectors
  annotation_df <- data.frame(
    placeholder = celltype_order # Placeholder to get correct number of rows
  )
  annotation_df[[celltype_column]] <- celltype_order
  annotation_df[[disease_column]] <- disease_order
  annotation_df[[subject_column]] <- subject_order
  annotation_df$placeholder <- NULL # Remove placeholder


  heatmap_top_annotation <- ComplexHeatmap::HeatmapAnnotation(
    df = annotation_df, # Provide the annotation data
    col = annotation_colors, # Provide the color list
    show_legend = c(TRUE, TRUE, FALSE), # Show cell_type, disease, hide subject by default
    # Name legends explicitly using column names
    annotation_legend_param = list(
      title_gp = grid::gpar(fontsize = 10),
      labels_gp = grid::gpar(fontsize = 8)
    ),
    annotation_name_gp = grid::gpar(fontsize = 10),
    simple_anno_size = grid::unit(0.3, "cm") # Adjust annotation bar height
  )

  # --- Create the heatmap object ---
  if (verbose) message("Generating heatmap...")
  ht <- ComplexHeatmap::Heatmap(
    heatmap_df_normalized,
    name = norm_legend_title, # Legend title based on normalization
    col = heatmap_color_scale, # Color scale based on normalization
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_column_names = FALSE,
    show_row_names = TRUE,
    row_names_gp = grid::gpar(fontsize = 8),
    top_annotation = heatmap_top_annotation, # Use the created annotation object
    column_split = celltype_order, # Split columns by cell type
    column_title = NULL, # No titles for splits needed if using annotation
    use_raster = TRUE, # Set to TRUE for large matrices for faster rendering & smaller file size
    # Adjust heatmap legend parameters
    heatmap_legend_param = list(
      title = norm_legend_title,
      legend_height = grid::unit(4, "cm"),
      title_gp = grid::gpar(fontsize = 10),
      labels_gp = grid::gpar(fontsize = 8)
    )
  )

  # --- Save Heatmap ---
  tryCatch({
    if(verbose) message("Saving heatmap to PDF: ", pdf_filename)
    grDevices::pdf(pdf_filename, height = 10, width = 16)
    ComplexHeatmap::draw(ht)
    grDevices::dev.off()
  }, error = function(e) {
    warning("Failed to save heatmap to PDF: ", e$message, call. = FALSE)
    # Ensure device is turned off if error occurred during drawing
    if (names(grDevices::dev.cur()) != "null device") grDevices::dev.off()
  })

  tryCatch({
    if(verbose) message("Saving heatmap to PNG: ", png_filename)
    grDevices::png(filename = png_filename, width = 1600, height = 1000, res = 100)
    ComplexHeatmap::draw(ht)
    grDevices::dev.off()
  }, error = function(e) {
    warning("Failed to save heatmap to PNG: ", e$message, call. = FALSE)
    # Ensure device is turned off if error occurred during drawing
    if (names(grDevices::dev.cur()) != "null device") grDevices::dev.off()
  })


  # Return the heatmap object invisibly
  if (verbose) message("Heatmap generation complete.")
  return(invisible(ht))
}

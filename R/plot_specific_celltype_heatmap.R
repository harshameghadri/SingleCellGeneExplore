# Ensure necessary packages are loaded in your R session before running:
# if (!requireNamespace("Seurat", quietly = TRUE)) install.packages("Seurat")
# if (!requireNamespace("SeuratObject", quietly = TRUE)) install.packages("SeuratObject") # For LayerData, Layers
# if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
# if (!requireNamespace("Matrix", quietly = TRUE)) install.packages("Matrix")
# if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) BiocManager::install("ComplexHeatmap")
# if (!requireNamespace("circlize", quietly = TRUE)) install.packages("circlize")
# if (!requireNamespace("viridis", quietly = TRUE)) install.packages("viridis")
# if (!requireNamespace("RColorBrewer", quietly = TRUE)) install.packages("RColorBrewer") # Optional
# if (!requireNamespace("stats", quietly = TRUE)) install.packages("stats") # for setNames, sd, etc. (usually base)
# if (!requireNamespace("grDevices", quietly = TRUE)) install.packages("grDevices") # for color palettes, pdf, png (usually base)
# if (!requireNamespace("utils", quietly = TRUE)) install.packages("utils") # for head (usually base)


#' Generate a Focused Pseudobulk Heatmap for Specific Cell Types
#'
#' This function creates a heatmap of average gene expression for specified
#' target cell types from a Seurat object. It calculates pseudobulk profiles
#' per subject and disease condition within the target cell types, normalizes
#' the expression, and plots a heatmap using ComplexHeatmap.
#' Compatible with Seurat v4 and v5 object structures (including layers).
#' Rows can be ordered to highlight genes upregulated in one condition (e.g., "PAH")
#' and downregulated in another (e.g., "CTRL").
#'
#' @param seurat_obj A Seurat object containing single-cell RNA sequencing data.
#' @param target_celltypes A character vector specifying the cell type(s) of interest
#'   to include in the heatmap. These must match values in the `celltype_column`.
#' @param genes_for_heatmap A character vector of gene names to be plotted.
#' @param celltype_column A string indicating the metadata column name in
#'   `seurat_obj@meta.data` that contains cell type annotations.
#' @param disease_column A string indicating the metadata column name for disease
#'   status or other primary condition for comparison.
#' @param subject_column A string, the metadata column name for subject or sample
#'   identifiers. Default is "orig.ident".
#' @param assay A string, the Seurat assay to use for gene expression data (e.g.,
#'   "RNA", "SCT"). Default is "RNA".
#' @param min.cells.per.subject An integer, the minimum number of cells required for a
#'   given cell type, disease, and subject combination to be included in
#'   pseudobulk calculation. Default is 3. (Set to 0 to include all regardless of cell count).
#' @param Cell_Colors A named list or vector providing colors for cell types. Names
#'   should match cell type names. If `NULL` or incomplete, colors are
#'   auto-generated. Default is `NULL`.
#' @param disease_colors A named list or vector for disease status colors. Names
#'   should match disease status labels. If `NULL` or incomplete, colors are
#'   auto-generated. Default is `NULL`.
#' @param subject_colors A named list or vector for subject colors. Shown in
#'   annotation if few subjects and colors provided/generated. Default is `NULL`.
#' @param normalization_method A string, either "zscore" (calculates Z-score per
#'   gene across samples) or "unity" (scales expression 0-1 per gene).
#'   Default is "zscore".
#' @param zscore_colors A `circlize::colorRamp2` function or a vector of colors
#'   used for the heatmap when `normalization_method` is "zscore".
#'   Default is `circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))`.
#' @param order_rows_by_condition A logical. If TRUE and `normalization_method` is "zscore",
#'   rows (genes) will be ordered based on their differential expression pattern
#'   (e.g., high in `condition1_for_ordering`, low in `condition2_for_ordering`).
#'   Hierarchical clustering of rows is disabled if this is TRUE. Default is `FALSE`.
#' @param condition1_for_ordering A string. If `order_rows_by_condition` is TRUE, this is
#'   the name of the primary condition (from `disease_column`) whose high expression
#'   should be prioritized at the top of the heatmap (e.g., "PAH").
#'   Required if `order_rows_by_condition` is TRUE. Default is `NULL`.
#' @param condition2_for_ordering A string. If `order_rows_by_condition` is TRUE, this is
#'   the name of the secondary/control condition (from `disease_column`) whose low expression,
#'   in conjunction with high expression in `condition1_for_ordering`,
#'   should be prioritized. Required if `order_rows_by_condition` is TRUE. Default is `NULL`.
#' @param row_names_fontsize A numeric value specifying the font size for row names
#'   (gene names) in the heatmap. Default is 8.
#' @param pdf_filename A string, the filename for saving the heatmap as a PDF.
#'   Set to `NULL` to disable PDF saving. Default is "specific_celltype_heatmap.pdf".
#' @param png_filename A string, the filename for saving the heatmap as a PNG.
#'   Set to `NULL` to disable PNG saving. Default is "specific_celltype_heatmap.png".
#' @param cores An integer, the number of cores to use for parallel processing.
#'   Uses `lapply` if `cores = 1` or on Windows. Default is 1.
#' @param verbose A logical, whether to print progress messages. Default is `TRUE`.
#'
#' @return Invisibly returns a `Heatmap` object generated by `ComplexHeatmap::Heatmap()`.
#'   This object can be drawn manually or further modified.
#'
#' @details
#' The function performs several steps including input validation, gene filtering,
#' pseudobulk calculation, and normalization.
#' If `order_rows_by_condition` is TRUE (and using z-score normalization),
#' genes are sorted to highlight differences between `condition1_for_ordering`
#' and `condition2_for_ordering`. Specifically, it calculates `mean_Zscore(condition1) - mean_Zscore(condition2)`
#' for each gene and sorts genes by this difference in descending order.
#'
#' @examples
#' \dontrun{
#' # Assuming 'pbmc_small' is a Seurat object
#' # Add example metadata
#' pbmc_small$cell_type_annotations <- sample(c("CD4 T", "B cell", "Monocyte"), ncol(pbmc_small), replace = TRUE)
#' pbmc_small$disease_status <- sample(c("ConditionA", "ConditionB"), ncol(pbmc_small), replace = TRUE)
#' pbmc_small$patient_id <- pbmc_small$orig.ident
#' example_genes <- rownames(pbmc_small)[1:20]
#'
#' # Heatmap with custom row ordering and larger row names
#' custom_ordered_heatmap <- plotSpecificCellTypeHeatmap(
#'   seurat_obj = pbmc_small,
#'   target_celltypes = "CD4 T",
#'   genes_for_heatmap = example_genes,
#'   celltype_column = "cell_type_annotations",
#'   disease_column = "disease_status",
#'   subject_column = "patient_id",
#'   normalization_method = "zscore",
#'   order_rows_by_condition = TRUE,
#'   condition1_for_ordering = "ConditionA", # e.g., "PAH"
#'   condition2_for_ordering = "ConditionB", # e.g., "CTRL"
#'   row_names_fontsize = 10, # Custom font size for gene names
#'   min.cells.per.subject = 0,
#'   pdf_filename = NULL,
#'   png_filename = "cd4_custom_order_heatmap.png"
#' )
#' # ComplexHeatmap::draw(custom_ordered_heatmap)
#' }
#'
#' @importFrom Seurat GetAssayData Assays DefaultAssay
#' @importFrom SeuratObject LayerData Layers
#' @importFrom Matrix rowMeans
#' @importFrom dplyr %>% filter arrange sym pull
#' @importFrom parallel mclapply
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation draw
#' @importFrom circlize colorRamp2
#' @importFrom viridis inferno
#' @importFrom grDevices pdf png dev.off rainbow
#' @importFrom stats setNames sd
#' @importFrom utils head
#' @import RColorBrewer
#'
#' @export
plotSpecificCellTypeHeatmap <- function(
    seurat_obj,
    target_celltypes,
    genes_for_heatmap,
    celltype_column,
    disease_column,
    subject_column = "orig.ident",
    assay = "RNA",
    min.cells.per.subject = 3,
    Cell_Colors = NULL,
    disease_colors = NULL,
    subject_colors = NULL,
    normalization_method = "zscore",
    zscore_colors = circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
    order_rows_by_condition = FALSE,
    condition1_for_ordering = NULL,
    condition2_for_ordering = NULL,
    row_names_fontsize = 8, # New parameter for row name font size
    pdf_filename = "specific_celltype_heatmap.pdf",
    png_filename = "specific_celltype_heatmap.png",
    cores = 1,
    verbose = TRUE
) {
  
  # --- Input Validation ---
  if (!inherits(seurat_obj, "Seurat")) {
    stop("seurat_obj must be a Seurat object.")
  }
  if (!is.character(target_celltypes) || length(target_celltypes) == 0) {
    stop("target_celltypes must be a non-empty character vector.")
  }
  if (!is.character(genes_for_heatmap) || length(genes_for_heatmap) == 0) {
    stop("genes_for_heatmap must be a non-empty character vector.")
  }
  if (!is.numeric(min.cells.per.subject) || min.cells.per.subject < 0 || length(min.cells.per.subject) != 1) {
    stop("min.cells.per.subject must be a single non-negative integer.")
  }
  if (order_rows_by_condition && (is.null(condition1_for_ordering) || is.null(condition2_for_ordering))) {
    stop("If order_rows_by_condition is TRUE, condition1_for_ordering and condition2_for_ordering must be specified.")
  }
  if (order_rows_by_condition && normalization_method != "zscore" && verbose) {
    message("Note: order_rows_by_condition is TRUE, but normalization_method is not 'zscore'. ",
            "Ordering will be based on the chosen normalization, which might not be ideal for differential pattern highlighting.")
  }
  if(!is.numeric(row_names_fontsize) || row_names_fontsize <= 0){
    warning("row_names_fontsize must be a positive number. Using default of 8.")
    row_names_fontsize <- 8
  }
  
  
  metadata_cols <- colnames(seurat_obj@meta.data)
  if (!celltype_column %in% metadata_cols) {
    stop("celltype_column '", celltype_column, "' not found in Seurat object metadata.")
  }
  if (!disease_column %in% metadata_cols) {
    stop("disease_column '", disease_column, "' not found in Seurat object metadata.")
  }
  if (!subject_column %in% metadata_cols) {
    stop("subject_column '", subject_column, "' not found in Seurat object metadata.")
  }
  if (!assay %in% Seurat::Assays(seurat_obj)) {
    stop("Assay '", assay, "' not found in Seurat object.")
  }
  
  available_celltypes_in_col <- unique(as.character(seurat_obj@meta.data[[celltype_column]]))
  if (!all(target_celltypes %in% available_celltypes_in_col)) {
    missing_types <- setdiff(target_celltypes, available_celltypes_in_col)
    stop("The following target_celltypes are not found in '", celltype_column, "': ",
         paste(missing_types, collapse = ", "))
  }
  
  assay_obj <- seurat_obj[[assay]]
  gene_matrix_source_type <- "" 
  
  gene_check_src_name <- "data"
  if (gene_check_src_name %in% names(assay_obj)) {
    available_genes <- rownames(Seurat::GetAssayData(seurat_obj, assay = assay, slot = gene_check_src_name))
    gene_matrix_source_type <- "slot"
  } else if (gene_check_src_name %in% SeuratObject::Layers(assay_obj)) {
    if(verbose) message("Accessing '", gene_check_src_name, "' as a layer in assay '", assay, "' for gene availability check (Seurat v5+).")
    available_genes <- rownames(SeuratObject::LayerData(seurat_obj, assay = assay, layer = gene_check_src_name))
    gene_matrix_source_type <- "layer"
  } else {
    gene_check_src_name <- "counts"
    if (gene_check_src_name %in% names(assay_obj)) {
      if(verbose) message("Using 'counts' slot for assay '", assay, "' gene availability check as 'data' was not found.")
      available_genes <- rownames(Seurat::GetAssayData(seurat_obj, assay = assay, slot = gene_check_src_name))
      gene_matrix_source_type <- "slot"
    } else if (gene_check_src_name %in% SeuratObject::Layers(assay_obj)) {
      if(verbose) message("Using 'counts' layer for assay '", assay, "' gene availability check as 'data' layer/slot was not found (Seurat v5+).")
      available_genes <- rownames(SeuratObject::LayerData(seurat_obj, assay = assay, layer = gene_check_src_name))
      gene_matrix_source_type <- "layer"
    } else {
      stop("Neither 'data' nor 'counts' found as a slot or layer in assay '", assay, "' for gene availability check.")
    }
  }
  
  valid_genes <- intersect(genes_for_heatmap, available_genes)
  missing_genes <- setdiff(genes_for_heatmap, available_genes)
  
  if(length(valid_genes) == 0) {
    stop("None of the specified genes_for_heatmap were found in assay '", assay, "' (checked '", gene_check_src_name, "' ", gene_matrix_source_type, ").")
  }
  if(length(missing_genes) > 0 && verbose) {
    message("The following genes were not found (assay '", assay, "', checked '", gene_check_src_name, "' ", gene_matrix_source_type, ") and will be skipped: \n",
            paste(missing_genes, collapse = ", "), "\n")
  }
  genes_for_heatmap <- valid_genes
  
  if(verbose) {
    message("Plotting heatmap for cell type(s): ", paste(target_celltypes, collapse=", "))
    message("Proceeding with ", length(genes_for_heatmap), " genes using assay '", assay, "'.")
    message("Using min.cells.per.subject = ", min.cells.per.subject)
    message("Normalization method: ", normalization_method)
    if(order_rows_by_condition) message("Ordering rows to highlight differences between '", condition1_for_ordering, "' and '", condition2_for_ordering, "'.")
  }
  
  seurat_obj$cellBarcode <- colnames(seurat_obj)
  seurat_obj@meta.data[[celltype_column]] <- factor(seurat_obj@meta.data[[celltype_column]])
  
  all_levels_in_data <- levels(seurat_obj@meta.data[[celltype_column]])
  cellTypes_to_iterate <- intersect(target_celltypes, all_levels_in_data)
  if (length(cellTypes_to_iterate) == 0) {
    stop("No valid target cell types remain after checking against data.")
  }
  if(verbose) message("Will process the following cell types: ", paste(cellTypes_to_iterate, collapse=", "))
  
  relevant_metadata_rows <- seurat_obj@meta.data[[celltype_column]] %in% cellTypes_to_iterate
  meta.data.sub <- seurat_obj@meta.data[relevant_metadata_rows, c(disease_column, celltype_column, subject_column, "cellBarcode"), drop = FALSE]
  
  if(nrow(meta.data.sub) == 0 && verbose){
    message("Warning: No cells found in metadata matching target_celltypes: ", paste(target_celltypes, collapse=", "))
  }
  
  myUnityNormalize <- function(x) {
    if(all(is.na(x))) return(x)
    min_x <- min(x, na.rm = TRUE); max_x <- max(x, na.rm = TRUE)
    range_x = max_x - min_x
    if (is.na(range_x) || range_x == 0) { return(rep(0, length(x))) }
    (x - min_x) / range_x
  }
  myZscoreNormalize <- function(x) {
    if(all(is.na(x))) return(x)
    mean_x <- mean(x, na.rm = TRUE); sd_x <- stats::sd(x, na.rm = TRUE)
    if (is.na(sd_x) || sd_x == 0) { return(rep(0, length(x))) }
    (x - mean_x) / sd_x
  }
  
  get.CT.DS.subj.vector <- function(currentCellType){
    tmp.meta.data.currentCT <- meta.data.sub[meta.data.sub[[celltype_column]] == currentCellType, , drop = FALSE]
    if(nrow(tmp.meta.data.currentCT) == 0) return(vector())
    
    diseases <- unique(as.character(tmp.meta.data.currentCT[[disease_column]]))
    subjects <- unique(as.character(tmp.meta.data.currentCT[[subject_column]]))
    tmp.CT.DS.subj.vector <- vector()
    filtered_out_details <- list()
    
    for(j in 1:length(diseases)){
      for(k in 1:length(subjects)){
        temp.cells.indices <- tmp.meta.data.currentCT[[disease_column]] == diseases[j] &
          tmp.meta.data.currentCT[[subject_column]] == subjects[k]
        temp.cells <- tmp.meta.data.currentCT$cellBarcode[temp.cells.indices]
        
        if (length(temp.cells) >= min.cells.per.subject) {
          tmp.CT.DS.subj.vector <- c(tmp.CT.DS.subj.vector,
                                     paste(currentCellType, diseases[j], subjects[k], sep = "__"))
        } else if (length(temp.cells) > 0) {
          filtered_out_details[[length(filtered_out_details) + 1]] <-
            paste0(currentCellType, "/", diseases[j], "/", subjects[k], " (", length(temp.cells), " cells)")
        }
      }
    }
    if(verbose && length(filtered_out_details) > 0) {
      cat("For ", currentCellType, ": Filtered out ", length(filtered_out_details),
          " combinations with fewer than ", min.cells.per.subject, " cells (e.g., ", filtered_out_details[[1]], ").\n", sep = "")
    }
    if(verbose && length(tmp.CT.DS.subj.vector) > 0) {
      cat("For ", currentCellType, ": Found ", length(tmp.CT.DS.subj.vector), " valid combinations.\n", sep="")
    } else if (verbose && length(tmp.CT.DS.subj.vector) == 0 && nrow(tmp.meta.data.currentCT) > 0){
      cat("For ", currentCellType, ": No valid combinations found meeting min.cells.per.subject criteria.\n", sep="")
    }
    return(tmp.CT.DS.subj.vector)
  }
  
  if (verbose) message("Generating valid combinations using ", cores, " core(s)...")
  if (cores > 1 && .Platform$OS.type != "windows") {
    celltype_disease_subject.list <- parallel::mclapply(cellTypes_to_iterate, get.CT.DS.subj.vector, mc.cores = cores)
  } else {
    if (cores > 1 && .Platform$OS.type == "windows" && verbose) {
      message("Note: mclapply not supported on Windows. Using sequential processing (lapply).")
    }
    celltype_disease_subject.list <- lapply(cellTypes_to_iterate, get.CT.DS.subj.vector)
  }
  celltype_disease_subject <- unlist(celltype_disease_subject.list)
  
  if(length(celltype_disease_subject) == 0) {
    stop("No cell type/disease/subject combinations remain after filtering. Check min.cells.per.subject or data.")
  }
  
  get.SubjectdiseaseCellTypeAvg <- function(combo_string){
    split_result <- strsplit(as.character(combo_string), "__")[[1]]
    temp.cell.type <- split_result[1]; temp.disease <- split_result[2]; temp.subject <- split_result[3]
    
    temp.cells.indices <- seurat_obj@meta.data[[celltype_column]] == temp.cell.type &
      seurat_obj@meta.data[[disease_column]] == temp.disease &
      seurat_obj@meta.data[[subject_column]] == temp.subject
    temp.cells <- seurat_obj@meta.data$cellBarcode[temp.cells.indices]
    temp.cells <- temp.cells[!is.na(temp.cells)]
    
    if (length(temp.cells) == 0) {
      if(verbose) cat("Warning: No cells found for combination '", combo_string, "' during average calculation step. This combo will be NA.\n", sep="")
      tmp.df <- as.data.frame(matrix(NA, nrow = length(genes_for_heatmap), ncol = 1))
      rownames(tmp.df) <- genes_for_heatmap; colnames(tmp.df) <- combo_string
      return(tmp.df)
    }
    
    assay_obj_expr <- seurat_obj[[assay]]
    expression_matrix <- NULL
    expr_data_src_name <- "data" 
    
    if (expr_data_src_name %in% names(assay_obj_expr)) {
      expression_matrix <- Seurat::GetAssayData(seurat_obj, assay = assay, slot = expr_data_src_name)
    } else if (expr_data_src_name %in% SeuratObject::Layers(assay_obj_expr)) {
      if(verbose && !exists("warned_expr_layer_data", envir = .GlobalEnv)) {
        message("Note: Accessing '", expr_data_src_name, "' as a layer in assay '", assay, "' for expression (Seurat v5+).")
        assign("warned_expr_layer_data", TRUE, envir = .GlobalEnv)
      }
      expression_matrix <- SeuratObject::LayerData(seurat_obj, assay = assay, layer = expr_data_src_name)
    } else {
      expr_data_src_name <- "counts" 
      if (expr_data_src_name %in% names(assay_obj_expr)) {
        if(verbose && !exists("warned_expr_slot_counts", envir = .GlobalEnv)) {
          message("Note: Using 'counts' slot from assay '", assay, "' for expression as 'data' not found.")
          assign("warned_expr_slot_counts", TRUE, envir = .GlobalEnv)
        }
        expression_matrix <- Seurat::GetAssayData(seurat_obj, assay = assay, slot = expr_data_src_name)
      } else if (expr_data_src_name %in% SeuratObject::Layers(assay_obj_expr)) {
        if(verbose && !exists("warned_expr_layer_counts", envir = .GlobalEnv)) {
          message("Note: Using 'counts' layer from assay '", assay, "' for expression as 'data' not found (Seurat v5+).")
          assign("warned_expr_layer_counts", TRUE, envir = .GlobalEnv)
        }
        expression_matrix <- SeuratObject::LayerData(seurat_obj, assay = assay, layer = expr_data_src_name)
      } else {
        stop("Neither 'data' nor 'counts' found as a slot or layer in assay '", assay, "' to get expression values.")
      }
    }
    
    genes_to_fetch <- intersect(genes_for_heatmap, rownames(expression_matrix))
    if(length(genes_to_fetch) == 0){ 
      if(verbose) message("Warning: No valid genes to fetch for combo ", combo_string)
      tmp.df <- as.data.frame(matrix(NA, nrow = length(genes_for_heatmap), ncol = 1))
      rownames(tmp.df) <- genes_for_heatmap; colnames(tmp.df) <- combo_string
      return(tmp.df)
    }
    
    if (length(temp.cells) > 1) {
      tmp.df <- as.data.frame(Matrix::rowMeans(expression_matrix[genes_to_fetch, temp.cells, drop = FALSE]))
    } else { 
      tmp.df <- as.data.frame(expression_matrix[genes_to_fetch, temp.cells, drop = FALSE])
    }
    colnames(tmp.df) <- combo_string
    return(tmp.df)
  }
  
  if (verbose) message("Calculating pseudo-bulk profiles for ", length(celltype_disease_subject) ," combinations using ", cores, " core(s)...")
  if (cores > 1 && .Platform$OS.type != "windows") {
    collapsed.mtx.list <- parallel::mclapply(celltype_disease_subject, get.SubjectdiseaseCellTypeAvg, mc.cores = cores)
  } else {
    collapsed.mtx.list <- lapply(celltype_disease_subject, get.SubjectdiseaseCellTypeAvg)
  }
  if(exists("warned_expr_layer_data", envir = .GlobalEnv)) rm("warned_expr_layer_data", envir = .GlobalEnv)
  if(exists("warned_expr_slot_counts", envir = .GlobalEnv)) rm("warned_expr_slot_counts", envir = .GlobalEnv)
  if(exists("warned_expr_layer_counts", envir = .GlobalEnv)) rm("warned_expr_layer_counts", envir = .GlobalEnv)
  
  collapsed.SubjectdiseaseCellTypeAvg.mtx <- do.call(cbind, collapsed.mtx.list)
  
  if (ncol(collapsed.SubjectdiseaseCellTypeAvg.mtx) > 0 && nrow(collapsed.SubjectdiseaseCellTypeAvg.mtx) > 0) {
    all_na_columns <- apply(collapsed.SubjectdiseaseCellTypeAvg.mtx, 2, function(col) all(is.na(col)))
    if (any(all_na_columns)) {
      num_removed <- sum(all_na_columns)
      if (verbose) {
        message("Removing ", num_removed, " combinations that resulted in all NA expression values (likely no cells found at averaging step).")
      }
      collapsed.SubjectdiseaseCellTypeAvg.mtx <- collapsed.SubjectdiseaseCellTypeAvg.mtx[, !all_na_columns, drop = FALSE]
      celltype_disease_subject <- celltype_disease_subject[!all_na_columns]
    }
  }
  
  if (ncol(collapsed.SubjectdiseaseCellTypeAvg.mtx) == 0) {
    stop("Pseudo-bulk matrix is empty after removing NA columns. All combinations might have had no cells or no expression data.")
  }
  
  if(verbose) {
    message("Final pseudo-bulk matrix for ", ncol(collapsed.SubjectdiseaseCellTypeAvg.mtx), " combinations.")
    message("Dimensions of pseudo-bulk matrix: ", paste(dim(collapsed.SubjectdiseaseCellTypeAvg.mtx), collapse = " x "))
  }
  
  heatmap_metadata <- data.frame(
    cell.ident = colnames(collapsed.SubjectdiseaseCellTypeAvg.mtx),
    temp_cell_type = sapply(strsplit(colnames(collapsed.SubjectdiseaseCellTypeAvg.mtx), "__"), `[`, 1),
    temp_disease = sapply(strsplit(colnames(collapsed.SubjectdiseaseCellTypeAvg.mtx), "__"), `[`, 2),
    temp_subject = sapply(strsplit(colnames(collapsed.SubjectdiseaseCellTypeAvg.mtx), "__"), `[`, 3),
    stringsAsFactors = FALSE
  )
  colnames(heatmap_metadata) <- c("cell.ident", celltype_column, disease_column, subject_column)
  
  actual_celltypes_in_heatmap <- unique(heatmap_metadata[[celltype_column]])
  ordered_celltype_levels <- intersect(target_celltypes, actual_celltypes_in_heatmap)
  if(length(ordered_celltype_levels) == 0 && length(actual_celltypes_in_heatmap) > 0){
    ordered_celltype_levels <- sort(actual_celltypes_in_heatmap)
    if(verbose) message("Warning: Target celltypes not in final heatmap data; ordering by sorted unique cell types found.")
  } else if (length(ordered_celltype_levels) == 0 && length(actual_celltypes_in_heatmap) == 0){
    stop("No cell types in final heatmap data to order by.")
  }
  heatmap_metadata[[celltype_column]] <- factor(heatmap_metadata[[celltype_column]], levels = ordered_celltype_levels)
  
  disease_levels_present <- sort(unique(heatmap_metadata[[disease_column]]))
  heatmap_metadata[[disease_column]] <- factor(heatmap_metadata[[disease_column]], levels = disease_levels_present)
  
  if (order_rows_by_condition) {
    if (!condition1_for_ordering %in% disease_levels_present) {
      stop("condition1_for_ordering '", condition1_for_ordering, "' not found in the disease levels of the data: ", paste(disease_levels_present, collapse=", "))
    }
    if (!condition2_for_ordering %in% disease_levels_present) {
      stop("condition2_for_ordering '", condition2_for_ordering, "' not found in the disease levels of the data: ", paste(disease_levels_present, collapse=", "))
    }
  }
  
  subject_levels_present <- sort(unique(heatmap_metadata[[subject_column]]))
  heatmap_metadata[[subject_column]] <- factor(heatmap_metadata[[subject_column]], levels = subject_levels_present)
  
  tryCatch({
    heatmap_metadata_ordered <- heatmap_metadata %>%
      dplyr::arrange(!!dplyr::sym(celltype_column), !!dplyr::sym(disease_column), !!dplyr::sym(subject_column))
  }, error = function(e){
    if(verbose) message("dplyr::arrange failed, using base R order. Error: ", e$message)
    heatmap_metadata_ordered <- heatmap_metadata[order(heatmap_metadata[[celltype_column]],
                                                       heatmap_metadata[[disease_column]],
                                                       heatmap_metadata[[subject_column]]), ]
  })
  
  cell_order <- heatmap_metadata_ordered$cell.ident
  celltype_order_factor <- heatmap_metadata_ordered[[celltype_column]]
  disease_order_factor <- heatmap_metadata_ordered[[disease_column]] 
  subject_order_factor <- heatmap_metadata_ordered[[subject_column]]
  
  genes_present_in_matrix <- intersect(genes_for_heatmap, rownames(collapsed.SubjectdiseaseCellTypeAvg.mtx))
  if(length(genes_present_in_matrix) == 0) {
    stop("No genes left to plot in the heatmap matrix.")
  }
  heatmap_df <- as.matrix(collapsed.SubjectdiseaseCellTypeAvg.mtx[genes_present_in_matrix, cell_order, drop = FALSE])
  
  if (verbose) message("Applying ", normalization_method, " normalization row-wise...")
  if (normalization_method == "unity") {
    heatmap_df_normalized <- t(apply(heatmap_df, MARGIN = 1, FUN = myUnityNormalize))
    norm_legend_title <- "Unity Norm. Expr."
    heatmap_color_scale <- viridis::inferno(256)
  } else if (normalization_method == "zscore") {
    heatmap_df_normalized <- t(apply(heatmap_df, MARGIN = 1, FUN = myZscoreNormalize))
    norm_legend_title <- "Z-score"
    if (is.function(zscore_colors)) {
      heatmap_color_scale <- zscore_colors
    } else if (is.vector(zscore_colors) && length(zscore_colors) >= 2) {
      if(length(zscore_colors) == 3 && all(sapply(zscore_colors, is.character))) {
        default_breaks <- tryCatch(environment(formals(plotSpecificCellTypeHeatmap)$zscore_colors)$breaks, error = function(e) c(-2,0,2))
        if(is.null(default_breaks) || length(default_breaks) != 3) default_breaks <- c(-2,0,2)
        heatmap_color_scale <- circlize::colorRamp2(default_breaks, zscore_colors)
      } else {
        heatmap_color_scale <- grDevices::colorRampPalette(zscore_colors)(256)
      }
    } else {
      if(verbose) message("Invalid zscore_colors. Using default blue-white-red (-2,0,2).")
      heatmap_color_scale <- circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
    }
  } else {
    stop("Invalid normalization_method. Choose 'unity' or 'zscore'.")
  }
  rownames(heatmap_df_normalized) <- rownames(heatmap_df); colnames(heatmap_df_normalized) <- colnames(heatmap_df)
  na_count <- sum(is.na(heatmap_df_normalized))
  if (na_count > 0 && verbose) {
    message("Replacing ", na_count, " NAs (e.g. from constant genes) with 0 after normalization.")
  }
  heatmap_df_normalized[is.na(heatmap_df_normalized)] <- 0
  
  cluster_rows_param <- TRUE 
  row_dend_width_param <- grid::unit(15, "mm")
  
  if (order_rows_by_condition) {
    if (verbose) message("Applying custom row ordering based on condition difference...")
    cluster_rows_param <- FALSE 
    row_dend_width_param <- grid::unit(0, "mm") 
    
    cols_cond1 <- colnames(heatmap_df_normalized)[disease_order_factor == condition1_for_ordering]
    cols_cond2 <- colnames(heatmap_df_normalized)[disease_order_factor == condition2_for_ordering]
    
    if (length(cols_cond1) == 0) warning("No columns found for condition1_for_ordering: '", condition1_for_ordering, "'. Row ordering might not be effective.")
    if (length(cols_cond2) == 0) warning("No columns found for condition2_for_ordering: '", condition2_for_ordering, "'. Row ordering might not be effective.")
    
    mean_expr_cond1 <- Matrix::rowMeans(heatmap_df_normalized[, cols_cond1, drop = FALSE], na.rm = TRUE)
    mean_expr_cond2 <- Matrix::rowMeans(heatmap_df_normalized[, cols_cond2, drop = FALSE], na.rm = TRUE)
    
    diff_metric <- mean_expr_cond1 - mean_expr_cond2
    
    ordered_gene_names <- names(sort(diff_metric, decreasing = TRUE, na.last = TRUE))
    
    heatmap_df_normalized <- heatmap_df_normalized[ordered_gene_names, , drop = FALSE]
  }
  
  generate_colors_for_annotation <- function(levels_vector, provided_colors, type_name) {
    if (length(levels_vector) == 0) return(stats::setNames(character(0), character(0)))
    final_colors_map <- stats::setNames(rep(NA_character_, length(levels_vector)), levels_vector)
    if (!is.null(provided_colors)) {
      for (level_name in names(provided_colors)) {
        if (level_name %in% levels_vector) {
          final_colors_map[level_name] <- provided_colors[[level_name]]
        }
      }
    }
    levels_needing_color <- names(final_colors_map[is.na(final_colors_map)])
    if (length(levels_needing_color) > 0) {
      if (verbose) {
        message("Auto-generating colors for ", length(levels_needing_color), " '", type_name, "' level(s): ",
                paste(utils::head(levels_needing_color,3), collapse=", "), if(length(levels_needing_color)>3) "...") 
      }
      auto_gen_colors <- character(0); num_to_gen <- length(levels_needing_color)
      if (requireNamespace("RColorBrewer", quietly = TRUE) && num_to_gen > 0) {
        pal_name <- "Set2"; max_cols <- RColorBrewer::brewer.pal.info[pal_name, "maxcolors"]
        if (num_to_gen <= max_cols && num_to_gen >=3 ) {
          auto_gen_colors <- RColorBrewer::brewer.pal(num_to_gen, pal_name)
        } else { 
          pal_name <- "Paired"; max_cols <- RColorBrewer::brewer.pal.info[pal_name, "maxcolors"]
          if (num_to_gen <= max_cols && num_to_gen >=3) {
            auto_gen_colors <- RColorBrewer::brewer.pal(num_to_gen, pal_name)
          } else if (num_to_gen < 3 && num_to_gen > 0) {
            auto_gen_colors <- RColorBrewer::brewer.pal(3, "Set2")[1:num_to_gen]
          }
        }
      }
      if (length(auto_gen_colors) != num_to_gen && num_to_gen > 0) {
        if (num_to_gen == 1) auto_gen_colors <- "grey50"
        else if (num_to_gen == 2) auto_gen_colors <- c("skyblue", "tomato")
        else auto_gen_colors <- grDevices::rainbow(num_to_gen)
      }
      names(auto_gen_colors) <- levels_needing_color
      final_colors_map[levels_needing_color] <- auto_gen_colors
    }
    return(final_colors_map[levels_vector])
  }
  
  ct_levels_for_color <- levels(celltype_order_factor)
  ds_levels_for_color <- levels(disease_order_factor)
  sb_levels_for_color <- levels(subject_order_factor)
  
  final_cell_colors    <- generate_colors_for_annotation(ct_levels_for_color, Cell_Colors, celltype_column)
  final_disease_colors <- generate_colors_for_annotation(ds_levels_for_color, disease_colors, disease_column)
  final_subject_colors <- generate_colors_for_annotation(sb_levels_for_color, subject_colors, subject_column)
  
  annotation_df_for_hm <- data.frame(row.names = cell_order) 
  if(length(ct_levels_for_color) > 0) annotation_df_for_hm[[celltype_column]] <- celltype_order_factor
  if(length(ds_levels_for_color) > 0) annotation_df_for_hm[[disease_column]] <- disease_order_factor
  
  annotation_colors_list <- list()
  if(length(final_cell_colors) > 0 && celltype_column %in% names(annotation_df_for_hm)) annotation_colors_list[[celltype_column]] <- final_cell_colors
  if(length(final_disease_colors) > 0 && disease_column %in% names(annotation_df_for_hm)) annotation_colors_list[[disease_column]] <- final_disease_colors
  
  show_subject_annotation_bar <- subject_column %in% names(heatmap_metadata_ordered) && length(sb_levels_for_color) > 0 && length(sb_levels_for_color) <= 20 
  if(show_subject_annotation_bar) { 
    annotation_df_for_hm[[subject_column]] <- subject_order_factor
    if(length(final_subject_colors) > 0) annotation_colors_list[[subject_column]] <- final_subject_colors
  }
  
  legends_to_show_map <- list()
  if(celltype_column %in% names(annotation_df_for_hm)) legends_to_show_map[[celltype_column]] <- (length(ct_levels_for_color) > 1)
  if(disease_column %in% names(annotation_df_for_hm)) legends_to_show_map[[disease_column]] <- (length(ds_levels_for_color) > 1)
  if(subject_column %in% names(annotation_df_for_hm) && show_subject_annotation_bar) legends_to_show_map[[subject_column]] <- (length(sb_levels_for_color) > 1)
  
  legends_to_show_vec <- unlist(legends_to_show_map[colnames(annotation_df_for_hm)])
  legends_to_show_vec[is.na(legends_to_show_vec)] <- FALSE
  
  heatmap_top_annotation <- ComplexHeatmap::HeatmapAnnotation(
    df = annotation_df_for_hm,
    col = annotation_colors_list,
    show_legend = legends_to_show_vec,
    annotation_legend_param = list(
      title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
      labels_gp = grid::gpar(fontsize = 8)
    ),
    annotation_name_gp = grid::gpar(fontsize = 10, fontface = "bold"),
    simple_anno_size = grid::unit(0.4, "cm"),
    gap = grid::unit(1, "mm")
  )
  
  if (verbose) message("Generating heatmap...")
  
  col_split_factor <- NULL 
  column_title_for_split <- NULL
  
  if (celltype_column %in% names(annotation_df_for_hm)) {
    current_col_split_factor <- annotation_df_for_hm[[celltype_column]]
    if(is.factor(current_col_split_factor) && length(levels(current_col_split_factor)) > 1) {
      col_split_factor <- current_col_split_factor
      column_title_for_split <- levels(col_split_factor)
    } 
  }
  
  ht <- ComplexHeatmap::Heatmap(
    heatmap_df_normalized, 
    name = norm_legend_title,
    col = heatmap_color_scale,
    cluster_rows = cluster_rows_param, 
    row_dend_width = row_dend_width_param, 
    cluster_columns = FALSE,
    show_column_names = FALSE,
    show_row_names = TRUE,
    row_names_gp = grid::gpar(fontsize = row_names_fontsize), # Use new font size parameter
    top_annotation = heatmap_top_annotation,
    column_split = col_split_factor,
    column_title = column_title_for_split,
    column_title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
    column_gap = grid::unit(2, "mm"),
    use_raster = TRUE, raster_quality = 4,
    heatmap_legend_param = list(
      title = norm_legend_title,
      legend_height = grid::unit(4, "cm"),
      title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
      labels_gp = grid::gpar(fontsize = 8)
    )
  )
  
  draw_heatmap_safely <- function(ht_obj, device_func, filename, ...) {
    tryCatch({
      actual_matrix <- ht_obj@matrix 
      
      if (deparse(substitute(device_func)) == "grDevices::pdf"){
        h_inches <- max(6, nrow(actual_matrix) * 0.12 + 5) 
        w_inches <- max(8, ncol(actual_matrix) * 0.15 + 6)
        device_func(filename, width = w_inches, height = h_inches, ...)
      } else if (deparse(substitute(device_func)) == "grDevices::png") {
        h_px <- max(600, nrow(actual_matrix) * 15 + 500)
        w_px <- max(800, ncol(actual_matrix) * 25 + 600)
        device_func(filename, width = w_px, height = h_px, ...)
      } else { 
        device_func(filename, ...) 
      }
      
      ComplexHeatmap::draw(ht_obj, heatmap_legend_side = "right", annotation_legend_side = "right")
      grDevices::dev.off()
      if(verbose) message("Heatmap successfully saved to: ", filename)
    }, error = function(e) {
      warning("Failed to save heatmap to ", filename, ": ", e$message, call. = FALSE)
      if (grDevices::dev.cur() != 1 && names(grDevices::dev.cur())[1] != "null device") { 
        grDevices::dev.off()
      }
    })
  }
  
  if (!is.null(pdf_filename)) {
    draw_heatmap_safely(ht, grDevices::pdf, pdf_filename)
  }
  if (!is.null(png_filename)) {
    draw_heatmap_safely(ht, grDevices::png, png_filename, res = 100)
  }
  
  if (verbose) message("Focused heatmap generation complete. Returning heatmap object.")
  return(invisible(ht))
}


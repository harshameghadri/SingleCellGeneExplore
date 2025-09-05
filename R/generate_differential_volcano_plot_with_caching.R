#' Generate Differential Expression Volcano Plot(s) with Caching
#'
#' This function performs differential gene expression analysis using Seurat's FindMarkers
#' for specified cell type(s) and generates volcano plot(s). It includes a mechanism
#' to cache FindMarkers results in the global environment to avoid re-computation.
#' Plot saving attempts cowplot::save_plot, then ggplot2::ggsave, then base R devices.
#' Allows customization of label font size, connector arrowhead type, x-axis breaks,
#' legend title, and parameters for label repulsion. Uses ggplot2::theme_classic().
#' Note: `repel_force` and `repel_box_padding` require EnhancedVolcano version >= 1.7.3.
#'
#' @param seurat_obj A Seurat object.
#' @param cell_type_to_analyze Character vector or NULL. Specifies cell type(s).
#' @param metadataColumn Character string. Metadata column for cell type annotations.
#'        Default is "cell.type.ident".
#' @param disease_ident_column Character string. Metadata column for DE comparison identities.
#'        Default is "disease.ident".
#' @param ident_1 Character string. First identity group for comparison. Default is "PAH".
#' @param ident_2 Character string. Second identity group for comparison. Default is "CTRL".
#' @param assay Character string. Assay for `FindMarkers`. Default is "RNA".
#' @param logfc_threshold_dm Numeric. `logfc.threshold` for `FindMarkers`. Default is 0.25.
#' @param min_pct Numeric. `min.pct` for `FindMarkers`. Default is 0.1.
#' @param test_use Character string. Test for `FindMarkers`. Default is "wilcox".
#' @param genes_to_label Character vector or NULL. Genes to label on the plot.
#' @param top_n_genes Integer. Number of top up/down genes to label if `genes_to_label` is NULL.
#'        Default is 10.
#' @param lab_font_size Numeric. Font size for gene labels. Default is 3.0.
#' @param connector_arrow_type Character. Type of arrowhead for connectors. Default is "open".
#' @param connector_length Numeric. Length of the connector lines from points to labels,
#'        interpreted in "lines" units (e.g., relative to text line height).
#'        Passed to `EnhancedVolcano`'s `lengthConnectors` argument after conversion to a grid::unit object.
#'        Default is 0.5.
#' @param repel_min_segment_length Numeric. Minimum segment length for ggrepel.
#'        Passed to `EnhancedVolcano`'s `min.segment.length` argument. Default is 0.1.
#' @param repel_force Numeric. Force of repulsion (requires EnhancedVolcano >= 1.7.3). Default is 1.0.
#'        If using an older version, this parameter will be ignored by EnhancedVolcano.
#' @param repel_box_padding Numeric. Padding around label bounding box (requires EnhancedVolcano >= 1.7.3). Default is 0.25.
#'        If using an older version, this parameter will be ignored by EnhancedVolcano.
#' @param x_axis_breaks Numeric vector or NULL. Specifies breaks for the x-axis. Default is NULL.
#' @param legend_title Character string. Title for the plot legend. Default is "Regulation".
#' @param volcpath Character string. Base directory for saving outputs.
#'        Default is "./VolcanoPlots_Cached".
#' @param fc_cutoff_plot Numeric. Absolute Log2FC cutoff for lines and coloring. Default is 1.
#' @param p_cutoff_plot Numeric. P-value cutoff line on plot. Default is 0.05.
#' @param color_up_regulated Character string. Color for up-regulated genes. Default is "gold".
#' @param color_down_regulated Character string. Color for down-regulated genes. Default is "royalblue".
#' @param color_non_significant Character string. Color for non-significant genes. Default is "grey".
#' @param save_deg_results Logical. Whether to save the DEG table. Default is TRUE.
#' @param use_cache Logical. Whether to use/save cached DE results. Default is TRUE.
#' @param force_recompute Logical. If TRUE, re-run `FindMarkers` even if cache exists. Default is FALSE.
#' @param cache_env Environment. Environment for caching. Default is `.GlobalEnv`.
#' @param png_width Numeric. Width for PNG plots in inches. Default is 12.
#' @param png_height Numeric. Height for PNG plots in inches. Default is 9.
#' @param png_dpi Numeric. DPI for PNG plots. Default is 300.
#' @param pdf_width Numeric. Width for PDF plots in inches. Default is 12.
#' @param pdf_height Numeric. Height for PDF plots in inches. Default is 9.
#'
#' @return Invisibly returns a list with paths to outputs and status for each cell type.
#'
#' @importFrom Seurat Idents<- FindMarkers DefaultAssay
#' @importFrom SeuratObject Assays Layers GetAssayData
#' @importFrom dplyr as_tibble select filter arrange desc slice_head pull everything
#' @importFrom EnhancedVolcano EnhancedVolcano
#' @importFrom ggplot2 ggsave theme_classic scale_x_continuous waiver labs
#' @importFrom cowplot save_plot
#' @importFrom grDevices png dev.off pdf
#' @importFrom grid unit # Import unit function from grid package
#' @importFrom utils write.table
#' @importFrom stringr str_wrap
#' @importFrom stats na.omit
#'
#' @export
generate_differential_volcano_plot_with_caching <- function(
    seurat_obj,
    cell_type_to_analyze = NULL,
    metadataColumn = "cell.type.ident",
    disease_ident_column = "disease.ident",
    ident_1 = "PAH",
    ident_2 = "CTRL",
    assay = "RNA",
    logfc_threshold_dm = 0.25,
    min_pct = 0.1,
    test_use = "wilcox",
    genes_to_label = NULL,
    top_n_genes = 10,
    lab_font_size = 3.0,
    connector_arrow_type = "open",
    connector_length = 0.5,
    repel_min_segment_length = 0.1,
    repel_force = 1.0, 
    repel_box_padding = 0.25, 
    x_axis_breaks = NULL,
    legend_title = "Regulation",
    volcpath = "./VolcanoPlots_Cached",
    fc_cutoff_plot = 1,
    p_cutoff_plot = 0.05,
    color_up_regulated = "gold",
    color_down_regulated = "royalblue",
    color_non_significant = "grey",
    save_deg_results = TRUE,
    use_cache = TRUE,
    force_recompute = FALSE,
    cache_env = .GlobalEnv,
    png_width = 12,
    png_height = 9,
    png_dpi = 300,
    pdf_width = 12,
    pdf_height = 9
) {
  
  # --- Input Validation & Package Checks ---
  if (!requireNamespace("Seurat", quietly = TRUE) || !requireNamespace("SeuratObject", quietly = TRUE)) stop("Packages 'Seurat' and 'SeuratObject' are required.")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' is required.")
  if (!requireNamespace("EnhancedVolcano", quietly = TRUE)) stop("Package 'EnhancedVolcano' is required.")
  if (!requireNamespace("stringr", quietly = TRUE)) stop("Package 'stringr' is required.")
  if (!requireNamespace("utils", quietly = TRUE)) stop("Package 'utils' is required.")
  if (!requireNamespace("stats", quietly = TRUE)) stop("Package 'stats' is required for na.omit.")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' (for ggsave, themes, scales, labs) is required.")
  if (!requireNamespace("cowplot", quietly = TRUE)) {
    warning("Package 'cowplot' not found. cowplot::save_plot will not be available as a primary saving method.")
  }
  if (!requireNamespace("grDevices", quietly = TRUE)) stop("Package 'grDevices' (for png, pdf, dev.off) is required.")
  if (!requireNamespace("grid", quietly = TRUE)) stop("Package 'grid' (for unit object) is required.") # Added grid check
  
  if (!is.logical(save_deg_results)) stop("`save_deg_results` must be a single logical value.")
  if (!is.logical(use_cache)) stop("`use_cache` must be a single logical value.")
  if (!is.logical(force_recompute)) stop("`force_recompute` must be a single logical value.")
  if (!is.environment(cache_env)) stop("`cache_env` must be an environment.")
  if (!is.numeric(lab_font_size) || lab_font_size <= 0) stop("`lab_font_size` must be a positive number.")
  if (!is.character(connector_arrow_type) || length(connector_arrow_type) != 1) stop("`connector_arrow_type` must be a single character string.")
  if (!is.numeric(fc_cutoff_plot) || fc_cutoff_plot < 0) stop("`fc_cutoff_plot` must be a non-negative number.")
  if (!is.null(x_axis_breaks) && !is.numeric(x_axis_breaks)) stop("`x_axis_breaks` must be NULL or a numeric vector.")
  if (!is.character(legend_title) || length(legend_title) != 1) stop("`legend_title` must be a single character string.")
  if (!is.numeric(connector_length) || connector_length < 0) stop("`connector_length` must be a non-negative number.")
  if (!is.numeric(repel_min_segment_length) || repel_min_segment_length < 0) stop("`repel_min_segment_length` must be a non-negative number.")
  if (!is.numeric(repel_force) || repel_force < 0) stop("`repel_force` must be a non-negative number.")
  if (!is.numeric(repel_box_padding) || repel_box_padding < 0) stop("`repel_box_padding` must be a non-negative number.")
  
  
  if (!metadataColumn %in% colnames(seurat_obj@meta.data)) stop(paste0("Metadata column '", metadataColumn, "' not found."))
  if (!disease_ident_column %in% colnames(seurat_obj@meta.data)) stop(paste0("Disease identity column '", disease_ident_column, "' not found."))
  
  unique_disease_idents <- unique(as.character(seurat_obj@meta.data[[disease_ident_column]]))
  if (!ident_1 %in% unique_disease_idents) stop(paste0("ident_1 '", ident_1, "' not found. Available: ", paste(unique_disease_idents, collapse=", ")))
  if (!ident_2 %in% unique_disease_idents) stop(paste0("ident_2 '", ident_2, "' not found. Available: ", paste(unique_disease_idents, collapse=", ")))
  
  all_available_cell_types <- unique(as.character(stats::na.omit(seurat_obj@meta.data[[metadataColumn]])))
  cell_types_to_process <- c()
  if (is.null(cell_type_to_analyze)) {
    message(paste0("Processing all unique cell types in '", metadataColumn, "'."))
    cell_types_to_process <- all_available_cell_types
    if (length(cell_types_to_process) == 0) stop(paste0("No cell types found in '", metadataColumn, "'."))
  } else {
    if (!is.character(cell_type_to_analyze)) stop("`cell_type_to_analyze` must be NULL or character vector.")
    missing_types <- setdiff(cell_type_to_analyze, all_available_cell_types)
    if (length(missing_types) > 0) warning(paste0("Skipping non-existent cell types: ", paste(missing_types, collapse = ", ")))
    cell_types_to_process <- intersect(cell_type_to_analyze, all_available_cell_types)
    if (length(cell_types_to_process) == 0) { message("No valid cell types to process."); return(invisible(NULL)) }
  }
  message(paste0("Will process: ", paste(cell_types_to_process, collapse = ", ")))
  all_results <- list()
  Seurat::DefaultAssay(seurat_obj) <- assay
  
  for (current_cell_type in cell_types_to_process) {
    message(paste0("\nProcessing cell type: ", current_cell_type))
    current_result_paths <- list(cell_type = current_cell_type, deg_table = NULL, png_plot = NULL, pdf_plot = NULL, status = "Skipped", cache_status = "Cache not used")
    sanitized_current_cell_type_for_name <- gsub("[^a-zA-Z0-9_.-]", "_", current_cell_type)
    sanitized_ident_1 <- gsub("[^a-zA-Z0-9_.-]", "_", ident_1)
    sanitized_ident_2 <- gsub("[^a-zA-Z0-9_.-]", "_", ident_2)
    sanitized_assay <- gsub("[^a-zA-Z0-9_.-]", "_", assay)
    cache_object_name <- paste("DECache", sanitized_current_cell_type_for_name, sanitized_ident_1, "vs", sanitized_ident_2, sanitized_assay, sep = "_")
    marker_diff <- NULL
    
    if (use_cache && !force_recompute && exists(cache_object_name, envir = cache_env)) {
      message(paste0("  Loading cached DE results '", cache_object_name, "'."))
      cached_data <- get(cache_object_name, envir = cache_env)
      if (is.data.frame(cached_data) && "avg_log2FC" %in% colnames(cached_data) && "p_val" %in% colnames(cached_data)) {
        marker_diff <- cached_data; current_result_paths$cache_status <- "Loaded from cache"
      } else {
        warning(paste0("  Cached '", cache_object_name, "' invalid. Recomputing."), call. = FALSE)
        current_result_paths$cache_status <- "Cache invalid, recomputing"
      }
    } else {
      current_result_paths$cache_status <- if (use_cache && force_recompute) "Cache ignored (force_recompute), recomputing"
      else if (use_cache) "Cache miss, computing"
      else "Cache disabled, computing"
      message(paste0("  ", current_result_paths$cache_status))
    }
    
    tryCatch({
      if (is.null(marker_diff)) {
        subset_obj <- subset(seurat_obj, cells = colnames(seurat_obj)[seurat_obj@meta.data[[metadataColumn]] == current_cell_type])
        if (ncol(subset_obj) == 0) { warning(paste0("No cells for '", current_cell_type, "'. Skipping DE."), call. = FALSE); all_results[[current_cell_type]] <- current_result_paths; next }
        Seurat::Idents(subset_obj) <- disease_ident_column
        ident_1_counts <- sum(Seurat::Idents(subset_obj) == ident_1); ident_2_counts <- sum(Seurat::Idents(subset_obj) == ident_2)
        if (ident_1_counts < 3 || ident_2_counts < 3) warning(paste0("Few cells for comparison in '", current_cell_type, "'. ", ident_1, "=", ident_1_counts, ", ", ident_2, "=", ident_2_counts,". Results may be unreliable."), call. = FALSE)
        if (ident_1_counts == 0 || ident_2_counts == 0) { warning(paste0("Not enough cells for comparison in '", current_cell_type, "'. Skipping DE."), call. = FALSE); all_results[[current_cell_type]] <- current_result_paths; next }
        message(paste0("  Performing FindMarkers for ", current_cell_type, ": ", ident_1, " vs ", ident_2))
        marker_diff_computed <- Seurat::FindMarkers(subset_obj, ident.1 = ident_1, ident.2 = ident_2, assay = assay, logfc.threshold = logfc_threshold_dm, test.use = test_use, min.pct = min_pct, verbose = FALSE)
        if (use_cache) {
          assign(cache_object_name, marker_diff_computed, envir = cache_env)
          message(paste0("  DE results saved to cache as '", cache_object_name, "'."))
          current_result_paths$cache_status <- paste(current_result_paths$cache_status, "and saved to cache")
        }
        marker_diff <- marker_diff_computed
      }
      if (is.null(marker_diff) || nrow(marker_diff) == 0) { warning(paste0("  No DEGs for '", current_cell_type, "'. Skipping plot."), call. = FALSE); all_results[[current_cell_type]] <- current_result_paths; next }
      if (!"genes" %in% colnames(marker_diff)) marker_diff$genes <- rownames(marker_diff)
      marker_diff <- dplyr::as_tibble(marker_diff)
      if ("genes" %in% colnames(marker_diff)) marker_diff <- dplyr::select(marker_diff, genes, dplyr::everything())
      p_col_for_plot <- "p_val"
      if (any(marker_diff[[p_col_for_plot]] == 0, na.rm = TRUE)) {
        min_nonzero <- min(marker_diff[[p_col_for_plot]][marker_diff[[p_col_for_plot]] > 0 & !is.na(marker_diff[[p_col_for_plot]])], na.rm = TRUE)
        marker_diff[[p_col_for_plot]][marker_diff[[p_col_for_plot]] == 0 & !is.na(marker_diff[[p_col_for_plot]])] <- if (is.finite(min_nonzero) && min_nonzero > 0) min_nonzero * 0.1 else .Machine$double.eps
      }
      sanitized_cell_type_name_for_path <- gsub("[^a-zA-Z0-9_.-]", "_", current_cell_type)
      cell_type_dir <- file.path(volcpath, sanitized_cell_type_name_for_path); if (!dir.exists(cell_type_dir)) dir.create(cell_type_dir, recursive = TRUE)
      if (save_deg_results) {
        output_deg_file <- file.path(cell_type_dir, paste0(sanitized_cell_type_name_for_path, ".marker.diff.", sanitized_ident_1, ".vs.", sanitized_ident_2,".",sanitized_assay, ".txt"))
        utils::write.table(marker_diff, file = output_deg_file, sep = "\t", row.names = FALSE, quote = FALSE)
        message(paste0("  DE results saved to: ", output_deg_file)); current_result_paths$deg_table <- output_deg_file
      } else { message(paste0("  Skipping save of DE results for '", current_cell_type, "'.")); current_result_paths$deg_table <- NULL }
      genes_for_volcano_label <- character(0)
      if (is.null(genes_to_label)) {
        sig_genes <- dplyr::filter(marker_diff, .data[[p_col_for_plot]] < p_cutoff_plot & abs(avg_log2FC) > fc_cutoff_plot)
        if (nrow(sig_genes) > 0) {
          top_up <- dplyr::slice_head(dplyr::arrange(sig_genes, dplyr::desc(avg_log2FC)), n = top_n_genes) %>% dplyr::pull(genes)
          top_down <- dplyr::slice_head(dplyr::arrange(sig_genes, avg_log2FC), n = top_n_genes) %>% dplyr::pull(genes)
          genes_for_volcano_label <- unique(c(top_up, top_down))
        }
        if (length(genes_for_volcano_label) == 0) message(paste0("  No genes for auto-labeling in '", current_cell_type, "' based on fc_cutoff_plot and p_cutoff_plot."))
      } else {
        genes_for_volcano_label <- intersect(genes_to_label, marker_diff$genes)
        if (length(genes_for_volcano_label) < length(genes_to_label)) warning(paste0("  Some provided genes_to_label not in DE results for '", current_cell_type, "'."), call. = FALSE)
      }
      if(is.null(genes_for_volcano_label)) genes_for_volcano_label <- character(0)
      
      keyvals <- ifelse(marker_diff$avg_log2FC < -fc_cutoff_plot & marker_diff[[p_col_for_plot]] < p_cutoff_plot, color_down_regulated, ifelse(marker_diff$avg_log2FC > fc_cutoff_plot & marker_diff[[p_col_for_plot]] < p_cutoff_plot, color_up_regulated, color_non_significant))
      keyvals[is.na(keyvals)] <- color_non_significant
      names(keyvals)[keyvals == color_up_regulated] <- paste0('Up (LFC > ', fc_cutoff_plot, ')')
      names(keyvals)[keyvals == color_non_significant] <- 'Non-sig'
      names(keyvals)[keyvals == color_down_regulated] <- paste0('Down (LFC < -', fc_cutoff_plot, ')')
      
      plot_title <- paste0(current_cell_type, ": ", ident_1, " vs ", ident_2, " (", assay, ")"); plot_subtitle <- paste0("LFC > |", fc_cutoff_plot, "| & P < ", p_cutoff_plot)
      
      current_x_breaks <- x_axis_breaks
      if (is.null(current_x_breaks)) {
        max_abs_lfc <- max(abs(marker_diff$avg_log2FC), na.rm = TRUE)
        if (is.finite(max_abs_lfc) && max_abs_lfc > 0) {
          limit <- ceiling(max_abs_lfc)
          if (limit <= 2) current_x_breaks <- seq(-limit, limit, by = 1)
          else if (limit <= 6) current_x_breaks <- unique(sort(c(seq(-limit, limit, by = 2), 0)))
          else current_x_breaks <- unique(sort(c(seq(-limit, limit, by = floor(limit/3)), 0)))
        } else { current_x_breaks <- ggplot2::waiver() }
      }
      
      ev_args <- list(
        toptable = marker_diff,
        lab = marker_diff$genes,
        x = 'avg_log2FC',
        y = p_col_for_plot,
        selectLab = genes_for_volcano_label,
        title = stringr::str_wrap(plot_title, 60),
        subtitle = stringr::str_wrap(plot_subtitle, 70),
        caption = paste0("LFC cut:", fc_cutoff_plot, "; P cut:", p_cutoff_plot),
        FCcutoff = fc_cutoff_plot,
        pCutoff = p_cutoff_plot,
        colCustom = keyvals,
        colAlpha = 0.8,
        legendPosition = 'right',
        legendLabSize = 10,
        legendIconSize = 4.0,
        drawConnectors = TRUE,
        widthConnectors = 0.5,
        typeConnectors = connector_arrow_type,
        arrowheads = TRUE,
        labSize = lab_font_size,
        lengthConnectors = grid::unit(connector_length, "lines"), # Convert to unit object
        min.segment.length = repel_min_segment_length,
        max.overlaps = Inf,
        pointSize = 2.0
      )
      
      volc_plot_base <- do.call(EnhancedVolcano::EnhancedVolcano, ev_args)
      volc_plot <- volc_plot_base +
        ggplot2::theme_classic() +
        ggplot2::scale_x_continuous(breaks = current_x_breaks) +
        ggplot2::labs(colour = legend_title)
      
      # --- PNG Saving with 3-Level Fallback ---
      png_filename <- file.path(cell_type_dir, paste0(sanitized_cell_type_name_for_path, "_", sanitized_ident_1, "_vs_", sanitized_ident_2,"_", sanitized_assay, ".Volcano.png"))
      png_saved_successfully <- FALSE
      if (requireNamespace("cowplot", quietly = TRUE)) {
        tryCatch({
          cowplot::save_plot(filename = png_filename, plot = volc_plot, base_width = png_width, base_height = png_height, dpi = png_dpi)
          message(paste0("  Volcano plot (PNG using cowplot::save_plot) saved to: ", png_filename))
          current_result_paths$png_plot <- png_filename
          png_saved_successfully <- TRUE
        }, error = function(e_cowplot_png) {
          warning(paste0("  cowplot::save_plot failed for PNG '", png_filename, "'. Error: ", e_cowplot_png$message, "\n  Attempting ggplot2::ggsave as fallback."), call. = FALSE)
        })
      } else { message("  cowplot package not available, skipping cowplot::save_plot for PNG. Attempting ggplot2::ggsave.") }
      if (!png_saved_successfully) {
        tryCatch({
          ggplot2::ggsave(filename = png_filename, plot = volc_plot, width = png_width, height = png_height, units = "in", dpi = png_dpi)
          message(paste0("  Volcano plot (PNG using ggplot2::ggsave) saved to: ", png_filename))
          current_result_paths$png_plot <- png_filename
          png_saved_successfully <- TRUE
        }, error = function(e_ggsave_png) {
          warning(paste0("  ggplot2::ggsave failed for PNG '", png_filename, "'. Error: ", e_ggsave_png$message, "\n  Attempting base R png() device as fallback."), call. = FALSE)
        })
      }
      if (!png_saved_successfully) {
        tryCatch({
          grDevices::png(filename = png_filename, width = png_width, height = png_height, units = "in", res = png_dpi)
          print(volc_plot)
          grDevices::dev.off()
          message(paste0("  Volcano plot (PNG using base R png()) saved to: ", png_filename))
          current_result_paths$png_plot <- png_filename
          png_saved_successfully <- TRUE
        }, error = function(e_base_png) {
          warning(paste0("  Base R png() device also failed for PNG '", png_filename, "'. Error: ", e_base_png$message, "\n  PNG plot NOT saved."), call. = FALSE)
        })
      }
      
      # --- PDF Saving with 3-Level Fallback ---
      pdf_filename <- file.path(cell_type_dir, paste0(sanitized_cell_type_name_for_path, "_", sanitized_ident_1, "_vs_", sanitized_ident_2,"_", sanitized_assay, ".Volcano.pdf"))
      pdf_saved_successfully <- FALSE
      if (requireNamespace("cowplot", quietly = TRUE)) {
        tryCatch({
          cowplot::save_plot(filename = pdf_filename, plot = volc_plot, base_width = pdf_width, base_height = pdf_height, device = "pdf")
          message(paste0("  Volcano plot (PDF using cowplot::save_plot) saved to: ", pdf_filename))
          current_result_paths$pdf_plot <- pdf_filename
          pdf_saved_successfully <- TRUE
        }, error = function(e_cowplot_pdf) {
          warning(paste0("  cowplot::save_plot failed for PDF '", pdf_filename, "'. Error: ", e_cowplot_pdf$message, "\n  Attempting ggplot2::ggsave as fallback."), call. = FALSE)
        })
      } else { message("  cowplot package not available, skipping cowplot::save_plot for PDF. Attempting ggplot2::ggsave.") }
      if (!pdf_saved_successfully) {
        tryCatch({
          ggplot2::ggsave(filename = pdf_filename, plot = volc_plot, width = pdf_width, height = pdf_height, units = "in", device = "pdf")
          message(paste0("  Volcano plot (PDF using ggplot2::ggsave) saved to: ", pdf_filename))
          current_result_paths$pdf_plot <- pdf_filename
          pdf_saved_successfully <- TRUE
        }, error = function(e_ggsave_pdf) {
          warning(paste0("  ggplot2::ggsave failed for PDF '", pdf_filename, "'. Error: ", e_ggsave_pdf$message, "\n  Attempting base R pdf() device as fallback."), call. = FALSE)
        })
      }
      if (!pdf_saved_successfully) {
        tryCatch({
          grDevices::pdf(file = pdf_filename, width = pdf_width, height = pdf_height)
          print(volc_plot)
          grDevices::dev.off()
          message(paste0("  Volcano plot (PDF using base R pdf()) saved to: ", pdf_filename))
          current_result_paths$pdf_plot <- pdf_filename
          pdf_saved_successfully <- TRUE
        }, error = function(e_base_pdf) {
          warning(paste0("  Base R pdf() device also failed for PDF '", pdf_filename, "'. Error: ", e_base_pdf$message, "\n  PDF plot NOT saved."), call. = FALSE)
        })
      }
      
      if(png_saved_successfully || pdf_saved_successfully) current_result_paths$status <- "Success (plot saved)" else current_result_paths$status <- "Success (DE done, plot saving failed)"
      
    }, error = function(e) {
      warning(paste0("An error occurred while processing cell type '", current_cell_type, "': ", e$message), call. = FALSE)
      current_result_paths$status <- paste("Error:", e$message)
    })
    all_results[[current_cell_type]] <- current_result_paths
  }
  message("\nFinished processing all specified cell types.")
  return(invisible(all_results))
}

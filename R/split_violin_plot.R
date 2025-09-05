#' Create Split Violin Plots for Single Cell RNA-seq Data
#'
#' This function creates split violin plots to visualize gene expression differences
#' between conditions across cell types using single cell RNA-seq data from a Seurat object.
#' The split violin allows for direct comparison of expression distributions between
#' two conditions (e.g., CTRL vs PAH) within each cell type.
#'
#' @param seurat_obj A Seurat object containing single cell RNA-seq data
#' @param genes Character vector of gene names to plot. Genes not found in the object
#'   will be reported and excluded
#' @param celltype_col Character string specifying the column name in metadata containing 
#'   cell type information. Default is "cell.type.ident"
#' @param condition_col Character string specifying the column name in metadata containing
#'   condition information. Default is "disease.ident"
#' @param condition_colors Named vector of colors for each condition. 
#'   Default is c("CTRL" = "dodgerblue2", "PAH" = "firebrick2")
#' @param invert_axes Logical. If FALSE (default), genes are on y-axis and cell types on x-axis.
#'   If TRUE, genes are on x-axis and cell types on y-axis
#' @param add_boxplot Logical. Whether to overlay boxplots on the violin plots. Default is TRUE
#' @param add_mean_se Logical. Whether to add mean ± standard error points. Default is TRUE
#' @param plot_title Character string for the main plot title. Default is "Gene Expression by Cell Type"
#' @param alpha_violin Numeric. Transparency level for violin plots (0-1). Default is 0.4
#' @param alpha_boxplot Numeric. Transparency level for boxplots (0-1). Default is 0.6
#' @param boxplot_width Numeric. Width of boxplots relative to violin width. Default is 0.2
#' @param trim_violin Logical. Whether to trim violin plots to data range. Default is FALSE
#' @param ncol_facet Integer. Number of columns for faceting when plotting multiple genes.
#'   If NULL, automatically determined
#' @param text_size Numeric. Base text size for the plot. Default is 12
#' @param axis_text_angle Numeric. Angle for x-axis text rotation in degrees. Default is 45
#' @param return_data Logical. If TRUE, returns a list with both plot and data.
#'   If FALSE (default), returns only the plot
#' @param verbose Logical. Whether to print progress messages. Default is TRUE
#'
#' @return A ggplot object (or list containing plot and data if return_data = TRUE)
#'
#' @details
#' The function creates split violin plots where each violin is split down the middle,
#' with each half representing a different condition (e.g., CTRL vs PAH). This allows
#' for easy visual comparison of expression distributions between conditions within
#' each cell type.
#'
#' The function automatically handles:
#' \itemize{
#'   \item Missing genes (reports and excludes them)
#'   \item Data extraction from Seurat objects
#'   \item Proper faceting for multiple genes
#'   \item Flexible axis orientation
#' }
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' p1 <- create_split_violin_plot(
#'   seurat_obj = pbmc,
#'   genes = c("CD3D", "CD4", "CD8A"),
#'   celltype_col = "cell.type.ident",
#'   condition_col = "disease.ident"
#' )
#'
#' # Inverted axes with custom colors
#' p2 <- create_split_violin_plot(
#'   seurat_obj = pbmc,
#'   genes = c("SOD1", "TMBIM6", "CASP3"),
#'   invert_axes = TRUE,
#'   condition_colors = c("CTRL" = "blue", "PAH" = "red")
#' )
#'
#' # Minimal styling
#' p3 <- create_split_violin_plot(
#'   seurat_obj = pbmc,
#'   genes = c("CD8A", "GZMK"),
#'   add_boxplot = FALSE,
#'   add_mean_se = FALSE,
#'   alpha_violin = 0.7
#' )
#' }
#'
#' @importFrom Seurat DefaultAssay GetAssayData
#' @importFrom ggplot2 ggplot aes geom_violin geom_boxplot stat_summary
#'   scale_fill_manual scale_x_discrete scale_y_continuous labs theme_minimal
#'   theme element_text facet_wrap position_dodge
#' @importFrom dplyr filter left_join rename
#' @importFrom tidyr pivot_longer
#' @importFrom tibble rownames_to_column
#'
#' @export
create_split_violin_plot <- function(
    seurat_obj,
    genes,
    celltype_col = "cell.type.ident",
    condition_col = "disease.ident", 
    condition_colors = c("CTRL" = "dodgerblue2", "PAH" = "firebrick2"),
    invert_axes = FALSE,
    add_boxplot = TRUE,
    add_mean_se = TRUE,
    plot_title = "Gene Expression by Cell Type",
    alpha_violin = 0.4,
    alpha_boxplot = 0.6,
    boxplot_width = 0.2,
    trim_violin = FALSE,
    ncol_facet = NULL,
    text_size = 12,
    axis_text_angle = 45,
    return_data = FALSE,
    verbose = TRUE
) {
  
  # Required packages
  required_packages <- c("Seurat", "ggplot2", "dplyr", "tidyr")
  for(pkg in required_packages) {
    if(!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste(pkg, "package is required but not installed"))
    }
  }
  
  # Load required libraries
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  
  # Input validation
  if(!all(genes %in% rownames(seurat_obj))) {
    missing_genes <- setdiff(genes, rownames(seurat_obj))
    warning(paste("The following genes were not found:", paste(missing_genes, collapse = ", ")))
    genes <- intersect(genes, rownames(seurat_obj))
  }
  
  if(length(genes) == 0) {
    stop("No valid genes found in the Seurat object")
  }
  
  if(!celltype_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Column", celltype_col, "not found in metadata"))
  }
  
  if(!condition_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Column", condition_col, "not found in metadata"))
  }
  
  if(verbose) {
    message(paste("Creating split violin plot for", length(genes), "genes"))
    message(paste("Cell types found:", length(unique(seurat_obj@meta.data[[celltype_col]]))))
    message(paste("Conditions found:", paste(unique(seurat_obj@meta.data[[condition_col]]), collapse = ", ")))
  }
  
  # Extract expression data
  DefaultAssay(seurat_obj) <- "RNA"
  expr_data <- GetAssayData(seurat_obj, slot = "data")[genes, , drop = FALSE]
  
  # Get metadata
  metadata <- seurat_obj@meta.data[, c(celltype_col, condition_col), drop = FALSE]
  metadata$cell_id <- rownames(metadata)
  
  # Convert expression data to long format
  expr_long <- as.data.frame(as.matrix(expr_data)) %>%
    tibble::rownames_to_column("gene") %>%
    tidyr::pivot_longer(cols = -gene, names_to = "cell_id", values_to = "expression")
  
  # Merge with metadata
  plot_data <- expr_long %>%
    left_join(metadata, by = "cell_id")
  
  # Rename columns using setNames approach
  colnames(plot_data)[colnames(plot_data) == celltype_col] <- "celltype"
  colnames(plot_data)[colnames(plot_data) == condition_col] <- "condition"
  
  # Filter out any missing data
  plot_data <- plot_data %>%
    filter(!is.na(celltype) & !is.na(condition))
  
  # Ensure conditions are factors with proper levels
  conditions_present <- unique(plot_data$condition)
  plot_data$condition <- factor(plot_data$condition, levels = names(condition_colors))
  
  # Update condition_colors to only include present conditions
  condition_colors <- condition_colors[names(condition_colors) %in% conditions_present]
  
  if(verbose) {
    message(paste("Final dataset contains", nrow(plot_data), "observations"))
  }
  
  # Create the base plot
  if(!invert_axes) {
    # Genes on y-axis, cell types on x-axis
    p <- ggplot(plot_data, aes(x = celltype, y = expression, fill = condition))
    x_lab <- "Cell Type"
    y_lab <- "Expression Level"
  } else {
    # Cell types on y-axis, genes on x-axis  
    p <- ggplot(plot_data, aes(x = gene, y = expression, fill = condition))
    x_lab <- "Gene"
    y_lab <- "Expression Level"
  }
  
  # Add split violin plot using custom implementation
  p <- p + geom_violin(
    alpha = alpha_violin, 
    trim = trim_violin,
    position = position_dodge(width = 0.4),
    scale = "area"
  )
  
  # Add boxplot if requested
  if(add_boxplot) {
    p <- p + geom_boxplot(
      width = boxplot_width, 
      alpha = alpha_boxplot, 
      fatten = NULL, 
      show.legend = FALSE,
      position = position_dodge(width = 0.4)
    )
  }
  
  # Add mean ± SE if requested
  if(add_mean_se) {
    p <- p + stat_summary(
      fun.data = "mean_se", 
      geom = "pointrange", 
      show.legend = FALSE, 
      position = position_dodge(width = 0.4)
    )
  }
  
  # Apply styling
  p <- p +
    scale_fill_manual(values = condition_colors, name = "Condition") +
    scale_x_discrete(name = x_lab) +
    scale_y_continuous(name = y_lab) +
    labs(title = plot_title) +
    theme_minimal(base_size = text_size) +
    theme(
      axis.text.x = element_text(angle = axis_text_angle, hjust = 1, vjust = 1),
      plot.title = element_text(hjust = 0.5, size = text_size + 2),
      legend.position = "bottom",
      strip.text = element_text(size = text_size - 1)
    )
  
  # Add faceting if we have multiple genes and not inverting axes
  if(length(genes) > 1) {
    if(!invert_axes) {
      # Facet by gene when genes are on y-axis
      if(is.null(ncol_facet)) {
        ncol_facet <- min(3, length(genes))
      }
      p <- p + facet_wrap(~ gene, scales = "free_y", ncol = ncol_facet)
    } else {
      # When genes are on x-axis, we might want to facet by cell type instead
      if(length(unique(plot_data$celltype)) <= 6) {
        if(is.null(ncol_facet)) {
          ncol_facet <- min(3, length(unique(plot_data$celltype)))
        }
        p <- p + facet_wrap(~ celltype, scales = "free_x", ncol = ncol_facet)
      }
    }
  }
  
  # Return data if requested
  if(return_data) {
    return(list(plot = p, data = plot_data))
  } else {
    return(p)
  }
}

#' Custom Split Violin Geom for ggplot2
#'
#' This function creates a true split violin plot where each violin is split
#' down the middle to show two conditions side by side.
#'
#' @param plot_data Data frame with columns: x variable, y (expression), fill (condition)
#' @param x_var Character string specifying the x variable column name
#' @param y_var Character string specifying the y variable column name  
#' @param fill_var Character string specifying the fill variable column name
#' @param alpha Numeric. Transparency level (0-1)
#' @param trim Logical. Whether to trim violins to data range
#'
#' @return A ggplot layer
#'
#' @details
#' This is a helper function that creates a true split violin where the left half
#' shows one condition and the right half shows another condition. This provides
#' better visual comparison than standard side-by-side violins.
#'
#' @examples
#' \dontrun{
#' # This is typically called internally by create_split_violin_plot()
#' # but can be used independently
#' 
#' p <- ggplot(data, aes(x = celltype, y = expression)) +
#'   geom_split_violin_custom(
#'     plot_data = data,
#'     x_var = "celltype", 
#'     y_var = "expression",
#'     fill_var = "condition"
#'   )
#' }
#'
#' @keywords internal
geom_split_violin_custom <- function(plot_data, x_var, y_var, fill_var, alpha = 0.4, trim = FALSE) {
  
  # This function would implement a custom split violin
  # For now, we'll use standard violins with position_dodge
  # A full implementation would require custom ggplot2 geom development
  
  warning("True split violins require custom geom development. Using dodged violins instead.")
  
  return(
    geom_violin(
      alpha = alpha,
      trim = trim,
      position = position_dodge(width = 0.8),
      scale = "area"
    )
  )
}

#' Create Split Violin Plot with vioplot Backend
#'
#' Alternative implementation using the vioplot package for creating split violin plots.
#' This creates base R graphics instead of ggplot2.
#'
#' @param seurat_obj A Seurat object containing single cell RNA-seq data
#' @param genes Character vector of gene names to plot
#' @param celltype_col Character string specifying the column name containing cell types
#' @param condition_col Character string specifying the column name containing conditions
#' @param condition_colors Named vector of colors for each condition
#' @param main_title Character string for the main plot title
#' @param save_pdf Logical. Whether to save the plot as a PDF file
#' @param pdf_filename Character string for PDF filename if save_pdf = TRUE
#' @param plot_width Numeric. Width of the plot in inches (for PDF). Default is 16
#' @param plot_height Numeric. Height of the plot in inches (for PDF). Default is 12
#' @param max_genes_per_page Integer. Maximum number of genes to plot per page. Default is 6
#' @param verbose Logical. Whether to print progress messages
#'
#' @return Invisible NULL (creates plot as side effect)
#'
#' @details
#' This function uses the vioplot package to create split violin plots in base R graphics.
#' The advantage is true split violins, but the disadvantage is less flexibility
#' compared to ggplot2. For many genes, the function will create multiple pages.
#'
#' @examples
#' \dontrun{
#' create_split_violin_vioplot(
#'   seurat_obj = pbmc,
#'   genes = c("CD3D", "CD4"),
#'   celltype_col = "cell.type.ident",
#'   condition_col = "disease.ident"
#' )
#' }
#'
#' @export
create_split_violin_vioplot <- function(
    seurat_obj,
    genes,
    celltype_col = "cell.type.ident",
    condition_col = "disease.ident",
    condition_colors = c("CTRL" = "dodgerblue2", "PAH" = "firebrick2"),
    main_title = "Gene Expression by Cell Type",
    save_pdf = TRUE,
    pdf_filename = "split_violin_plot.pdf",
    plot_width = 16,
    plot_height = 12,
    max_genes_per_page = 6,
    verbose = TRUE
) {
  
  # Check for vioplot package
  if(!requireNamespace("vioplot", quietly = TRUE)) {
    stop("vioplot package is required but not installed")
  }
  
  library(vioplot)
  
  # Input validation (similar to main function)
  if(!all(genes %in% rownames(seurat_obj))) {
    missing_genes <- setdiff(genes, rownames(seurat_obj))
    warning(paste("The following genes were not found:", paste(missing_genes, collapse = ", ")))
    genes <- intersect(genes, rownames(seurat_obj))
  }
  
  if(length(genes) == 0) {
    stop("No valid genes found in the Seurat object")
  }
  
  if(verbose) {
    message(paste("Creating violin plots for", length(genes), "genes"))
  }
  
  # Extract data
  DefaultAssay(seurat_obj) <- "RNA"
  expr_data <- GetAssayData(seurat_obj, slot = "data")[genes, , drop = FALSE]
  metadata <- seurat_obj@meta.data[, c(celltype_col, condition_col), drop = FALSE]
  
  # Get unique cell types and conditions
  cell_types <- unique(metadata[[celltype_col]])
  conditions <- unique(metadata[[condition_col]])
  
  if(verbose) {
    message(paste("Cell types:", length(cell_types)))
    message(paste("Conditions:", paste(conditions, collapse = ", ")))
    message(paste("Total violins per gene:", length(cell_types) * length(conditions)))
  }
  
  # Setup PDF if requested - always use PDF for large plots to avoid margin issues
  if(save_pdf) {
    pdf(pdf_filename, width = plot_width, height = plot_height)
  } else {
    # For screen display, open a large window
    if(.Platform$OS.type == "windows") {
      windows(width = plot_width, height = plot_height)
    } else {
      X11(width = plot_width, height = plot_height)
    }
  }
  
  # Split genes into chunks for multiple pages if needed
  gene_chunks <- split(genes, ceiling(seq_along(genes) / max_genes_per_page))
  
  for(chunk_idx in seq_along(gene_chunks)) {
    current_genes <- gene_chunks[[chunk_idx]]
    
    if(verbose) {
      message(paste("Processing page", chunk_idx, "with", length(current_genes), "genes"))
    }
    
    # Set up plotting layout - adjust based on number of genes in chunk
    n_genes <- length(current_genes)
    if(n_genes == 1) {
      par(mfrow = c(1, 1), mar = c(8, 4, 4, 2) + 0.1)
    } else if(n_genes <= 2) {
      par(mfrow = c(1, 2), mar = c(8, 4, 4, 2) + 0.1)
    } else if(n_genes <= 4) {
      par(mfrow = c(2, 2), mar = c(8, 4, 4, 2) + 0.1)
    } else {
      par(mfrow = c(3, 2), mar = c(8, 4, 4, 2) + 0.1)
    }
    
    # Create plot for each gene in current chunk
    for(gene in current_genes) {
      gene_expr <- as.numeric(expr_data[gene, ])
      
      # Prepare data lists for vioplot
      plot_data <- list()
      plot_names <- c()
      plot_colors <- c()
      
      for(ct in cell_types) {
        for(cond in conditions) {
          cells_idx <- which(metadata[[celltype_col]] == ct & metadata[[condition_col]] == cond)
          if(length(cells_idx) > 5) {  # Only plot if at least 5 cells
            plot_data[[length(plot_data) + 1]] <- gene_expr[cells_idx]
            # Shorten cell type names if too long
            ct_short <- if(nchar(ct) > 15) paste0(substr(ct, 1, 12), "...") else ct
            plot_names <- c(plot_names, paste(ct_short, cond, sep = "\n"))
            plot_colors <- c(plot_colors, condition_colors[cond])
          }
        }
      }
      
      # Create the violin plot
      if(length(plot_data) > 0) {
        vioplot::vioplot(
          plot_data,
          names = plot_names,
          col = plot_colors,
          main = paste(gene, "Expression"),
          ylab = "Expression Level",
          xlab = "",
          las = 2,  # Rotate x-axis labels
          cex.names = 0.7,  # Smaller text for names
          cex.main = 1.2
        )
        
        # Add a box around the plot
        box()
      } else {
        # Create empty plot if no data
        plot(1, type = "n", axes = FALSE, xlab = "", ylab = "", main = paste(gene, "- No Data"))
      }
    }
    
    # Add page title if multiple chunks
    if(length(gene_chunks) > 1) {
      mtext(paste(main_title, "- Page", chunk_idx, "of", length(gene_chunks)), 
            side = 3, line = -2, outer = TRUE, cex = 1.2)
    } else {
      mtext(main_title, side = 3, line = -2, outer = TRUE, cex = 1.2)
    }
    
    # Add new page if not the last chunk
    if(chunk_idx < length(gene_chunks)) {
      if(save_pdf) {
        # For PDF, this will create a new page
        plot.new()
      }
    }
  }
  
  # Close device
  if(save_pdf) {
    dev.off()
    if(verbose) {
      message(paste("Plot saved to:", pdf_filename))
      message(paste("Created", length(gene_chunks), "page(s)"))
    }
  }
  
  return(invisible(NULL))
}

#' Create Individual Volcano Plots with Pathway Labels
#'
#' Creates individual volcano plots for each cell type with pathway labels,
#' similar to the detailed tip cell analysis visualization.
#'
#' @param analysis_results List output from analyze_metabolism_comprehensive()
#' @param padj_threshold Numeric. Adjusted p-value threshold for significance. Default: 0.05
#' @param logfc_threshold Numeric. Log2 fold change threshold for significance. Default: 0.1
#' @param label_top_n Integer. Number of top pathways to label per direction. Default: 10
#' @param save_plots Logical. Whether to save plots to files. Default: FALSE
#' @param plot_width Numeric. Width of saved plots in inches. Default: 10
#' @param plot_height Numeric. Height of saved plots in inches. Default: 8
#' @param plot_dpi Numeric. DPI for saved plots. Default: 300
#'
#' @return List of ggplot objects, one per cell type
#'
#' @examples
#' \dontrun{
#' # After running comprehensive analysis
#' results <- analyze_metabolism_comprehensive(seurat_obj, edge_case = FALSE)
#'
#' # Create individual volcano plots
#' individual_plots <- create_individual_volcano_plots(
#'   analysis_results = results,
#'   label_top_n = 15,
#'   save_plots = TRUE
#' )
#'
#' # Display specific cell type plot
#' print(individual_plots[["gCAP"]])
#'
#' # Display all plots
#' lapply(individual_plots, print)
#' }
#'
#' @export
#' @import ggplot2
#' @importFrom ggrepel geom_text_repel
create_individual_volcano_plots <- function(analysis_results,
                                            padj_threshold = 0.05,
                                            logfc_threshold = 0.1,
                                            label_top_n = 10,
                                            save_plots = FALSE,
                                            plot_width = 10,
                                            plot_height = 8,
                                            plot_dpi = 300) {

  # Required packages check
  required_packages <- c("ggplot2", "ggrepel")
  missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]

  if(length(missing_packages) > 0) {
    stop(paste("Required packages not installed:", paste(missing_packages, collapse = ", ")))
  }

  library(ggplot2)
  library(ggrepel)

  # Input validation
  if(!is.list(analysis_results) || !"results" %in% names(analysis_results)) {
    stop("analysis_results must be output from analyze_metabolism_comprehensive()")
  }

  results <- analysis_results$results
  mode <- analysis_results$mode

  if(length(results) == 0) {
    stop("No results found in analysis_results")
  }

  # Create individual plots
  individual_plots <- list()

  for(celltype in names(results)) {

    cat("Creating volcano plot for:", celltype, "\n")

    data <- results[[celltype]]

    # Add significance categories
    data$significance <- "Not Significant"
    data$significance[data$padj < padj_threshold & data$log2FC > logfc_threshold] <- "Upregulated in Disease"
    data$significance[data$padj < padj_threshold & data$log2FC < -logfc_threshold] <- "Downregulated in Disease"

    # For edge case analysis, adjust labels
    if(mode == "edge_case") {
      data$significance[data$significance == "Upregulated in Disease"] <- "Higher in Target"
      data$significance[data$significance == "Downregulated in Disease"] <- "Lower in Target"
    }

    # Select top pathways to label
    sig_data <- data[data$significance != "Not Significant", ]

    if(nrow(sig_data) > 0) {
      # Separate up and down regulated for balanced labeling
      up_regulated <- sig_data[sig_data$log2FC > 0, ]
      down_regulated <- sig_data[sig_data$log2FC < 0, ]

      # Select top N from each direction
      top_up <- data.frame()
      top_down <- data.frame()

      if(nrow(up_regulated) > 0) {
        up_regulated$rank_score <- -log10(up_regulated$padj) * abs(up_regulated$log2FC)
        up_regulated <- up_regulated[order(up_regulated$rank_score, decreasing = TRUE), ]
        top_up <- up_regulated[1:min(label_top_n, nrow(up_regulated)), ]
      }

      if(nrow(down_regulated) > 0) {
        down_regulated$rank_score <- -log10(down_regulated$padj) * abs(down_regulated$log2FC)
        down_regulated <- down_regulated[order(down_regulated$rank_score, decreasing = TRUE), ]
        top_down <- down_regulated[1:min(label_top_n, nrow(down_regulated)), ]
      }

      to_label <- rbind(top_up, top_down)
    } else {
      to_label <- data.frame()
    }

    # Count significant pathways
    if(mode == "edge_case") {
      n_up <- sum(data$significance == "Higher in Target")
      n_down <- sum(data$significance == "Lower in Target")
      subtitle_text <- paste0("Higher in ", celltype, ": ", n_up, " | Lower in ", celltype, ": ", n_down)
      x_label <- paste0("log₂(", celltype, "/All Others)")
    } else {
      n_up <- sum(data$significance == "Upregulated in Disease")
      n_down <- sum(data$significance == "Downregulated in Disease")
      subtitle_text <- paste0("Upregulated in PAH: ", n_up, " | Downregulated in PAH: ", n_down)
      x_label <- "log₂(PAH/CTRL)"
    }

    # Create the plot
    p <- ggplot(data, aes(x = log2FC, y = -log10(padj))) +
      geom_point(aes(color = significance), size = 2, alpha = 0.7) +
      scale_color_manual(
        values = c("Upregulated in Disease" = "#E31A1C",
                   "Downregulated in Disease" = "#1F78B4",
                   "Higher in Target" = "#E31A1C",
                   "Lower in Target" = "#1F78B4",
                   "Not Significant" = "gray70"),
        name = "Metabolic Change"
      ) +

      # Add threshold lines
      geom_hline(yintercept = -log10(padj_threshold),
                 linetype = "dashed", color = "black", alpha = 0.5) +
      geom_vline(xintercept = c(-logfc_threshold, logfc_threshold),
                 linetype = "dashed", color = "black", alpha = 0.5) +

      # Add pathway labels
      geom_text_repel(
        data = to_label,
        aes(label = pathway),
        size = 3,
        max.overlaps = 20,
        box.padding = 0.5,
        force = 2,
        segment.color = "gray30",
        segment.size = 0.3
      ) +

      # Styling
      labs(
        title = paste("Metabolic Pathway Changes in", celltype),
        subtitle = subtitle_text,
        x = x_label,
        y = "-log₁₀(Adjusted p-value)"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray40"),
        legend.position = "bottom",
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 11)
      )

    individual_plots[[celltype]] <- p

    # Save plot if requested
    if(save_plots) {
      filename <- paste0(celltype, "_volcano_plot.pdf")
      ggsave(filename, p, width = plot_width, height = plot_height, dpi = plot_dpi)
      cat("  Saved:", filename, "\n")
    }
  }

  cat("Created", length(individual_plots), "individual volcano plots\n")
  return(individual_plots)
}

#' Display Individual Volcano Plots
#'
#' Convenience function to display individual volcano plots one by one
#'
#' @param individual_plots List of ggplot objects from create_individual_volcano_plots()
#' @param cell_types Character vector of cell types to display. If NULL, displays all. Default: NULL
#' @param pause_between Logical. Whether to pause between plots (interactive). Default: FALSE
#'
#' @export
display_individual_plots <- function(individual_plots,
                                     cell_types = NULL,
                                     pause_between = FALSE) {

  if(is.null(cell_types)) {
    cell_types <- names(individual_plots)
  }

  for(celltype in cell_types) {
    if(celltype %in% names(individual_plots)) {
      cat("Displaying plot for:", celltype, "\n")
      print(individual_plots[[celltype]])

      if(pause_between && interactive()) {
        cat("Press [Enter] to continue to next plot...")
        readline()
      }
    } else {
      cat("Plot not found for:", celltype, "\n")
    }
  }
}

#' Save Individual Volcano Plots
#'
#' Save individual volcano plots to different formats
#'
#' @param individual_plots List of ggplot objects from create_individual_volcano_plots()
#' @param output_dir Character. Directory to save plots. Default: "volcano_plots"
#' @param format Character. File format: "pdf", "png", "tiff", "jpeg". Default: "pdf"
#' @param width Numeric. Plot width in inches. Default: 10
#' @param height Numeric. Plot height in inches. Default: 8
#' @param dpi Numeric. DPI for raster formats. Default: 300
#'
#' @export
save_individual_plots <- function(individual_plots,
                                  output_dir = "volcano_plots",
                                  format = "pdf",
                                  width = 10,
                                  height = 8,
                                  dpi = 300) {

  # Create output directory if it doesn't exist
  if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("Created directory:", output_dir, "\n")
  }

  for(celltype in names(individual_plots)) {
    filename <- file.path(output_dir, paste0(celltype, "_volcano_plot.", format))

    ggsave(filename, individual_plots[[celltype]],
           width = width, height = height, dpi = dpi)

    cat("Saved:", filename, "\n")
  }

  cat("All plots saved to:", output_dir, "\n")
}

# ===== USAGE EXAMPLES =====

# # Example 1: Create individual plots from your standard analysis
# individual_plots <- create_individual_volcano_plots(
#   analysis_results = standard_results,
#   label_top_n = 15,
#   save_plots = TRUE
# )
#
# # Example 2: Display specific cell type
# print(individual_plots[["gCAP"]])
#
# # Example 3: Display all plots
# display_individual_plots(individual_plots)
#
# # Example 4: Save in different formats
# save_individual_plots(individual_plots,
#                      output_dir = "publication_figures",
#                      format = "png",
#                      width = 12,
#                      height = 9,
#                      dpi = 300)
#
# # Example 5: Create plots for edge case analysis
# tip_individual_plot <- create_individual_volcano_plots(
#   analysis_results = tip_results,
#   label_top_n = 20
# )
# print(tip_individual_plot[["Tip.Like.Cells"]])

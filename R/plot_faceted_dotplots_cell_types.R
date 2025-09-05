#' Create Faceted Dot Plot for Cell Type Distribution by Condition
#'
#' @description
#' This function creates a faceted dot plot that visualizes cell type distribution across 
#' multiple conditions (e.g., disease states). The plot represents three dimensions of data:
#' cell type, absolute count (dot size), and relative proportion (dot color).
#'
#' @param seurat_obj A Seurat object containing the single-cell data
#' @param cell_type_column Character string specifying the column in the meta.data that 
#'   contains cell type information. Default is "cell.type".
#' @param condition_column Character string specifying the column in the meta.data that 
#'   contains condition information (e.g., disease status). Default is "disease".
#' @param text_angle Numeric value specifying the angle of x-axis text labels. Default is 45.
#' @param max_dot_size Numeric value specifying the maximum size of dots. Default is 15.
#' @param ctrl_condition Character string specifying the value in condition_column that represents
#'   control samples. Default is "CTRL".
#' @param disease_condition Character string specifying the value in condition_column that represents
#'   disease samples. Default is "PAH".
#' @param ctrl_color Character string specifying the color for control condition. Default is "dodgerblue".
#' @param disease_color Character string specifying the color for disease condition. Default is "firebrick2".
#' @param sort_by Character string specifying how to sort cell types. Options are "count" (default) 
#'   or "name". If "count", cell types are sorted by total cell count.
#' @param show_labels Logical indicating whether to show count and percentage labels. Default is TRUE.
#' @param title Character string specifying the plot title. Default is "Cell Type Distribution by Condition".
#' @param subtitle Character string specifying the plot subtitle. Default is NULL.
#' @param x_lab Character string specifying the x-axis label. Default is "Cell Type".
#' @param y_lab Character string specifying the y-axis label. Default is the value of condition_column.
#'
#' @return A ggplot object representing the faceted dot plot
#'
#' @details
#' The plot shows each cell type on the x-axis and condition on the y-axis. For each 
#' cell type/condition combination, a dot is drawn where:
#' - The size of the dot represents the absolute number of cells
#' - The color of the dot represents the condition (CTRL = dodgerblue, PAH = firebrick2)
#' - Text labels show both the count and percentage values
#'
#' This visualization is particularly useful for identifying shifts in cell type 
#' distribution between conditions (e.g., disease vs. control).
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' plot_faceted_dotplots_cell_types(seurat_obj)
#'
#' # Customized plot
#' plot_faceted_dotplots_cell_types(
#'   seurat_obj,
#'   cell_type_column = "CellType",
#'   condition_column = "Treatment",
#'   ctrl_condition = "Control",
#'   disease_condition = "Treated",
#'   max_dot_size = 20,
#'   title = "Cell Type Changes After Treatment"
#' )
#'
#' # For only stromal cells
#' stromal_cells <- subset(seurat_obj, lineage == "Stromal")
#' plot_faceted_dotplots_cell_types(stromal_cells)
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_point scale_size_area scale_color_manual 
#'    geom_text labs theme_minimal theme element_text element_rect
#' @importFrom dplyr group_by summarise mutate ungroup arrange
#' @importFrom rlang sym
#'
#' @export
plot_faceted_dotplots_cell_types <- function(seurat_obj, 
                                             cell_type_column = "cell.type",
                                             condition_column = "disease",
                                             text_angle = 45,
                                             max_dot_size = 15,
                                             ctrl_condition = "CTRL",
                                             disease_condition = "PAH",
                                             ctrl_color = "dodgerblue",
                                             disease_color = "firebrick2",
                                             sort_by = "count", # "count" or "name"
                                             show_labels = TRUE,
                                             title = "Cell Type Distribution by Condition",
                                             subtitle = NULL,
                                             x_lab = "Cell Type",
                                             y_lab = NULL) {
  
  # Required packages
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is needed for this function to work. Please install it.",
         call. = FALSE)
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  # Set y-axis label to condition_column if not specified
  if (is.null(y_lab)) {
    y_lab <- condition_column
  }
  
  # Extract metadata
  meta_data <- seurat_obj@meta.data
  
  # Verify columns exist in metadata
  if (!cell_type_column %in% colnames(meta_data)) {
    stop(paste("Column", cell_type_column, "not found in Seurat object metadata."))
  }
  if (!condition_column %in% colnames(meta_data)) {
    stop(paste("Column", condition_column, "not found in Seurat object metadata."))
  }
  
  # Count cells by cell type and condition
  count_data <- meta_data %>%
    dplyr::group_by(!!rlang::sym(cell_type_column), !!rlang::sym(condition_column)) %>%
    dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
    dplyr::ungroup()
  
  # Calculate percentages and totals
  count_data <- count_data %>%
    dplyr::group_by(!!rlang::sym(cell_type_column)) %>%
    dplyr::mutate(total = sum(count),
                  percentage = count / total * 100) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(scaled_count = count / max(count) * 100)
  
  # Sort by total count or name
  if (sort_by == "count") {
    cell_type_totals <- count_data %>%
      dplyr::group_by(!!rlang::sym(cell_type_column)) %>%
      dplyr::summarise(total = sum(count), .groups = "drop") %>%
      dplyr::arrange(dplyr::desc(total))
    
    count_data[[cell_type_column]] <- factor(count_data[[cell_type_column]], 
                                             levels = cell_type_totals[[cell_type_column]])
  } else if (sort_by == "name") {
    count_data[[cell_type_column]] <- factor(count_data[[cell_type_column]], 
                                             levels = sort(unique(count_data[[cell_type_column]])))
  }
  
  # Create custom color mapping
  color_values <- c()
  color_values[ctrl_condition] <- ctrl_color
  color_values[disease_condition] <- disease_color
  
  # Create the plot
  p <- ggplot2::ggplot(count_data, 
                       ggplot2::aes(x = !!rlang::sym(cell_type_column), 
                                    y = !!rlang::sym(condition_column))) +
    ggplot2::geom_point(ggplot2::aes(size = count, color = !!rlang::sym(condition_column))) +
    ggplot2::scale_size_area(max_size = max_dot_size) +
    ggplot2::scale_color_manual(values = color_values,
                                name = "Condition")
  
  # Add text labels if requested
  if (show_labels) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = paste0(count, "\n(", round(percentage), "%)")),
      vjust = 1.5, 
      size = 3
    )
  }
  
  # Add labels and theme
  p <- p + ggplot2::labs(
    title = title,
    subtitle = subtitle,
    x = x_lab,
    y = y_lab,
    size = "Cell Count"
  ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = text_angle, hjust = 1),
      legend.position = "right",
      panel.border = ggplot2::element_rect(color = "lightgray", fill = NA),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )
  
  return(p)
}

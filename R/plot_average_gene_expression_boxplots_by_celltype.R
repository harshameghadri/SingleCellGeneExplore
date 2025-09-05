#' Plot Average Gene Expression Boxplots or Raincloud Plots by Cell Type
#'
#' @description
#' Generates boxplots or raincloud plots comparing the average expression of specified genes
#' across different cell types. Data is expected in a list of tibbles, where each
#' list element corresponds to a cell type. Plots are grouped and colored by
#' a 'disease' condition. Individual data points can be overlaid as jittered points,
#' colored by subject. Statistical significance can be displayed as stars.
#' Plots can be generated as individual plots or as faceted plots with flexible
#' Y-axis scaling. Plots are saved to PDF and/or PNG files and can be displayed
#' in the R graphics device.
#'
#' @param data_list A named list of tibbles. The names of the list elements are
#'   assumed to represent cell types (any ".APS" suffix will be removed). Each
#'   tibble within the list must contain at least the following columns:
#'   \itemize{
#'     \item `gene`: Character vector of gene names.
#'     \item `average_expression`: Numeric vector of expression values.
#'     \item `disease`: Character or factor vector indicating the condition
#'       (e.g., "CTRL", "PAH"). Used for boxplot fill color.
#'     \item `subject`: Character or factor vector indicating the subject ID.
#'       Used for jitter point color when `jitter = TRUE`.
#'   }
#' @param genes_of_interest A character vector specifying the genes for which
#'   to generate plots. Plots will only be created for genes found in the data.
#' @param subject_colors A named character vector where names are subject IDs
#'   (matching the 'subject' column in `data_list` tibbles) and values are
#'   valid R color strings (e.g., `c("S1" = "red", "S2" = "blue")`). Used to
#'   color jitter points when `jitter = TRUE`.
#' @param filter_column A character string specifying the column name to use for
#'   filtering cell types. Default is "lineage.type". Can be any column present
#'   in the tibbles (e.g., "lineage.type.2", "cell_category", etc.).
#' @param filter_values A character vector specifying which values from the
#'   `filter_column` to include in plots. If NULL (default), all cell types
#'   are included. Examples: c("Immune"), c("Myeloid", "Lymphoid"), 
#'   c("CD4.T.Cells", "B.Cells").
#' @param significance_data A named list of tibbles containing statistical results.
#'   Names should match cell type names (with or without ".APS" suffix). Each
#'   tibble should contain columns: `gene` and the column specified in 
#'   `significance_column`. If NULL (default), no significance stars are added.
#' @param significance_column A character string specifying which p-value column
#'   to use for significance testing. Default is "adjusted_p_value.BH".
#'   Alternative: "p_value".
#' @param significance_thresholds A named numeric vector specifying p-value
#'   thresholds for significance stars. Default is 
#'   `c("***" = 0.001, "**" = 0.01, "*" = 0.05)`.
#' @param plot_title_prefix A character string used as the prefix for the
#'   output PDF and PNG filenames. Default is "Gene_Expression_by_Cell_Type".
#' @param jitter Logical. If `TRUE` (default), adds individual data points to the
#'   plots, colored by subject according to `subject_colors`.
#' @param plot_raincloud Logical. If `TRUE`, creates raincloud plots instead of
#'   traditional boxplots. Raincloud plots show violin plots (distribution), 
#'   box plots (quartiles), and individual data points. Default is `FALSE`.
#' @param path Optional. A character string specifying the directory path where
#'   the output plot files should be saved. If `NULL` (default), plots are
#'   saved in the current working directory.
#' @param save_pdf Logical. If `TRUE` (default), saves plots as a PDF file.
#' @param save_png Logical. If `TRUE`, saves plots as a PNG file. Default is `FALSE`.
#' @param view_plot Logical or numeric. If `TRUE`, displays all plots. If `FALSE`,
#'   no plots are displayed. If numeric, displays that many plot pages.
#'   Default is `TRUE`.
#' @param plots_per_page Integer. Number of plots to display per page/panel.
#'   Default is 2. Can be set to any number (e.g., 6 for 6 plots per page).
#' @param x_axis_text_size Numeric. Font size for x-axis text. Default is 10.
#' @param y_axis_text_size Numeric. Font size for y-axis text. Default is 10.
#' @param x_axis_angle Numeric. Angle for x-axis text rotation in degrees.
#'   Default is 45.
#' @param y_axis_limits A numeric vector of length 2 specifying custom y-axis
#'   limits c(min, max). If NULL (default), limits are calculated automatically.
#'   When `faceted_free_scales = TRUE`, this parameter's behavior depends on
#'   `respect_manual_limits`.
#' @param adaptive_y_limits Logical. If `TRUE` (default), each gene gets its own
#'   optimized Y-axis limits. If `FALSE`, uses global limits across all genes.
#'   Ignored when `y_axis_limits` is specified or `faceted_free_scales = TRUE`.
#' @param outlier_method Character. Method for handling outliers when calculating
#'   adaptive Y-limits. Options: "iqr" (default), "percentile", "none".
#' @param outlier_factor Numeric. Factor for IQR-based outlier detection.
#'   Default is 1.5. Larger values are more permissive of outliers.
#' @param percentile_range Numeric vector of length 2. Percentile range to use
#'   when `outlier_method = "percentile"`. Default is c(0.05, 0.95).
#' @param faceted_free_scales Logical. If `TRUE`, creates faceted plots instead
#'   of individual plots. Allows for flexible Y-axis scaling per gene.
#'   Default is `FALSE`.
#' @param facet_scaling Character. Controls facet scaling behavior when
#'   `faceted_free_scales = TRUE`. Options: "auto" (default), "free_y", 
#'   "free_x", "fixed". "auto" intelligently chooses based on expression range
#'   diversity.
#' @param respect_manual_limits Logical. When `TRUE` (default) and 
#'   `y_axis_limits` is specified, manual limits override facet free scaling.
#'   When `FALSE`, facet free scaling takes precedence over manual limits.
#' @param disease_colors A named character vector for disease condition colors.
#'   Default is `c("CTRL" = "dodgerblue", "PAH" = "firebrick")`. Names should
#'   match the values in the `disease` column of your data.
#' @param verbose Logical. If `TRUE`, prints detailed progress messages through
#'   all processing stages. Default is `FALSE`.
#'
#' @return A list containing:
#'   \itemize{
#'     \item `plots`: List of ggplot objects created
#'     \item `genes_plotted`: Character vector of genes successfully plotted
#'     \item `genes_missing`: Character vector of genes not found in data
#'     \item `cell_types_used`: Character vector of cell types included in plots
#'     \item `plot_info`: List with plotting parameters used
#'   }
#'
#' @importFrom dplyr bind_rows mutate filter %>% distinct pull left_join group_by summarise
#' @importFrom stringr str_remove
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_jitter geom_violin geom_point
#'   position_jitterdodge position_dodge theme_classic labs theme element_text
#'   scale_fill_manual scale_color_manual coord_cartesian annotate ylim facet_wrap
#'   element_rect
#' @importFrom patchwork wrap_plots
#' @importFrom grDevices pdf png dev.off
#' @importFrom purrr map map_dfr
#' @importFrom stats quantile sd
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic usage with default settings
#' result <- Plot_Average_Gene_Expression_BoxPlots_byCellType(
#'   data_list = my_data_list,
#'   genes_of_interest = c("TP53", "BRCA1"),
#'   subject_colors = my_subject_colors
#' )
#'
#' # Filter by immune lineage with significance stars
#' result <- Plot_Average_Gene_Expression_BoxPlots_byCellType(
#'   data_list = my_data_list,
#'   genes_of_interest = c("TP53", "BRCA1", "CD4", "CD8A"),
#'   subject_colors = my_subject_colors,
#'   filter_column = "lineage.type",
#'   filter_values = "Immune",
#'   significance_data = my_stats_data,
#'   verbose = TRUE
#' )
#'
#' # Raincloud plots with faceted free scales for genes with different ranges
#' result <- Plot_Average_Gene_Expression_BoxPlots_byCellType(
#'   data_list = my_data_list,
#'   genes_of_interest = c("TEK", "FOXO1", "FOXO3", "FOXF1"),
#'   subject_colors = my_subject_colors,
#'   plot_raincloud = TRUE,
#'   faceted_free_scales = TRUE,
#'   facet_scaling = "free_y",
#'   plots_per_page = 4,
#'   save_png = TRUE
#' )
#'
#' # Custom Y-axis limits with adaptive scaling
#' result <- Plot_Average_Gene_Expression_BoxPlots_byCellType(
#'   data_list = my_data_list,
#'   genes_of_interest = c("GENE1", "GENE2"),
#'   subject_colors = my_subject_colors,
#'   adaptive_y_limits = TRUE,
#'   outlier_method = "percentile",
#'   percentile_range = c(0.1, 0.9)
#' )
#' }
Plot_Average_Gene_Expression_BoxPlots_byCellType.2 <- function(
    data_list,
    genes_of_interest,
    subject_colors,
    filter_column = "lineage.type",
    filter_values = NULL,
    significance_data = NULL,
    significance_column = "adjusted_p_value.BH",
    significance_thresholds = c("***" = 0.001, "**" = 0.01, "*" = 0.05),
    plot_title_prefix = "Gene_Expression_by_Cell_Type",
    jitter = TRUE,
    plot_raincloud = FALSE,
    path = NULL,
    save_pdf = TRUE,
    save_png = FALSE,
    view_plot = TRUE,
    plots_per_page = 2,
    x_axis_text_size = 10,
    y_axis_text_size = 10,
    x_axis_angle = 45,
    y_axis_limits = NULL,
    adaptive_y_limits = TRUE,
    outlier_method = "iqr",
    outlier_factor = 1.5,
    percentile_range = c(0.05, 0.95),
    faceted_free_scales = FALSE,
    facet_scaling = "auto",
    respect_manual_limits = TRUE,
    disease_colors = c("CTRL" = "dodgerblue", "PAH" = "firebrick"),
    verbose = FALSE) {
  
  # Load required libraries
  if (verbose) cat("Stage 1: Loading required libraries...\n")
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(stringr)
  library(purrr)
  
  # Initialize return object
  result_object <- list(
    plots = list(),
    genes_plotted = character(),
    genes_missing = character(),
    cell_types_used = character(),
    plot_info = list()
  )
  
  # Stage 2: Data preparation
  if (verbose) cat("Stage 2: Preparing and combining data...\n")
  
  # Combine all tibbles into a single dataframe
  data <- bind_rows(data_list, .id = "cell_type") %>%
    mutate(cell_type_clean = str_remove(cell_type, "\\.APS.*$"))
  
  if (verbose) {
    cat("  - Combined", length(data_list), "cell type datasets\n")
    cat("  - Total rows in combined data:", nrow(data), "\n")
  }
  
  # Stage 3: Apply filters
  if (verbose) cat("Stage 3: Applying filters...\n")
  
  if (!is.null(filter_values)) {
    if (!filter_column %in% colnames(data)) {
      stop("Filter column '", filter_column, "' not found in data. Available columns: ", 
           paste(colnames(data), collapse = ", "))
    }
    
    original_rows <- nrow(data)
    data <- data %>% filter(.data[[filter_column]] %in% filter_values)
    
    if (nrow(data) == 0) {
      stop("No data found for ", filter_column, " values: ", paste(filter_values, collapse = ", "))
    }
    
    if (verbose) {
      cat("  - Filtered by", filter_column, "for values:", paste(filter_values, collapse = ", "), "\n")
      cat("  - Rows after filtering:", nrow(data), "(removed", original_rows - nrow(data), "rows)\n")
      cat("  - Cell types included:", paste(unique(data$cell_type_clean), collapse = ", "), "\n")
    }
  }
  
  result_object$cell_types_used <- unique(data$cell_type_clean)
  
  # Stage 4: Filter genes
  if (verbose) cat("Stage 4: Filtering genes of interest...\n")
  
  data_to_plot <- data %>%
    filter(gene %in% genes_of_interest)
  
  # Get genes actually present in data
  genes_present <- data_to_plot %>%
    distinct(gene) %>%
    pull(gene)
  
  # Maintain order from genes_of_interest
  genes_present <- intersect(genes_of_interest, genes_present)
  genes_not_present <- setdiff(genes_of_interest, genes_present)
  
  result_object$genes_plotted <- genes_present
  result_object$genes_missing <- genes_not_present
  
  if (verbose) {
    cat("  - Genes requested:", length(genes_of_interest), "\n")
    cat("  - Genes found in data:", length(genes_present), "\n")
    if (length(genes_not_present) > 0) {
      cat("  - Genes missing:", paste(genes_not_present, collapse = ", "), "\n")
    }
  }
  
  if (length(genes_present) == 0) {
    if (verbose) cat("Stage 4 ERROR: No genes found in data to plot.\n")
    return(result_object)
  }
  
  # Stage 5: Prepare significance data
  if (verbose && !is.null(significance_data)) {
    cat("Stage 5: Preparing significance annotations...\n")
  }
  
  # Function to get significance stars
  get_significance_stars <- function(p_value, thresholds = significance_thresholds) {
    if (is.na(p_value)) return("")
    
    for (i in 1:length(thresholds)) {
      if (p_value <= thresholds[i]) {
        return(names(thresholds)[i])
      }
    }
    return("")
  }
  
  # Function to calculate smart Y-axis limits
  calculate_smart_ylimits <- function(expression_data, method = "iqr", factor = 1.5, percentiles = c(0.05, 0.95)) {
    if (method == "iqr") {
      Q1 <- quantile(expression_data, 0.25, na.rm = TRUE)
      Q3 <- quantile(expression_data, 0.75, na.rm = TRUE)
      IQR <- Q3 - Q1
      
      y_min <- max(0, Q1 - factor * IQR)
      y_max <- Q3 + factor * IQR
      
    } else if (method == "percentile") {
      y_min <- quantile(expression_data, percentiles[1], na.rm = TRUE)
      y_max <- quantile(expression_data, percentiles[2], na.rm = TRUE)
      
    } else {
      y_min <- min(expression_data, na.rm = TRUE)
      y_max <- max(expression_data, na.rm = TRUE)
    }
    
    # Add small buffer
    buffer <- (y_max - y_min) * 0.1
    return(c(y_min - buffer, y_max + buffer))
  }
  
  significance_annotations <- NULL
  if (!is.null(significance_data)) {
    if (!significance_column %in% colnames(significance_data[[1]])) {
      warning("Significance column '", significance_column, "' not found in significance data. ",
              "Available columns: ", paste(colnames(significance_data[[1]]), collapse = ", "))
    } else {
      sig_data_combined <- map_dfr(names(significance_data), function(cell_name) {
        cell_clean <- str_remove(cell_name, "\\.APS.*$")
        
        if ("gene" %in% colnames(significance_data[[cell_name]]) && 
            significance_column %in% colnames(significance_data[[cell_name]])) {
          
          significance_data[[cell_name]] %>%
            dplyr::select(gene, !!significance_column) %>%
            mutate(
              cell_type_clean = cell_clean,
              significance_stars = map_chr(.data[[significance_column]], get_significance_stars)
            ) %>%
            filter(significance_stars != "")
        } else {
          NULL
        }
      })
      
      significance_annotations <- sig_data_combined
      
      if (verbose && !is.null(significance_annotations)) {
        cat("  - Significance annotations prepared for", 
            length(unique(significance_annotations$cell_type_clean)), "cell types\n")
        cat("  - Total significant associations:", nrow(significance_annotations), "\n")
      }
    }
  }
  
  # Stage 6: Organize plots
  if (verbose) cat("Stage 6: Organizing plots into groups...\n")
  
  n_genes <- length(genes_present)
  genes_to_group <- genes_present
  remainder <- n_genes %% plots_per_page
  
  if (remainder != 0) {
    genes_to_group <- c(genes_present, rep(NA, plots_per_page - remainder))
  }
  
  if (length(genes_to_group) > 0) {
    gene_groups <- matrix(genes_to_group, ncol = plots_per_page, byrow = TRUE)
  } else {
    gene_groups <- matrix(nrow = 0, ncol = plots_per_page)
  }
  
  if (verbose) {
    cat("  - Total plot groups:", nrow(gene_groups), "\n")
    cat("  - Plots per page:", plots_per_page, "\n")
    cat("  - Faceted free scales:", faceted_free_scales, "\n")
  }
  
  # Stage 7: Setup output paths
  if (verbose) cat("Stage 7: Setting up output paths...\n")
  
  if (is.null(path)) {
    path_dir <- getwd()
  } else {
    path_dir <- path
    if (!dir.exists(path_dir)) {
      if (verbose) cat("  - Creating directory:", path_dir, "\n")
      dir.create(path_dir, recursive = TRUE)
    }
  }
  
  safe_prefix <- gsub("[^a-zA-Z0-9_.-]", "_", plot_title_prefix)
  if (!is.null(filter_values)) {
    safe_prefix <- paste0(safe_prefix, "_", gsub("[^a-zA-Z0-9_.-]", "_", paste(filter_values, collapse = "_")))
  }
  if (plot_raincloud) {
    safe_prefix <- paste0(safe_prefix, "_raincloud")
  }
  if (faceted_free_scales) {
    safe_prefix <- paste0(safe_prefix, "_faceted")
  }
  
  pdf_file_name <- file.path(path_dir, paste0(safe_prefix, ".pdf"))
  png_file_name <- file.path(path_dir, paste0(safe_prefix, ".png"))
  
  if (verbose) {
    if (save_pdf) cat("  - PDF output:", pdf_file_name, "\n")
    if (save_png) cat("  - PNG output:", png_file_name, "\n")
  }
  
  # Store plot info
  result_object$plot_info <- list(
    plot_type = if(plot_raincloud) "raincloud" else "boxplot",
    faceted = faceted_free_scales,
    facet_scaling = if(faceted_free_scales) facet_scaling else "none",
    plots_per_page = plots_per_page,
    adaptive_y_limits = adaptive_y_limits,
    manual_y_limits = y_axis_limits
  )
  
  # Stage 8: Generate plots
  if (verbose) cat("Stage 8: Generating plots...\n")
  
  all_plots <- list()
  
  if (nrow(gene_groups) > 0) {
    for (i in 1:nrow(gene_groups)) {
      if (verbose) cat("  - Creating plot group", i, "of", nrow(gene_groups), "\n")
      
      # Get non-NA genes for this group
      group_genes <- gene_groups[i, !is.na(gene_groups[i, ])]
      
      if (faceted_free_scales && length(group_genes) > 1) {
        # Create faceted plot
        if (verbose) cat("    - Creating faceted plot for genes:", paste(group_genes, collapse = ", "), "\n")
        
        group_data <- data_to_plot %>%
          filter(gene %in% group_genes) %>%
          filter(!is.na(average_expression))
        
        if (nrow(group_data) > 0) {
          # Determine facet scaling behavior
          if (!is.null(y_axis_limits) && respect_manual_limits) {
            final_scales <- "free_x"
            use_manual_limits <- TRUE
            if (verbose) cat("      - Using manual Y limits with free X scales\n")
            
          } else if (facet_scaling == "auto") {
            # Smart auto-detection based on expression range diversity
            gene_ranges <- group_data %>%
              group_by(gene) %>%
              summarise(range_size = max(average_expression, na.rm = TRUE) - 
                          min(average_expression, na.rm = TRUE), .groups = "drop")
            
            range_cv <- sd(gene_ranges$range_size) / mean(gene_ranges$range_size)
            final_scales <- if (range_cv > 1.0) "free_y" else "free_x"
            use_manual_limits <- FALSE
            if (verbose) cat("      - Auto-detected scaling:", final_scales, "(CV =", round(range_cv, 2), ")\n")
            
          } else {
            final_scales <- facet_scaling
            use_manual_limits <- !is.null(y_axis_limits) && facet_scaling %in% c("free_x", "fixed")
            if (verbose) cat("      - Using specified scaling:", final_scales, "\n")
          }
          
          # Create base plot
          if (plot_raincloud) {
            p <- ggplot(group_data, aes(x = cell_type_clean, y = average_expression)) +
              geom_violin(aes(fill = disease), 
                          alpha = 0.5,
                          position = position_dodge(width = 0.9),
                          trim = FALSE, scale = "width") +
              geom_boxplot(aes(fill = disease),
                           alpha = 0.7, width = 0.2,
                           position = position_dodge(width = 0.9),
                           outlier.shape = NA) +
              {if(jitter) geom_point(aes(color = subject),
                                     position = position_jitterdodge(jitter.width = 0.1, 
                                                                     dodge.width = 0.9,
                                                                     jitter.height = 0),
                                     size = 1.5, alpha = 0.8)}
          } else {
            p <- ggplot(group_data, aes(x = cell_type_clean, y = average_expression)) +
              geom_boxplot(aes(fill = disease), alpha = 0.5, outlier.shape = NA,
                           position = position_dodge(width = 0.8)) +
              {if(jitter) geom_jitter(aes(color = subject),
                                      position = position_jitterdodge(jitter.width = 0.2, 
                                                                      dodge.width = 0.8),
                                      size = 1, alpha = 0.7)}
          }
          
          # Add faceting
          p <- p + facet_wrap(~ gene, scales = final_scales, 
                              ncol = min(plots_per_page, length(group_genes)))
          
          # Apply styling
          p <- p +
            theme_classic() +
            labs(title = paste("Gene Expression:", paste(group_genes, collapse = ", ")),
                 x = "Cell Type",
                 y = "Average Expression") +
            theme(
              axis.text.x = element_text(angle = x_axis_angle, hjust = 1, size = x_axis_text_size),
              axis.text.y = element_text(size = y_axis_text_size),
              axis.title = element_text(size = 12),
              plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
              legend.position = "right",
              strip.text = element_text(size = 12, face = "bold"),
              strip.background = element_rect(fill = "lightgray", alpha = 0.5)
            ) +
            scale_fill_manual(name = "Condition", values = disease_colors, drop = FALSE) +
            scale_color_manual(name = "Subject", values = subject_colors, na.value = "grey50")
          
          # Apply Y-axis limits if needed
          if (use_manual_limits) {
            p <- p + coord_cartesian(ylim = y_axis_limits)
          }
          
          # Add significance stars
          if (!is.null(significance_annotations)) {
            for (gene in group_genes) {
              gene_sig <- significance_annotations %>% filter(gene == !!gene)
              
              if (nrow(gene_sig) > 0) {
                for (k in 1:nrow(gene_sig)) {
                  cell_type_name <- gene_sig$cell_type_clean[k]
                  stars <- gene_sig$significance_stars[k]
                  
                  gene_cell_data <- group_data %>% 
                    filter(gene == !!gene, cell_type_clean == cell_type_name)
                  
                  if (nrow(gene_cell_data) > 0) {
                    if (use_manual_limits) {
                      star_y <- y_axis_limits[2] * 0.95
                    } else {
                      star_y <- max(gene_cell_data$average_expression, na.rm = TRUE) * 1.1
                    }
                    
                    p <- p + annotate("text", 
                                      x = cell_type_name, 
                                      y = star_y,
                                      label = stars, 
                                      size = 3.5, 
                                      vjust = 0)
                  }
                }
              }
            }
          }
          
          all_plots[[length(all_plots) + 1]] <- p
        }
        
      } else {
        # Create individual plots
        plot_list <- list()
        
        for (j in 1:plots_per_page) {
          gene <- gene_groups[i, j]
          if (!is.na(gene)) {
            if (verbose) cat("    - Processing gene:", gene, "\n")
            
            gene_data <- data_to_plot %>%
              filter(gene == !!gene) %>%
              filter(!is.na(average_expression))
            
            if (nrow(gene_data) > 0) {
              # Calculate Y-axis limits
              if (is.null(y_axis_limits)) {
                if (adaptive_y_limits) {
                  y_limits <- calculate_smart_ylimits(gene_data$average_expression, 
                                                      method = outlier_method, 
                                                      factor = outlier_factor,
                                                      percentiles = percentile_range)
                  if (verbose) cat("      - Adaptive Y limits:", round(y_limits, 3), "\n")
                } else {
                  # Global limits
                  all_expression <- data_to_plot$average_expression[!is.na(data_to_plot$average_expression)]
                  y_limits <- calculate_smart_ylimits(all_expression, method = "none")
                }
              } else {
                y_limits <- y_axis_limits
                if (verbose) cat("      - Using manual Y limits:", y_limits, "\n")
              }
              
              # Create base plot
              if (plot_raincloud) {
                # Check if gghalves is available, otherwise use custom flat violin
                if (requireNamespace("gghalves", quietly = TRUE)) {
                  library(gghalves)
                  
                  p <- ggplot(gene_data, aes(x = cell_type_clean, y = average_expression, fill = disease)) +
                    # Half violin (cloud) on the left
                    gghalves::geom_half_violin(
                      aes(fill = disease),
                      position = position_nudge(x = -0.3),
                      side = "l",  # left side only
                      alpha = 0.7,
                      trim = FALSE
                    ) +
                    # Box plot in the center-right
                    geom_boxplot(
                      aes(fill = disease),
                      position = position_nudge(x = -0.1),
                      width = 0.15,
                      alpha = 0.8,
                      outlier.shape = NA
                    ) +
                    # Points (rain) on the right
                    {if(jitter) geom_point(
                      aes(color = subject),
                      position = position_jitter(width = 0.1, height = 0),
                      size = 1.5,
                      alpha = 0.8
                    )}
                  
                } else {
                  # Fallback: Create custom flat violin or use alternative
                  warning("gghalves package not available. Using alternative raincloud implementation.")
                  
                  p <- ggplot(gene_data, aes(x = cell_type_clean, y = average_expression)) +
                    # Create half violin effect using stat_density and coord manipulation
                    stat_density(
                      aes(fill = disease),
                      geom = "polygon",
                      position = position_nudge(x = -0.3),
                      alpha = 0.7,
                      trim = FALSE
                    ) +
                    geom_boxplot(
                      aes(fill = disease),
                      position = position_nudge(x = -0.1),
                      width = 0.15,
                      alpha = 0.8,
                      outlier.shape = NA
                    ) +
                    {if(jitter) geom_point(
                      aes(color = subject),
                      position = position_jitter(width = 0.1, height = 0),
                      size = 1.5,
                      alpha = 0.8
                    )}
                }
              }
              
              # Add styling
              p <- p +
                theme_classic() +
                labs(title = gene,
                     x = "Cell Type",
                     y = "Average Expression") +
                theme(
                  axis.text.x = element_text(angle = x_axis_angle, hjust = 1, size = x_axis_text_size),
                  axis.text.y = element_text(size = y_axis_text_size),
                  axis.title = element_text(size = 12),
                  plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
                  legend.position = "right"
                ) +
                scale_fill_manual(name = "Condition",
                                  values = disease_colors,
                                  drop = FALSE) +
                scale_color_manual(name = "Subject",
                                   values = subject_colors,
                                   na.value = "grey50")
              
              # Apply Y-axis limits
              if (is.null(y_axis_limits)) {
                p <- p + coord_cartesian(ylim = y_limits)
              } else {
                p <- p + ylim(y_axis_limits)
              }
              
              # Add significance stars
              if (!is.null(significance_annotations)) {
                gene_sig <- significance_annotations %>%
                  filter(gene == !!gene)
                
                if (nrow(gene_sig) > 0) {
                  for (k in 1:nrow(gene_sig)) {
                    cell_type_name <- gene_sig$cell_type_clean[k]
                    stars <- gene_sig$significance_stars[k]
                    
                    cell_data <- gene_data %>% filter(cell_type_clean == cell_type_name)
                    if (nrow(cell_data) > 0) {
                      star_y <- if (is.null(y_axis_limits)) {
                        max(cell_data$average_expression, na.rm = TRUE) + (y_limits[2] - y_limits[1]) * 0.05
                      } else {
                        y_axis_limits[2] * 0.95
                      }
                      
                      p <- p + annotate("text", 
                                        x = cell_type_name, 
                                        y = star_y,
                                        label = stars, 
                                        size = 4, 
                                        vjust = 0)
                    }
                  }
                }
              }
              
              plot_list[[length(plot_list) + 1]] <- p
            } else {
              if (verbose) cat("      - No data found for gene:", gene, "\n")
            }
          }
        }
        
        # Combine individual plots in group
        if (length(plot_list) > 0) {
          if (plots_per_page == 1) {
            combined_plot <- plot_list[[1]]
          } else {
            combined_plot <- wrap_plots(plot_list, ncol = min(plots_per_page, length(plot_list)))
          }
          all_plots[[length(all_plots) + 1]] <- combined_plot
        }
      }
    }
  }
  
  result_object$plots <- all_plots
  
  if (verbose) cat("  - Generated", length(all_plots), "plot groups\n")
  
  # Stage 9: Save files
  if (verbose) cat("Stage 9: Saving output files...\n")
  
  if (length(all_plots) > 0) {
    # Calculate dimensions based on plot type and faceting
    if (faceted_free_scales) {
      pdf_width <- min(20, 4 * plots_per_page)
      pdf_height <- if(plot_raincloud) 12 else 10
    } else {
      pdf_width <- min(20, 5 * plots_per_page)
      pdf_height <- if(plot_raincloud) 10 else 8
    }
    
    # Save PDF
    if (save_pdf) {
      tryCatch({
        pdf(pdf_file_name, width = pdf_width, height = pdf_height)
        for(plot_obj in all_plots) {
          print(plot_obj)
        }
        dev.off()
        
        if (verbose) cat("  - PDF saved successfully:", basename(pdf_file_name), "\n")
      }, error = function(e) {
        warning("Failed to save PDF: ", e$message)
        if(names(dev.cur()) == "pdf") dev.off()
      })
    }
    
    # Save PNG
    if (save_png) {
      tryCatch({
        png(png_file_name, width = pdf_width, height = pdf_height, units = "in", res = 300)
        for(plot_obj in all_plots) {
          print(plot_obj)
        }
        dev.off()
        
        if (verbose) cat("  - PNG saved successfully:", basename(png_file_name), "\n")
      }, error = function(e) {
        warning("Failed to save PNG: ", e$message)
        if(names(dev.cur()) == "png") dev.off()
      })
    }
  }
  
  # Stage 10: Display plots
  if (verbose) cat("Stage 10: Displaying plots...\n")
  
  if (length(all_plots) > 0) {
    if (is.logical(view_plot) && view_plot) {
      if (interactive()) {
        if (verbose) cat("  - Displaying all", length(all_plots), "plot groups\n")
        for(plot_obj in all_plots) {
          print(plot_obj)
        }
      } else {
        if (verbose) cat("  - Non-interactive session: plots saved but not displayed\n")
      }
    } else if (is.numeric(view_plot) && view_plot > 0) {
      plots_to_show <- min(view_plot, length(all_plots))
      if (interactive()) {
        if (verbose) cat("  - Displaying first", plots_to_show, "plot groups\n")
        for(i in 1:plots_to_show) {
          print(all_plots[[i]])
        }
      } else {
        if (verbose) cat("  - Non-interactive session: plots saved but not displayed\n")
      }
    } else {
      if (verbose) cat("  - Plot display disabled by user\n")
    }
  }
  
  # Stage 11: Final summary
  if (verbose) {
    cat("Stage 11: Process complete!\n")
    cat("Summary:\n")
    cat("  - Genes plotted:", length(result_object$genes_plotted), "\n")
    cat("  - Cell types used:", length(result_object$cell_types_used), "\n")
    cat("  - Plot groups created:", length(result_object$plots), "\n")
    cat("  - Plot type:", result_object$plot_info$plot_type, "\n")
    cat("  - Faceted plots:", result_object$plot_info$faceted, "\n")
    if (save_pdf) cat("  - PDF file:", basename(pdf_file_name), "\n")
    if (save_png) cat("  - PNG file:", basename(png_file_name), "\n")
  }
  
  return(invisible(result_object))
}

# Helper function to create example data for testing
#' Create Example Data for Plot Function Testing
#'
#' @description
#' Creates sample data in the format expected by Plot_Average_Gene_Expression_BoxPlots_byCellType
#' for testing and demonstration purposes.
#'
#' @param n_subjects_ctrl Integer. Number of control subjects. Default is 5.
#' @param n_subjects_pah Integer. Number of PAH subjects. Default is 5.
#' @param genes Character vector. Gene names to include. Default includes common genes.
#' @param cell_types Character vector. Cell type names to include.
#' @param include_lineage Logical. Whether to include lineage.type columns. Default is TRUE.
#'
#' @return A list containing:
#'   \itemize{
#'     \item `data_list`: Named list of tibbles ready for plotting function
#'     \item `subject_colors`: Named vector of subject colors
#'     \item `significance_data`: Sample significance data (optional)
#'   }
#'
#' @export
#'
#' @examples
#' # Create example data
#' example_data <- create_example_gene_expression_data()
#' 
#' # Use with plotting function
#' result <- Plot_Average_Gene_Expression_BoxPlots_byCellType(
#'   data_list = example_data$data_list,
#'   genes_of_interest = c("TP53", "BRCA1"),
#'   subject_colors = example_data$subject_colors,
#'   verbose = TRUE
#' )
create_example_gene_expression_data <- function(
    n_subjects_ctrl = 5,
    n_subjects_pah = 5,
    genes = c("TP53", "BRCA1", "EGFR", "TEK", "FOXO1", "FOXO3", "CD4", "CD8A"),
    cell_types = c("CD4.T.Cells", "CD8.T.Cells", "B.Cells", "NK.Cells", "Alveolar.Macrophages"),
    include_lineage = TRUE) {
  
  library(dplyr)
  library(tidyr)
  
  # Create subject IDs
  subjects_ctrl <- paste0("CTRL_", sprintf("%02d", 1:n_subjects_ctrl))
  subjects_pah <- paste0("PAH_", sprintf("%02d", 1:n_subjects_pah))
  all_subjects <- c(subjects_ctrl, subjects_pah)
  
  # Create subject colors
  n_subjects <- length(all_subjects)
  subject_colors <- rainbow(n_subjects)
  names(subject_colors) <- all_subjects
  
  # Function to create data for one cell type
  create_cell_type_data <- function(cell_type, subjects_ctrl, subjects_pah, genes) {
    # Create base expression patterns with some biological realism
    base_expression <- case_when(
      grepl("TEK", genes) ~ rnorm(length(genes), mean = 0.1, sd = 0.05),
      grepl("FOXO", genes) ~ rnorm(length(genes), mean = 1.5, sd = 0.3),
      grepl("CD", genes) ~ rnorm(length(genes), mean = 0.8, sd = 0.2),
      TRUE ~ rnorm(length(genes), mean = 1.0, sd = 0.25)
    )
    
    # Add disease effect
    disease_effect <- case_when(
      grepl("TP53|BRCA1", genes) ~ 0.5,  # Upregulated in PAH
      grepl("TEK", genes) ~ -0.02,       # Slightly downregulated
      TRUE ~ 0.1
    )
    
    # Create data
    ctrl_data <- expand_grid(subject = subjects_ctrl, gene = genes) %>%
      mutate(
        disease = "CTRL",
        average_expression = pmax(0, rep(base_expression, length(subjects_ctrl)) + 
                                    rnorm(n(), mean = 0, sd = 0.1))
      )
    
    pah_data <- expand_grid(subject = subjects_pah, gene = genes) %>%
      mutate(
        disease = "PAH",
        average_expression = pmax(0, rep(base_expression + disease_effect, length(subjects_pah)) + 
                                    rnorm(n(), mean = 0, sd = 0.1))
      )
    
    combined_data <- bind_rows(ctrl_data, pah_data)
    
    # Add lineage information if requested
    if (include_lineage) {
      lineage_type <- case_when(
        grepl("T.Cells|B.Cells|NK", cell_type) ~ "Immune",
        grepl("Macrophages", cell_type) ~ "Immune",
        TRUE ~ "Other"
      )
      
      lineage_type_2 <- case_when(
        grepl("T.Cells|B.Cells|NK", cell_type) ~ "Lymphoid",
        grepl("Macrophages", cell_type) ~ "Myeloid",
        TRUE ~ lineage_type
      )
      
      combined_data <- combined_data %>%
        mutate(
          lineage.type = lineage_type,
          lineage.type.2 = lineage_type_2
        )
    }
    
    return(combined_data)
  }
  
  # Create data list
  data_list <- map(cell_types, ~create_cell_type_data(.x, subjects_ctrl, subjects_pah, genes)) %>%
    setNames(paste0(cell_types, ".APS"))
  
  # Create example significance data
  significance_data <- map(cell_types, function(cell_type) {
    tibble(
      gene = genes,
      p_value = runif(length(genes), 0, 0.1),
      adjusted_p_value.BH = p.adjust(p_value, method = "BH")
    )
  }) %>%
    setNames(paste0(cell_types, ".APS"))
  
  return(list(
    data_list = data_list,
    subject_colors = subject_colors,
    significance_data = significance_data
  ))
}

# Example usage and testing function
#' Test Plot Function with Various Configurations
#'
#' @description
#' Runs the plotting function with different parameter combinations to demonstrate
#' functionality and help with testing.
#'
#' @param save_outputs Logical. Whether to save plot outputs. Default is FALSE.
#' @param output_dir Character. Directory for saving test outputs.
#'
#' @export
test_plot_function <- function(save_outputs = FALSE, output_dir = tempdir()) {
  
  cat("Creating example data...\n")
  example_data <- create_example_gene_expression_data()
  
  test_genes <- c("TP53", "TEK", "FOXO1", "CD4")
  
  # Test 1: Basic boxplots
  cat("\nTest 1: Basic boxplots\n")
  result1 <- Plot_Average_Gene_Expression_BoxPlots_byCellType(
    data_list = example_data$data_list,
    genes_of_interest = test_genes,
    subject_colors = example_data$subject_colors,
    plot_title_prefix = "Test1_Basic",
    save_pdf = save_outputs,
    save_png = FALSE,
    view_plot = 1,
    path = if(save_outputs) output_dir else NULL,
    verbose = TRUE
  )
  
  # Test 2: Raincloud plots with significance
  cat("\nTest 2: Raincloud plots with significance\n")
  result2 <- Plot_Average_Gene_Expression_BoxPlots_byCellType(
    data_list = example_data$data_list,
    genes_of_interest = test_genes,
    subject_colors = example_data$subject_colors,
    significance_data = example_data$significance_data,
    plot_raincloud = TRUE,
    plot_title_prefix = "Test2_Raincloud",
    save_pdf = save_outputs,
    view_plot = 1,
    path = if(save_outputs) output_dir else NULL,
    verbose = TRUE
  )
  
  # Test 3: Faceted free scales
  cat("\nTest 3: Faceted plots with free Y scales\n")
  result3 <- Plot_Average_Gene_Expression_BoxPlots_byCellType(
    data_list = example_data$data_list,
    genes_of_interest = test_genes,
    subject_colors = example_data$subject_colors,
    faceted_free_scales = TRUE,
    facet_scaling = "free_y",
    plots_per_page = 4,
    plot_title_prefix = "Test3_Faceted",
    save_pdf = save_outputs,
    view_plot = 1,
    path = if(save_outputs) output_dir else NULL,
    verbose = TRUE
  )
  
  # Test 4: Lineage filtering
  cat("\nTest 4: Lineage filtering\n")
  result4 <- Plot_Average_Gene_Expression_BoxPlots_byCellType(
    data_list = example_data$data_list,
    genes_of_interest = test_genes,
    subject_colors = example_data$subject_colors,
    filter_column = "lineage.type",
    filter_values = "Immune",
    adaptive_y_limits = TRUE,
    plot_title_prefix = "Test4_Immune_Only",
    save_pdf = save_outputs,
    view_plot = 1,
    path = if(save_outputs) output_dir else NULL,
    verbose = TRUE
  )
  
  cat("\nAll tests completed!\n")
  if (save_outputs) {
    cat("Outputs saved to:", output_dir, "\n")
    cat("Files created:", paste(list.files(output_dir, pattern = "Test.*\\.pdf"), collapse = ", "), "\n")
  }
  
  return(list(
    test1 = result1,
    test2 = result2,
    test3 = result3,
    test4 = result4
  ))
}
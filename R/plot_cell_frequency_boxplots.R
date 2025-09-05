#' Create Boxplots with Significance Stars for Cell Type Frequencies
#'
#' @description
#' This function generates boxplots showing the frequency distribution of cell types between
#' conditions (e.g., disease vs control), with significance stars indicating statistical
#' differences based on Wilcoxon rank-sum tests (unpaired) or Wilcoxon signed-rank tests (paired).
#' Relies on `wilcox.test` internal handling for small sample sizes.
#'
#' @param cell_type_column Character string specifying the column in meta.data containing cell types.
#' @param condition_column Character string specifying the column in meta.data containing disease/control status.
#' @param control_color Character string specifying color for control (default: "dodgerblue2").
#' @param control_level Character string specifying the control level in condition_column (default: "CTRL").
#' @param custom_x_axis_order Optional character vector. If provided, it dictates the specific order of
#'   cell types on the x-axis, overriding `sort_alphabetically`.
#'   Cell types from the data that are included in this vector will be ordered as specified.
#'   Cell types present in the data but NOT in `custom_x_axis_order` will be appended to the end,
#'   sorted alphabetically (A-Z).
#'   Cell types specified in `custom_x_axis_order` but NOT present in the data will be ignored,
#'   and a warning will be issued. (Default: `NULL`).
#' @param disease_color Character string specifying color for disease (default: "firebrick").
#' @param disease_level Character string specifying the disease level in condition_column (default: "PAH").
#' @param dpi Numeric value specifying resolution for PNG output (default: 300).
#' @param height Numeric value specifying figure height in inches (default: 6).
#' @param output_dir Character string specifying directory to save output files (default: current directory).
#' @param output_prefix Character string for the output file names (without extension).
#' @param paired Logical indicating whether samples should be treated as paired (default: FALSE). Requires samples to be present in both conditions.
#' @param plot_title Character string specifying the plot title (default: "Cell Type Frequencies by Condition").
#' @param sample_column Character string specifying the column in meta.data containing sample IDs. Used for pairing if `paired = TRUE`.
#' @param save_results Logical indicating whether to save the significance results to a file (default: TRUE).
#' @param seurat_obj A Seurat object containing the single-cell data.
#' @param sort_alphabetically Logical or character string defining the sorting order of cell types on the
#'   x-axis if `custom_x_axis_order` is `NULL`.
#'   Accepted values:
#'   \itemize{
#'     \item `TRUE`: Sorts cell types alphabetically (A-Z).
#'     \item `FALSE` (default): Also sorts cell types alphabetically (A-Z), maintaining consistency with previous default behavior.
#'     \item `"INVERSE"`: Sorts cell types in reverse alphabetical order (Z-A).
#'   }
#'   If an invalid value is provided, a warning is issued, and it defaults to `FALSE` (A-Z sorting).
#' @param star_height Numeric value specifying the absolute y-coordinate for placing significance stars (default: 95).
#'   Ensure this value is appropriate for the y-axis scale, especially when `y_max_limit` is set.
#' @param width Numeric value specifying figure width in inches (default: 10).
#' @param x_angle Numeric value specifying the angle of x-axis labels (default: 45).
#' @param y_max_limit Optional numeric value specifying the maximum limit for the y-axis.
#'   If NULL (default), the y-axis limit is determined dynamically based on data and `star_height`.
#'   If set, `star_height` should ideally be less than `y_max_limit`.
#'
#' @return A ggplot object of the boxplot with significance stars.
#'
#' @details
#' This function calculates the frequency of each cell type within each sample.
#' It then performs Wilcoxon tests (unpaired rank-sum or paired signed-rank)
#' to compare conditions for each cell type (provided the type exists in both conditions)
#' and adds significance stars (* p<0.05, ** p<0.01, *** p<0.001)
#' or "ns" (not significant) to the plot. It saves the plot as both PDF and PNG files,
#' and optionally saves the test results. For paired tests, only samples present under
#' both conditions are included in the test for a given cell type.
#' The order of cell types on the x-axis can be controlled using `sort_alphabetically`
#' or `custom_x_axis_order`.
#' Note: This version removed explicit minimum sample size checks before running Wilcoxon tests,
#' relying on `stats::wilcox.test`'s internal handling.
#'
#' @examples
#' \dontrun{
#' # Assuming PAH.obj.Int.Stromal.4.clean is a Seurat object
#' # and necessary packages (dplyr, ggplot2, rlang, tidyr if paired) are installed.
#'
#' # Example with default y-axis scaling and default A-Z sorting
#' plot_cell_frequency_boxplots(
#'   seurat_obj = PAH.obj.Int.Stromal.4.clean,
#'   cell_type_column = "stromal.cell.type.final",
#'   condition_column = "disease.ident",
#'   sample_column = "orig.ident",
#'   output_prefix = "Stromal_Cell_Frequencies_DefaultSort",
#'   output_dir = "./results",
#'   control_level = "Control",
#'   disease_level = "Disease",
#'   paired = FALSE
#' )
#'
#' # Example with reverse alphabetical sorting (Z-A)
#' plot_cell_frequency_boxplots(
#'   seurat_obj = PAH.obj.Int.Stromal.4.clean,
#'   cell_type_column = "stromal.cell.type.final",
#'   condition_column = "disease.ident",
#'   sample_column = "orig.ident",
#'   output_prefix = "Stromal_Cell_Frequencies_ReverseSort",
#'   output_dir = "./results",
#'   control_level = "Control",
#'   disease_level = "Disease",
#'   sort_alphabetically = "INVERSE",
#'   paired = FALSE
#' )
#'
#' # Example with a custom x-axis order
#' custom_order <- c("Adventitial Fibroblast", "Pericyte", "Smooth Muscle Cell") # Example types
#' plot_cell_frequency_boxplots(
#'   seurat_obj = PAH.obj.Int.Stromal.4.clean,
#'   cell_type_column = "stromal.cell.type.final",
#'   condition_column = "disease.ident",
#'   sample_column = "orig.ident",
#'   output_prefix = "Stromal_Cell_Frequencies_CustomOrder",
#'   output_dir = "./results",
#'   control_level = "Control",
#'   disease_level = "Disease",
#'   custom_x_axis_order = custom_order,
#'   paired = FALSE
#' )
#'
#' # Example with a custom y-axis maximum limit of 50%
#' # Note: You might need to adjust star_height if setting y_max_limit.
#' # For y_max_limit = 50, a star_height of e.g., 47 might be appropriate.
#' plot_cell_frequency_boxplots(
#'   seurat_obj = PAH.obj.Int.Stromal.4.clean,
#'   cell_type_column = "stromal.cell.type.final",
#'   condition_column = "disease.ident",
#'   sample_column = "orig.ident",
#'   output_prefix = "Stromal_Cell_Frequencies_CustomYlim",
#'   output_dir = "./results",
#'   control_level = "Control",
#'   disease_level = "Disease",
#'   y_max_limit = 50,
#'   star_height = 47, # Adjusted star_height for the new y_max_limit
#'   paired = FALSE
#' )
#' }
#'
#' @importFrom dplyr group_by summarise mutate select case_when n filter rename bind_rows
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_text scale_fill_manual theme_classic
#'   scale_y_continuous labs theme element_text ggsave position_dodge element_blank
#' @importFrom stats wilcox.test as.formula
#' @importFrom rlang sym
#' @importFrom tidyr pivot_wider
#' @importFrom utils write.table head str
#'
#' @export
plot_cell_frequency_boxplots <- function(seurat_obj,
                                         cell_type_column,
                                         condition_column,
                                         sample_column,
                                         output_prefix,
                                         output_dir = ".",
                                         width = 10,
                                         height = 6,
                                         dpi = 300,
                                         control_level = "CTRL",
                                         disease_level = "PAH",
                                         control_color = "dodgerblue2",
                                         disease_color = "firebrick",
                                         star_height = 95,
                                         y_max_limit = NULL,
                                         paired = FALSE,
                                         x_angle = 45,
                                         plot_title = "Cell Type Frequencies by Condition",
                                         save_results = TRUE,
                                         sort_alphabetically = FALSE, # New parameter
                                         custom_x_axis_order = NULL) { # New parameter

  # --- Input Validation and Setup ---

  # Check required packages
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is needed. Please install it.", call. = FALSE)
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is needed. Please install it.", call. = FALSE)
  }
  if (!requireNamespace("rlang", quietly = TRUE)) {
    stop("Package 'rlang' is needed. Please install it.", call. = FALSE)
  }
  if (paired && !requireNamespace("tidyr", quietly = TRUE)) {
    stop("Package 'tidyr' is needed for paired tests. Please install it.", call. = FALSE)
  }

  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    message("Created output directory: ", output_dir)
  }

  # Verify columns exist in metadata
  meta_data <- seurat_obj@meta.data
  required_cols <- c(cell_type_column, condition_column, sample_column)
  missing_cols <- required_cols[!required_cols %in% colnames(meta_data)]
  if (length(missing_cols) > 0) {
    stop("Column(s) not found in Seurat object metadata: ", paste(missing_cols, collapse = ", "))
  }

  # Ensure condition levels exist
  if (!control_level %in% meta_data[[condition_column]]) {
    stop("Control level '", control_level, "' not found in column '", condition_column, "'.")
  }
  if (!disease_level %in% meta_data[[condition_column]]) {
    stop("Disease level '", disease_level, "' not found in column '", condition_column, "'.")
  }

  # Filter metadata to only include specified conditions
  meta_data <- meta_data[meta_data[[condition_column]] %in% c(control_level, disease_level), ]
  if(nrow(meta_data) == 0) {
    stop("No cells found matching the specified control/disease levels.")
  }

  # --- Frequency Calculation ---
  freq_data <- meta_data %>%
    dplyr::group_by(!!rlang::sym(sample_column),
                    !!rlang::sym(condition_column),
                    !!rlang::sym(cell_type_column)) %>%
    dplyr::summarise(n_cell_type_sample = dplyr::n(), .groups = "drop") %>%
    dplyr::group_by(!!rlang::sym(sample_column), !!rlang::sym(condition_column)) %>%
    dplyr::mutate(n_total_sample = sum(n_cell_type_sample)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(n_total_sample > 0) %>% # Avoid division by zero
    dplyr::mutate(freq = (n_cell_type_sample / n_total_sample) * 100) %>%
    dplyr::select(!!rlang::sym(sample_column),
                  !!rlang::sym(condition_column),
                  !!rlang::sym(cell_type_column),
                  freq)

  # Get unique cell types present in the filtered data for statistics and ordering
  all_cell_types_present <- unique(as.character(freq_data[[cell_type_column]]))
  if (length(all_cell_types_present) == 0) {
    stop("No cell types found after frequency calculation. Check input data and filters.")
  }

  # --- Significance Testing ---
  add_significance_stars <- function(p.value) {
    dplyr::case_when(
      is.na(p.value) ~ "ns",
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE            ~ "ns"
    )
  }

  significance_results <- list()

  for (ct in all_cell_types_present) { # Use all cell types found in data
    ct_data <- freq_data[as.character(freq_data[[cell_type_column]]) == ct, ]
    conditions_present <- unique(ct_data[[condition_column]])
    p_val <- NA

    if (length(conditions_present) == 2) {
      if (paired) {
        paired_data <- ct_data %>%
          tidyr::pivot_wider(id_cols = !!rlang::sym(sample_column),
                             names_from = !!rlang::sym(condition_column),
                             values_from = freq) %>%
          dplyr::filter(!is.na(!!rlang::sym(control_level)) & !is.na(!!rlang::sym(disease_level)))

        if (nrow(paired_data) > 0) { # Rely on wilcox.test for handling small N
          test_result <- tryCatch({
            stats::wilcox.test(paired_data[[control_level]], paired_data[[disease_level]], paired = TRUE)
          }, warning = function(w) {
            message("Warning during paired Wilcoxon test for '", ct, "': ", w$message)
            stats::wilcox.test(paired_data[[control_level]], paired_data[[disease_level]], paired = TRUE) # Re-run to get result despite warning
          }, error = function(e) {
            warning("Error during paired Wilcoxon test for '", ct, "': ", e$message)
            NULL
          })
          if (!is.null(test_result)) p_val <- test_result$p.value
        } else {
          warning("No samples found with data for *both* conditions for cell type '", ct, "' after pivoting. Skipping paired test.")
          p_val <- NA
        }
      } else { # Unpaired
        # Rely on wilcox.test for handling small N and ties
        formula_str <- paste("freq ~", rlang::sym(condition_column))
        test_result <- tryCatch({
          stats::wilcox.test(stats::as.formula(formula_str), data = ct_data)
        }, warning = function(w) {
          message("Warning during unpaired Wilcoxon test for '", ct, "': ", w$message)
          stats::wilcox.test(stats::as.formula(formula_str), data = ct_data) # Re-run
        }, error = function(e) {
          warning("Error during unpaired Wilcoxon test for '", ct, "': ", e$message)
          NULL
        })
        if (!is.null(test_result)) p_val <- test_result$p.value
      }
    } else {
      warning("Cell type '", ct, "' not present in both conditions ('", control_level, "', '", disease_level, "'). Skipping significance test.")
      p_val <- NA
    }

    significance_results[[ct]] <- data.frame(
      cell_type = ct,
      p_value = p_val,
      significance = add_significance_stars(p_val),
      stringsAsFactors = FALSE
    )
  }

  significance_df <- dplyr::bind_rows(significance_results)
  if (nrow(significance_df) > 0) {
    significance_df <- significance_df %>%
      dplyr::rename(!!rlang::sym(cell_type_column) := cell_type)
  } else {
    message("No significance results generated (significance_df is empty).")
    significance_df <- data.frame(
      col1 = character(),
      p_value = numeric(),
      significance = character(),
      stringsAsFactors = FALSE
    )
    names(significance_df)[1] <- cell_type_column
  }

  message("--- Significance Data Frame (Head) ---")
  print(utils::head(significance_df))
  message("--- Structure of Significance Data Frame ---")
  utils::str(significance_df)
  message("--------------------------------------")

  if (save_results && nrow(significance_df) > 0 && any(!is.na(significance_df$p_value))) {
    results_file <- file.path(output_dir, paste0(output_prefix, "_significance_results.txt"))
    tryCatch({
      write.table(significance_df, file = results_file, sep = "\t", row.names = FALSE, quote = FALSE)
      message("Significance results saved to: ", results_file)
    }, error = function(e) {
      warning("Could not save significance results: ", e$message)
    })
  } else if (save_results) {
    message("Significance results file not saved (no valid p-values, results empty, or saving not requested).")
  }

  # --- X-axis Cell Type Ordering ---
  plot_levels <- character(0)

  if (length(all_cell_types_present) > 0) {
    if (!is.null(custom_x_axis_order)) {
      if (!is.character(custom_x_axis_order)) {
        warning("'custom_x_axis_order' is not a character vector. Ignoring and using default A-Z sort.")
        plot_levels <- sort(all_cell_types_present, decreasing = FALSE)
      } else {
        # Cell types from custom order that are present in the data
        ordered_part <- custom_x_axis_order[custom_x_axis_order %in% all_cell_types_present]
        # Cell types in data but not in custom order, sorted A-Z
        remaining_part <- sort(setdiff(all_cell_types_present, custom_x_axis_order), decreasing = FALSE)
        plot_levels <- c(ordered_part, remaining_part)

        # Warn about custom types not found in data
        custom_not_in_data <- setdiff(custom_x_axis_order, all_cell_types_present)
        if (length(custom_not_in_data) > 0) {
          warning("The following cell types in 'custom_x_axis_order' are not present in the data and will be ignored: ",
                  paste(custom_not_in_data, collapse = ", "), ".")
        }
        if (length(plot_levels) == 0 && length(all_cell_types_present) > 0) {
          warning("Custom order resulted in no overlapping cell types for plotting. Defaulting to A-Z sort.")
          plot_levels <- sort(all_cell_types_present, decreasing = FALSE)
        }
      }
    } else { # Use sort_alphabetically
      valid_sort_options <- c(TRUE, FALSE, "INVERSE")
      if (!(sort_alphabetically %in% valid_sort_options ||
            (is.logical(sort_alphabetically) && length(sort_alphabetically) == 1) )) { # check if it's single TRUE/FALSE
        warning(paste0("Invalid value for 'sort_alphabetically': '",
                       as.character(sort_alphabetically),
                       "'. Must be TRUE, FALSE, or 'INVERSE'. Defaulting to FALSE (A-Z sort)."))
        sort_alphabetically_eff <- FALSE
      } else {
        sort_alphabetically_eff <- sort_alphabetically
      }

      if (identical(sort_alphabetically_eff, "INVERSE")) {
        plot_levels <- sort(all_cell_types_present, decreasing = TRUE) # Z-A
      } else { # TRUE or FALSE (default)
        plot_levels <- sort(all_cell_types_present, decreasing = FALSE) # A-Z
      }
    }
  } else {
    # This case should be caught by the earlier stop if all_cell_types_present is empty.
    # If somehow reached, plot_levels remains empty.
    message("No cell types available for ordering (all_cell_types_present is empty).")
  }

  # Apply factor levels if plot_levels were determined
  if (length(plot_levels) > 0) {
    # Filter data to include only cell types that are in plot_levels
    # This ensures that if custom_x_axis_order was very restrictive, we only plot what's intended.
    freq_data <- freq_data[as.character(freq_data[[cell_type_column]]) %in% plot_levels, ]

    if (nrow(freq_data) > 0) {
      freq_data[[cell_type_column]] <- factor(as.character(freq_data[[cell_type_column]]), levels = plot_levels)
    } # If freq_data becomes empty, ggplot will handle it (empty plot)

    if(nrow(significance_df) > 0 && cell_type_column %in% names(significance_df)) {
      significance_df[[cell_type_column]] <- as.character(significance_df[[cell_type_column]]) # Ensure character before filtering
      significance_df <- significance_df[significance_df[[cell_type_column]] %in% plot_levels, ]
      if (nrow(significance_df) > 0) {
        significance_df[[cell_type_column]] <- factor(significance_df[[cell_type_column]], levels = plot_levels)
      }
    }
  } else if (length(all_cell_types_present) > 0) {
    # Fallback if plot_levels is empty but there were cell types. Should ideally not happen with current logic.
    warning("Cell type ordering resulted in no levels, but data was present. Defaulting to A-Z sort of available types for plot.")
    plot_levels <- sort(all_cell_types_present, decreasing = FALSE)
    freq_data[[cell_type_column]] <- factor(as.character(freq_data[[cell_type_column]]), levels = plot_levels)
    if(nrow(significance_df) > 0 && cell_type_column %in% names(significance_df)) {
      significance_df[[cell_type_column]] <- as.character(significance_df[[cell_type_column]])
      significance_df <- significance_df[significance_df[[cell_type_column]] %in% plot_levels, ]
      if(nrow(significance_df) > 0) {
        significance_df[[cell_type_column]] <- factor(significance_df[[cell_type_column]], levels = plot_levels)
      }
    }
  }


  # --- Plotting ---
  condition_colors <- setNames(c(control_color, disease_color),
                               c(control_level, disease_level))

  current_max_freq <- if (nrow(freq_data) > 0 && sum(!is.na(freq_data$freq)) > 0) {
    max(freq_data$freq, na.rm = TRUE)
  } else {
    0
  }

  y_for_stars <- star_height
  plot_y_axis_upper_limit <- NULL

  if (!is.null(y_max_limit) && is.numeric(y_max_limit)) {
    plot_y_axis_upper_limit <- y_max_limit
    if (y_for_stars >= plot_y_axis_upper_limit * 0.98 && y_for_stars > 0) {
      warning(paste0("Parameter 'star_height' (", star_height, ") is close to or above 'y_max_limit' (",
                     y_max_limit, "). Stars might be clipped or overlap. ",
                     "Consider adjusting 'star_height' to be lower, e.g., around ",
                     round(plot_y_axis_upper_limit * 0.95, 1), "."))
    }
  } else {
    plot_y_axis_upper_limit <- max(current_max_freq * 1.1, y_for_stars * 1.05, 10, na.rm = TRUE)
    if (!is.finite(plot_y_axis_upper_limit)) {
      plot_y_axis_upper_limit <- max(100, y_for_stars * 1.05, na.rm = TRUE) # Sensible fallback
      if(!is.finite(plot_y_axis_upper_limit)) plot_y_axis_upper_limit <- 100 # Absolute fallback
    }
  }

  # Check if freq_data has rows and the cell_type_column is not all NA after factoring
  if(nrow(freq_data) == 0 || all(is.na(freq_data[[cell_type_column]]))) {
    warning("No data available for plotting after processing and ordering. An empty plot will be generated.")
    p <- ggplot2::ggplot() +
      ggplot2::labs(title = paste(plot_title, "(No data to display)"), x="Cell Type", y="Frequency (%)") +
      ggplot2::theme_minimal()
  } else {
    p <- ggplot2::ggplot(freq_data,
                         ggplot2::aes(x = !!rlang::sym(cell_type_column),
                                      y = freq,
                                      fill = factor(!!rlang::sym(condition_column), levels = c(control_level, disease_level)))) +
      ggplot2::geom_boxplot(outlier.shape = NA,
                            position = ggplot2::position_dodge(width = 0.8),
                            width = 0.7) +
      ggplot2::scale_fill_manual(values = condition_colors, name = "Condition") +
      ggplot2::scale_y_continuous(limits = c(0, plot_y_axis_upper_limit),
                                  expand = ggplot2::expansion(mult = c(0, 0.05))) + # expand slightly for top
      ggplot2::labs(
        title = plot_title,
        x = "Cell Type",
        y = "Frequency (%)"
      ) +
      ggplot2::theme_classic() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = x_angle, hjust = 1, vjust = 1, size = 10),
        axis.text.y = ggplot2::element_text(size = 10),
        axis.title = ggplot2::element_text(size = 12),
        plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.title = ggplot2::element_text(size = 11),
        legend.text = ggplot2::element_text(size = 10),
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank()
      )

    if (nrow(significance_df) > 0 && any(!is.na(significance_df$p_value))) {
      plot_sig_df <- significance_df[!is.na(significance_df$significance) &
                                       !is.na(significance_df[[cell_type_column]]) & # Ensure cell type is not NA
                                       significance_df$significance != "ns", ] # Only plot actual stars
      if(nrow(plot_sig_df) > 0) {
        p <- p + ggplot2::geom_text(
          data = plot_sig_df,
          ggplot2::aes(
            x = !!rlang::sym(cell_type_column),
            y = y_for_stars,
            label = significance,
            group = NULL,
            fill = NULL
          ),
          inherit.aes = FALSE, # Important to prevent issues with main aes
          size = 5,
          vjust = 0.5
        )
      }
    }
  }

  # --- Save Plots ---
  pdf_file <- file.path(output_dir, paste0(output_prefix, "_boxplot.pdf"))
  png_file <- file.path(output_dir, paste0(output_prefix, "_boxplot.png"))

  tryCatch({
    ggplot2::ggsave(pdf_file, p, width = width, height = height, device = "pdf")
    message("Plot saved to: ", pdf_file)
  }, error = function(e) {
    warning("Could not save PDF plot: ", e$message)
  })

  tryCatch({
    ggplot2::ggsave(png_file, p, width = width, height = height, dpi = dpi, device = "png", bg = "white")
    message("Plot saved to: ", png_file)
  }, error = function(e) {
    warning("Could not save PNG plot: ", e$message)
  })

  # --- Return Plot Object ---
  return(p)
}

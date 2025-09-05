#' Perform and Plot Cell Cycle Analysis Results
#'
#' @description
#' Performs cell cycle scoring on a Seurat object (if not already done), conducts
#' various statistical tests comparing cell cycle phases across cell types and
#' disease conditions, generates a series of plots visualizing the results, and
#' saves both the plots and statistical outputs.
#'
#' @param seurat_obj A Seurat object. Should contain an assay with expression data
#'   (e.g., "RNA"). Metadata should include columns specified by `cell_type_column`
#'   and `disease_column`.
#' @param cell_type_column Character string. The name of the column in
#'   `seurat_obj@meta.data` containing cell type annotations.
#' @param disease_column Character string. The name of the column in
#'   `seurat_obj@meta.data` containing disease status annotations (must have exactly two levels, e.g., "PAH", "CTRL").
#' @param save_path Character string. Path to the directory where plots and
#'   statistical results will be saved. The directory will be created if it
#'   doesn't exist.
#' @param file_prefix Character string. A prefix to be added to all saved output files.
#' @param plot_width Numeric. Width of the saved plots in inches.
#' @param plot_height Numeric. Height of the saved plots in inches.
#' @param reduction_type Character string. Dimensionality reduction to use for plotting
#'   (e.g., "umap", "tsne"). Defaults to "umap".
#' @param s_genes Character vector. Genes defining the S phase. Defaults to
#'   `Seurat::cc.genes.updated.2019$s.genes`.
#' @param g2m_genes Character vector. Genes defining the G2/M phase. Defaults to
#'   `Seurat::cc.genes.updated.2019$g2m.genes`.
#' @param dot_plot_genes Character vector. Genes to include in the cell cycle dot plot
#'   (Plot 10). Defaults to `c("PCNA", "MKI67", "TOP2A", "CCNB1", "CDK1", "MCM5")`.
#' @param colors_use Optional. A character vector of colors to use for discrete plots
#'   (e.g., cell cycle phases, conditions). If `NULL` (default), attempts to use
#'   `scCustomize::scCustomize_Palette`. If `scCustomize` is not available,
#'   uses default `ggplot2` colors (`scales::hue_pal`).
#'
#' @return Returns the input `seurat_obj`, potentially updated with cell cycle scores
#'   ("S.Score", "G2M.Score", "Phase" columns in metadata).
#'
#' @details
#' The function performs the following steps:
#' 1.  **Cell Cycle Scoring:** Scores cells for S and G2/M phases using `Seurat::CellCycleScoring` if "Phase" column is not already present in metadata.
#' 2.  **Metadata Validation:** Checks for required columns and ensures the disease column has exactly two levels. Converts relevant columns to factors.
#' 3.  **Statistical Analysis:**
#'     * Overall Chi-squared tests (Cell Type vs Phase, Disease vs Phase).
#'     * Pairwise Fisher's exact tests (Cell Type vs Phase).
#'     * Proportion test (or Fisher's if appropriate) for Disease vs Cycling status (G1 vs S/G2M).
#'     * Odds ratio calculation for Disease vs Cycling status.
#'     * Logistic regression modeling Cycling status ~ Disease * Cell Type (optional, may fail with sparse data).
#'     * Pairwise Fisher's exact tests for Disease vs Cycling status within each Cell Type.
#' 4.  **Save Statistics:** Saves the results of statistical tests to text and CSV files.
#' 5.  **Generate Plots:** Creates various plots including UMAPs, violin plots, ridge plots, bar plots, feature plots, a dot plot, density plots, and a logistic regression visualization plot.
#' 6.  **Save Plots:** Saves generated plots as both PDF and PNG files.
#'
#' It relies on several packages listed in the `@importFrom` tags. Ensure these are installed.
#' The `scCustomize` package is suggested for enhanced color palettes but not required.
#'
#' @importFrom Seurat CellCycleScoring DimPlot VlnPlot RidgePlot FeatureScatter DotPlot FeaturePlot GetAssayData Assays DefaultAssay
#' @importFrom ggplot2 ggplot aes geom_bar geom_density geom_bar geom_point theme element_text ylab xlab ggtitle scale_fill_manual scale_color_manual theme_minimal ggsave position_dodge position_jitterdodge unit
#' @importFrom dplyr %>% mutate select filter group_by summarise arrange bind_rows rename all_of pull
#' @importFrom tidyr pivot_longer
#' @importFrom rstatix pairwise_fisher_test
#' @importFrom epitools oddsratio
#' @importFrom broom tidy
#' @importFrom modelr add_predictions
#' @importFrom stats chisq.test fisher.test prop.test glm p.adjust as.formula sd na.omit
#' @importFrom utils capture.output write.csv head tail
#' @importFrom scales hue_pal
#' @importFrom grDevices pdf png dev.off dev.cur
#' @importFrom methods is
#'
#' @export
plot_cell_cycle_analysis <- function(seurat_obj,
                                     cell_type_column, # New parameter for cell type column name
                                     disease_column,   # New parameter for disease column name
                                     save_path,
                                     file_prefix,
                                     plot_width,
                                     plot_height,
                                     reduction_type = "umap",
                                     s_genes = Seurat::cc.genes.updated.2019$s.genes,
                                     g2m_genes = Seurat::cc.genes.updated.2019$g2m.genes,
                                     dot_plot_genes = c("PCNA", "MKI67", "TOP2A", "CCNB1", "CDK1", "MCM5"),
                                     colors_use = NULL) { # Default colors_use to NULL

  # No library() calls inside the function; dependencies handled by @importFrom

  message("--- Starting Cell Cycle Analysis Function ---")

  # --- Dependency Check & Color Setup ---
  use_scCustomize_colors <- FALSE
  if (is.null(colors_use)) {
    message("Checking for optional package 'scCustomize' for color palettes...")
    if (requireNamespace("scCustomize", quietly = TRUE)) {
      message("scCustomize found. Will attempt to use its palettes.")
      use_scCustomize_colors <- TRUE
    } else {
      message("scCustomize not found. Using default ggplot2 colors.")
    }
  } else {
    message("Using user-provided 'colors_use'.")
  }

  # Define default phase colors (will be overridden if scCustomize used later)
  default_phase_colors <- scales::hue_pal()(3) # G1, S, G2M
  names(default_phase_colors) <- c("G1", "S", "G2M")


  # --- Perform Cell Cycle Scoring (if not already done) ---
  message("Checking for existing cell cycle scoring...")
  if (!"Phase" %in% colnames(seurat_obj@meta.data)) {
    message("Performing cell cycle scoring...")
    # Check if genes are present in the object
    genes_present <- rownames(GetAssayData(seurat_obj, assay = DefaultAssay(seurat_obj)))
    s_genes_present <- intersect(s_genes, genes_present)
    g2m_genes_present <- intersect(g2m_genes, genes_present)

    if(length(s_genes_present) < 5 || length(g2m_genes_present) < 5) {
      warning("Too few S or G2M genes found in the object. Cell cycle scoring might be unreliable or fail.", call. = FALSE)
    }

    seurat_obj <- tryCatch({
      Seurat::CellCycleScoring(
        object = seurat_obj,
        s.features = s_genes_present, # Use only genes present
        g2m.features = g2m_genes_present, # Use only genes present
        set.ident = FALSE # Do not change the primary identity
      )
    }, error = function(e) {
      stop("Cell cycle scoring failed. Error: ", e$message, call. = FALSE)
    })
    message("Cell cycle scoring complete.")
  } else {
    message("Cell cycle scoring already found in metadata. Skipping scoring.")
  }

  # --- Metadata Validation and Formatting ---
  message("Checking and formatting metadata columns...")
  required_cols <- c(cell_type_column, disease_column, "Phase", "S.Score", "G2M.Score") # Added score columns check
  if (!all(required_cols %in% colnames(seurat_obj@meta.data))) {
    missing_req <- setdiff(required_cols, colnames(seurat_obj@meta.data))
    stop(paste("Seurat object metadata must contain the following columns:",
               paste(required_cols, collapse = ", "),
               "\nMissing columns:", paste(missing_req, collapse=", ")))
  }

  # Check that the disease column has exactly two levels
  disease_levels <- unique(na.omit(seurat_obj@meta.data[[disease_column]]))
  if (length(disease_levels) != 2) {
    stop(paste("The disease column '", disease_column, "' must have exactly two unique non-NA values (e.g., PAH and CTRL). Found:",
               paste(disease_levels, collapse=", ")))
  }

  # Ensure cell type and disease columns are factors
  seurat_obj@meta.data[[cell_type_column]] <- as.factor(seurat_obj@meta.data[[cell_type_column]])
  seurat_obj@meta.data[[disease_column]] <- as.factor(seurat_obj@meta.data[[disease_column]])

  # Ensure Phase is a factor with consistent ordering
  seurat_obj@meta.data <- seurat_obj@meta.data %>%
    dplyr::mutate(Phase = factor(Phase, levels = c("G1", "S", "G2M")))

  message("Metadata formatting complete.")

  # --- Setup Save Directory ---
  message("Checking/creating save directory: ", save_path)
  if (!dir.exists(save_path)) {
    dir.create(save_path, recursive = TRUE)
    message("Created directory: ", save_path)
  } else {
    message("Save directory already exists.")
  }


  # --- Perform Statistical Analysis ---
  message("Performing statistical analysis...")

  # Convert metadata to a data frame for easier manipulation
  metadata_df <- as.data.frame(seurat_obj@meta.data)

  # 1. Overall Chi-squared Test (Cell Type vs Phase)
  message("  Performing overall Chi-squared test (Cell Type vs Phase)...")
  cell_type_phase_counts <- table(metadata_df[[cell_type_column]], metadata_df$Phase)
  overall_celltype_phase_test <- tryCatch({
    stats::chisq.test(cell_type_phase_counts, simulate.p.value = TRUE)
  }, error = function(e) {
    message("    Warning: Could not perform overall Chi-squared test for Cell Type vs Phase. Error: ", e$message)
    return(NULL)
  })
  message("  Overall Chi-squared test (Cell Type vs Phase) complete.")


  # 2. Pairwise Fisher's Exact Tests (Cell Type vs Phase)
  message("  Performing pairwise Fisher's exact tests (Cell Type vs Phase)...")
  pairwise_celltype_phase_results <- tryCatch({
    # Use .data pronoun for robustness
    rstatix::pairwise_fisher_test(metadata_df, formula = Phase ~ .data[[cell_type_column]], p.adjust.method = "BH")
  }, error = function(e) {
    message("    Warning: Could not perform pairwise Fisher's exact tests for Cell Type vs Phase. Error: ", e$message)
    return(NULL)
  })
  message("  Pairwise Fisher's exact tests (Cell Type vs Phase) complete.")


  # 3. Overall Chi-squared Test (Disease vs Phase)
  message("  Performing overall Chi-squared test (Disease vs Phase)...")
  disease_phase_counts <- table(metadata_df[[disease_column]], metadata_df$Phase)
  overall_disease_phase_test <- tryCatch({
    stats::chisq.test(disease_phase_counts, simulate.p.value = TRUE)
  }, error = function(e) {
    message("    Warning: Could not perform overall Chi-squared test for Disease vs Phase. Error: ", e$message)
    return(NULL)
  })
  message("  Overall Chi-squared test (Disease vs Phase) complete.")


  # 4. Proportion Test / Fisher's Exact Test (Disease vs Cycling Status)
  message("  Performing proportion/Fisher test (Disease vs Cycling Status)...")
  metadata_df <- metadata_df %>%
    dplyr::mutate(Is_Cycling = ifelse(Phase %in% c("S", "G2M"), "Cycling", "G1")) %>%
    dplyr::mutate(Is_Cycling = as.factor(Is_Cycling))

  disease_cycling_counts <- table(metadata_df[[disease_column]], metadata_df$Is_Cycling)
  disease_cycling_test <- tryCatch({
    if (nrow(disease_cycling_counts) == 2 && ncol(disease_cycling_counts) == 2 && all(disease_cycling_counts > 0)) {
      stats::prop.test(disease_cycling_counts) # Use prop.test for 2x2 with counts > 0
    } else if (nrow(disease_cycling_counts) == 2 && ncol(disease_cycling_counts) == 2) {
      message("    Note: Using Fisher's exact test for Disease vs Cycling due to low counts or zero cells in a category.")
      stats::fisher.test(disease_cycling_counts) # Use Fisher if counts are low/zero
    } else {
      message("    Warning: Contingency table for Disease vs Cycling is not 2x2. Skipping proportion/Fisher test.")
      return(NULL)
    }
  }, error = function(e) {
    message("    Warning: Could not perform proportion/Fisher test for Disease vs Cycling. Error: ", e$message)
    return(NULL)
  })
  message("  Proportion/Fisher test (Disease vs Cycling Status) complete.")


  # 5. Odds Ratio (Disease vs Cycling Status)
  message("  Calculating Odds Ratio (Disease vs Cycling Status)...")
  odds_ratio_cycling <- tryCatch({
    if (nrow(disease_cycling_counts) == 2 && ncol(disease_cycling_counts) == 2) {
      # Ensure table has correct dimensions before passing to oddsratio
      epitools::oddsratio(disease_cycling_counts, method = "wald")
    } else {
      message("    Warning: Contingency table for Disease vs Cycling is not 2x2. Skipping odds ratio calculation.")
      return(NULL)
    }
  }, error = function(e) {
    message("    Warning: Could not calculate odds ratio for Disease vs Cycling. Error: ", e$message)
    return(NULL)
  })
  message("  Odds Ratio calculation complete.")


  # 6. Logistic Regression (Cycling Status ~ Disease * Cell Type) - Might fail with low counts
  message("  Attempting to fit Logistic Regression model (Cycling Status ~ Disease * Cell Type)...")
  logistic_model_interaction <- NULL # Initialize
  if (length(unique(metadata_df$Is_Cycling)) > 1) {
    logistic_model_interaction <- tryCatch({
      formula_str <- paste("Is_Cycling ~", disease_column, "*", cell_type_column)
      stats::glm(stats::as.formula(formula_str), data = metadata_df, family = "binomial")
    }, error = function(e) {
      message("    Warning: Could not fit logistic regression model with interaction. Error: ", e$message)
      return(NULL)
    })
  } else {
    message("    Warning: No variation in cycling status. Skipping logistic regression.")
  }
  message("  Logistic Regression model fitting complete (check warnings/errors).")


  # 7. Pairwise Fisher's Exact Tests (Disease vs Cycling) within each Cell Type
  message("  Performing pairwise Fisher's exact tests for Disease vs Cycling within each Cell Type...")

  cell_types <- levels(metadata_df[[cell_type_column]])
  pairwise_disease_cycling_by_celltype_results_list <- list()

  for (cell_type in cell_types) {
    message("    Processing cell type: ", cell_type)
    subset_df <- metadata_df %>%
      dplyr::filter(.data[[cell_type_column]] == cell_type)

    # Check for sufficient variation within the subset
    if (length(unique(subset_df[[disease_column]])) == 2 && length(unique(subset_df$Is_Cycling)) == 2) {
      cell_type_disease_cycling_counts <- table(subset_df[[disease_column]], subset_df$Is_Cycling)

      # Ensure the table is actually 2x2 before testing
      if(nrow(cell_type_disease_cycling_counts) == 2 && ncol(cell_type_disease_cycling_counts) == 2) {
        fisher_test_result <- tryCatch({
          stats::fisher.test(cell_type_disease_cycling_counts)
        }, error = function(e) {
          message("      Warning: Fisher's test failed for ", cell_type, ". Error: ", e$message)
          return(NULL)
        })

        if (!is.null(fisher_test_result)) {
          pairwise_disease_cycling_by_celltype_results_list[[cell_type]] <- data.frame(
            Cell_Type = cell_type,
            p_value = fisher_test_result$p.value,
            odds_ratio_est = if("estimate" %in% names(fisher_test_result)) fisher_test_result$estimate else NA_real_, # Fisher OR estimate might not always be present
            conf_int_low = if("conf.int" %in% names(fisher_test_result)) fisher_test_result$conf.int[1] else NA_real_,
            conf_int_high = if("conf.int" %in% names(fisher_test_result)) fisher_test_result$conf.int[2] else NA_real_,
            stringsAsFactors = FALSE
          )
          message("      Fisher's test successful for ", cell_type)
        }
      } else {
        message("    Skipping Fisher's test for ", cell_type, ": Table is not 2x2 after subsetting.")
      }
    } else {
      message("    Skipping Fisher's test for ", cell_type, ": Not enough variation in disease or cycling status.")
    }
  }

  # Combine results and apply multiple testing correction
  pairwise_disease_cycling_by_celltype_results <- dplyr::bind_rows(pairwise_disease_cycling_by_celltype_results_list)

  if (nrow(pairwise_disease_cycling_by_celltype_results) > 0) {
    pairwise_disease_cycling_by_celltype_results <- pairwise_disease_cycling_by_celltype_results %>%
      dplyr::mutate(p_adj = stats::p.adjust(p_value, method = "BH")) %>%
      dplyr::arrange(p_adj) # Arrange by adjusted p-value
    message("Pairwise disease vs cycling tests within cell types complete. Results combined and adjusted.")
  } else {
    message("No pairwise disease vs cycling tests could be performed within cell types.")
    pairwise_disease_cycling_by_celltype_results <- NULL # Ensure it's NULL if empty
  }


  # --- Save Statistical Results to Files ---
  message("Saving statistical results...")

  # Helper to save statistical output (captures print output)
  save_stat_output <- function(test_result, file_suffix, save_path, file_prefix) {
    if (!is.null(test_result)) {
      # Use tryCatch for saving robustness
      tryCatch({
        # Use base::print to ensure S3 method dispatch works correctly if needed
        output_text <- utils::capture.output(base::print(test_result))
        file_name <- file.path(save_path, paste0(file_prefix, file_suffix))
        # Use base::writeLines explicitly, although not strictly necessary
        base::writeLines(output_text, con = file_name)
        message("  Saved ", file_name)
      }, error = function(e){
        message("    Warning: Failed to save ", file_suffix, ". Error: ", e$message)
      })
    } else {
      message("  Skipping saving ", file_suffix, ": Test result is NULL.")
    }
  }

  # Helper to save data frame results
  save_df_output <- function(df_result, file_suffix, save_path, file_prefix) {
    if (!is.null(df_result) && nrow(df_result) > 0) {
      tryCatch({
        file_name <- file.path(save_path, paste0(file_prefix, file_suffix))
        utils::write.csv(df_result, file = file_name, row.names = FALSE)
        message("  Saved ", file_name)
      }, error = function(e){
        message("    Warning: Failed to save ", file_suffix, ". Error: ", e$message)
      })
    } else {
      message("  Skipping saving ", file_suffix, ": Results are NULL or empty.")
    }
  }


  # Save overall Chi-squared tests
  save_stat_output(overall_celltype_phase_test, "_overall_celltype_phase_chisq_test.txt", save_path, file_prefix)
  save_stat_output(overall_disease_phase_test, "_overall_disease_phase_chisq_test.txt", save_path, file_prefix)

  # Save pairwise Fisher's results (Cell Type vs Phase)
  save_df_output(pairwise_celltype_phase_results, "_pairwise_celltype_phase_fisher_tests.csv", save_path, file_prefix)

  # Save Disease vs Cycling test
  save_stat_output(disease_cycling_test, "_disease_cycling_proportion_test.txt", save_path, file_prefix)

  # Save Odds Ratio
  save_stat_output(odds_ratio_cycling, "_disease_cycling_odds_ratio.txt", save_path, file_prefix)

  # Save logistic regression summary and tidy results
  if (!is.null(logistic_model_interaction)) {
    save_stat_output(summary(logistic_model_interaction), "_logistic_regression_summary.txt", save_path, file_prefix)
    tidy_logistic_results <- broom::tidy(logistic_model_interaction)
    save_df_output(tidy_logistic_results, "_logistic_regression_tidy.csv", save_path, file_prefix)
  } else {
    message("  Skipping saving logistic regression results: Model is NULL.")
  }

  # Save pairwise Disease vs Cycling results by Cell Type
  save_df_output(pairwise_disease_cycling_by_celltype_results, "_pairwise_disease_cycling_by_celltype.csv", save_path, file_prefix)

  message("Statistical results saving complete.")


  # --- Generate and Save Plots ---
  message("Generating and saving plots...")

  # Determine number of phases and conditions for palettes
  num_phases <- length(levels(seurat_obj@meta.data$Phase))
  num_conditions <- length(levels(seurat_obj@meta.data[[disease_column]]))

  # Get colors for discrete plots
  message("  Determining colors for discrete plots...")
  if (!is.null(colors_use)) {
    # Use provided colors if available
    phase_colors <- colors_use[1:num_phases]
    condition_colors <- colors_use[1:num_conditions] # Assuming enough colors provided
    names(phase_colors) <- levels(seurat_obj@meta.data$Phase)
    names(condition_colors) <- levels(seurat_obj@meta.data[[disease_column]])
    message("    Using provided 'colors_use'.")
  } else if (use_scCustomize_colors) {
    # Try using scCustomize
    phase_colors <- scCustomize::scCustomize_Palette(num_groups = num_phases)
    condition_colors <- scCustomize::scCustomize_Palette(num_groups = num_conditions)
    names(phase_colors) <- levels(seurat_obj@meta.data$Phase)
    names(condition_colors) <- levels(seurat_obj@meta.data[[disease_column]])
    message("    Using scCustomize palettes.")
  } else {
    # Fallback to ggplot defaults
    phase_colors <- default_phase_colors # Use pre-defined default
    condition_colors <- scales::hue_pal()(num_conditions)
    names(condition_colors) <- levels(seurat_obj@meta.data[[disease_column]])
    message("    Using default ggplot2 palettes.")
  }
  message("  Color palettes determined.")


  # Plot 1: UMAP/t-SNE colored by Cell Cycle Phase
  message("  Generating Plot 1: UMAP colored by Phase...")
  p1_umap_phase <- Seurat::DimPlot(seurat_obj, reduction = reduction_type, group.by = "Phase", label = TRUE, cols = phase_colors) +
    ggplot2::ggtitle("UMAP colored by Cell Cycle Phase")
  message("  Plot 1 generated.")

  # Plot 2: Distribution of S.Score by Cell Cycle Phase (Violin Plot)
  message("  Generating Plot 2: S.Score Violin Plot...")
  p2_s_score_vln <- Seurat::VlnPlot(seurat_obj, features = "S.Score", group.by = "Phase", cols = phase_colors, pt.size = 0) +
    ggplot2::ggtitle("S Phase Score Distribution by Phase")
  message("  Plot 2 generated.")

  # Plot 3: Distribution of G2M.Score by Cell Cycle Phase (Violin Plot)
  message("  Generating Plot 3: G2M.Score Violin Plot...")
  p3_g2m_score_vln <- Seurat::VlnPlot(seurat_obj, features = "G2M.Score", group.by = "Phase", cols = phase_colors, pt.size = 0) +
    ggplot2::ggtitle("G2M Phase Score Distribution by Phase")
  message("  Plot 3 generated.")

  # Plot 4: Ridge plots for S.Score and G2M.Score by Phase
  message("  Generating Plot 4: Ridge Plot...")
  p4_scores_ridge <- Seurat::RidgePlot(seurat_obj, features = c("S.Score", "G2M.Score"), group.by = "Phase", cols = phase_colors) +
    ggplot2::ggtitle("Cell Cycle Score Distribution by Phase (Ridge Plot)")
  message("  Plot 4 generated.")

  # Plot 5: FeatureScatter of S.Score vs G2M.Score
  message("  Generating Plot 5: FeatureScatter (S vs G2M Score)...")
  p5_feature_scatter <- Seurat::FeatureScatter(seurat_obj, feature1 = "S.Score", feature2 = "G2M.Score", group.by = "Phase", cols = phase_colors) +
    ggplot2::ggtitle("S Score vs G2M Score colored by Phase")
  message("  Plot 5 generated.")

  # Plot 6: Cell Cycle Phase Distribution within Cell Types (Bar Plot)
  message("  Generating Plot 6: Cell Cycle Phase Distribution by Cell Type (Bar Plot)...")
  p6_celltype_phase_barplot <- ggplot2::ggplot(seurat_obj@meta.data, ggplot2::aes(x = .data[[cell_type_column]], fill = Phase)) +
    ggplot2::geom_bar(position = "fill") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::ylab("Proportion of Cells") +
    ggplot2::xlab(cell_type_column) + # Set x-axis label dynamically
    ggplot2::ggtitle(paste0("Cell Cycle Phase Distribution within ", cell_type_column, "\n(Overall Chi-sq p = ",
                            format_p_value_helper(overall_celltype_phase_test$p.value), ")")) + # Use helper
    ggplot2::scale_fill_manual(values = phase_colors) # Fill by Phase, use phase palette
  message("  Plot 6 generated.")


  # Plot 7: Cell Cycle Phase Distribution by Condition (Bar Plot)
  message("  Generating Plot 7: Cell Cycle Phase Distribution by Condition (Bar Plot)...")
  p7_condition_phase_barplot <- ggplot2::ggplot(seurat_obj@meta.data, ggplot2::aes(x = .data[[disease_column]], fill = Phase)) +
    ggplot2::geom_bar(position = "fill") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::ylab("Proportion of Cells") +
    ggplot2::xlab(disease_column) + # Set x-axis label dynamically
    ggplot2::ggtitle(paste0("Cell Cycle Phase Distribution by ", disease_column, "\n(Overall Chi-sq p = ",
                            format_p_value_helper(overall_disease_phase_test$p.value), ")")) + # Use helper
    ggplot2::scale_fill_manual(values = phase_colors) # Fill by Phase, use phase palette
  message("  Plot 7 generated.")


  # Plot 8: FeaturePlot of S.Score on UMAP/t-SNE
  message("  Generating Plot 8: S.Score FeaturePlot...")
  p8_umap_s_score <- Seurat::FeaturePlot(seurat_obj, features = "S.Score", reduction = reduction_type) +
    ggplot2::ggtitle("UMAP colored by S Phase Score")
  message("  Plot 8 generated.")

  # Plot 9: FeaturePlot of G2M.Score on UMAP/t-SNE
  message("  Generating Plot 9: G2M.Score FeaturePlot...")
  p9_umap_g2m_score <- Seurat::FeaturePlot(seurat_obj, features = "G2M.Score", reduction = reduction_type) +
    ggplot2::ggtitle("UMAP colored by G2M Phase Score")
  message("  Plot 9 generated.")

  # Plot 10: Dot Plot of key cell cycle genes
  message("  Generating Plot 10: Dot Plot of Key Cell Cycle Genes...")
  genes_for_dot_plot_present <- intersect(dot_plot_genes, rownames(seurat_obj))

  if (length(genes_for_dot_plot_present) > 0) {
    p10_dotplot_genes <- Seurat::DotPlot(seurat_obj, features = genes_for_dot_plot_present, group.by = "Phase", cols = "RdYlBu") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
      ggplot2::ggtitle("Expression of Key Cell Cycle Genes by Phase")
    message("  Plot 10 generated.")
  } else {
    p10_dotplot_genes <- NULL
    message("  Skipping Plot 10: None of the specified cell cycle genes for the dot plot were found.")
  }


  # Plot 11: Stacked Density Plot of S and G2M Scores by Cell Type
  message("  Generating Plot 11: Stacked Density Plot of S and G2M Scores by Cell Type...")
  metadata_df_long <- seurat_obj@meta.data %>%
    dplyr::select(dplyr::all_of(c(cell_type_column, "S.Score", "G2M.Score", "Phase"))) %>%
    tidyr::pivot_longer(cols = c(S.Score, G2M.Score), names_to = "Score_Type", values_to = "Score")

  # Use as.formula for facet_wrap
  p11_density_scores_by_celltype <- ggplot2::ggplot(metadata_df_long, ggplot2::aes(x = Score, fill = Score_Type)) +
    ggplot2::geom_density(alpha = 0.7) +
    ggplot2::facet_wrap(stats::as.formula(paste("~", cell_type_column)), scales = "free_y") +
    ggplot2::ggtitle(paste0("Density Plot of S and G2M Scores by ", cell_type_column)) +
    ggplot2::theme_minimal() +
    ggplot2::scale_fill_manual(values = if(use_scCustomize_colors) tryCatch(scCustomize::scCustomize_Palette(num_groups = 2), error=function(e) scales::hue_pal()(2)) else scales::hue_pal()(2)) # Added tryCatch for safety
  message("  Plot 11 generated.")


  # Plot 12: Visualize Logistic Regression Interaction Effect
  message("  Generating Plot 12: Logistic Regression Interaction Plot...")
  p12_logistic_interaction <- NULL # Initialize
  if (!is.null(logistic_model_interaction)) {
    predicted_data <- tryCatch({
      metadata_df %>%
        modelr::add_predictions(logistic_model_interaction, type = "response") %>%
        dplyr::rename(Predicted_Cycling_Prob = pred)
    }, error = function(e){
      message("    Warning: Failed to add predictions from logistic model. Skipping plot 12. Error: ", e$message)
      return(NULL)
    })

    if (!is.null(predicted_data)) {
      summary_predicted_data <- predicted_data %>%
        dplyr::group_by(.data[[cell_type_column]], .data[[disease_column]]) %>%
        dplyr::summarise(Mean_Predicted_Cycling_Prob = mean(Predicted_Cycling_Prob, na.rm = TRUE), .groups = 'drop')

      p12_logistic_interaction <- ggplot2::ggplot(summary_predicted_data,
                                                  ggplot2::aes(x = .data[[cell_type_column]],
                                                               y = Mean_Predicted_Cycling_Prob,
                                                               fill = .data[[disease_column]])) +
        ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge(width = 0.8), width = 0.7) +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
        ggplot2::ylab("Mean Predicted Probability of Being Cycling (S or G2M)") +
        ggplot2::xlab(cell_type_column) +
        ggplot2::ggtitle(paste0("Predicted Cycling Probability by ", cell_type_column, " and ", disease_column)) +
        ggplot2::scale_fill_manual(values = condition_colors)
      message("  Plot 12 generated.")
    }
  } else {
    message("  Skipping Plot 12: Logistic regression model was not successfully fitted.")
  }


  # --- Save plots ---
  message("Saving generated plots...")
  plot_list <- list(
    p1_umap_phase = p1_umap_phase,
    p2_s_score_vln = p2_s_score_vln,
    p3_g2m_score_vln = p3_g2m_score_vln,
    p4_scores_ridge = p4_scores_ridge,
    p5_feature_scatter = p5_feature_scatter,
    p6_celltype_phase_barplot = p6_celltype_phase_barplot,
    p7_condition_phase_barplot = p7_condition_phase_barplot,
    p8_umap_s_score = p8_umap_s_score,
    p9_umap_g2m_score = p9_umap_g2m_score,
    p10_dotplot_genes = p10_dotplot_genes, # Will be NULL if skipped
    p11_density_scores_by_celltype = p11_density_scores_by_celltype,
    p12_logistic_interaction = p12_logistic_interaction # Will be NULL if skipped
  )

  # Remove any NULL plots from the list before saving
  plot_list <- Filter(Negate(is.null), plot_list)


  for (plot_name in names(plot_list)) {
    current_plot <- plot_list[[plot_name]]
    tryCatch({
      # Save as PDF
      pdf_filename_full <- file.path(save_path, paste0(file_prefix, "_", plot_name, ".pdf"))
      ggplot2::ggsave(filename = pdf_filename_full, plot = current_plot, width = plot_width, height = plot_height)

      # Save as PNG
      png_filename_full <- file.path(save_path, paste0(file_prefix, "_", plot_name, ".png"))
      ggplot2::ggsave(filename = png_filename_full, plot = current_plot, width = plot_width, height = plot_height, dpi = 300)
      message("  Saved ", plot_name)
    }, error = function(e) {
      message("    Warning: Failed to save plot '", plot_name, "'. Error: ", e$message)
    })
  }

  message("Plots and statistical results saved successfully to ", save_path)

  # Return the updated Seurat object
  message("--- Cell Cycle Analysis Function Complete ---")
  return(seurat_obj)
}


#' Format P-value Helper
#'
#' @description Internal helper function to format p-values for display. Not exported.
#'
#' @param p Numeric p-value.
#'
#' @return Character string representation of the formatted p-value.
#'
#' @keywords internal
#' @noRd
format_p_value_helper <- function(p) {
  # Check if p is NULL or NA first
  if (is.null(p) || is.na(p)) return("N/A")
  # Proceed with formatting if p is a valid number
  if (p < 0.001) return("< 0.001")
  # Use sprintf for controlled formatting
  if (p < 0.01) return(sprintf("%.3f", p)) # 3 decimal places
  if (p < 0.05) return(sprintf("%.2f", p)) # 2 decimal places
  return(sprintf("%.2f", p)) # Default to 2 decimal places
}


# Example Usage (remove or comment out in final package script):
# Assuming 'your_seurat_object' is your Seurat object
#
# updated_seurat_object <- plot_cell_cycle_analysis(
#   seurat_obj = your_seurat_object,
#   cell_type_column = "cell_type_annotation", # Specify your cell type column name
#   disease_column = "disease_status",      # Specify your disease column name (e.g., "CTRL", "PAH")
#   save_path = "results/cell_cycle",       # Specify save directory
#   file_prefix = "experiment1",            # Specify file prefix
#   plot_width = 10,                        # Specify plot width
#   plot_height = 8                         # Specify plot height
# )

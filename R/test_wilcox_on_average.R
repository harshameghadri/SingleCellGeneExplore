#' Perform Wilcoxon Test and Calculate Fold Change on Averaged Expression Data
#'
#' @description
#' This function takes a data frame or tibble containing expression data, groups it by gene,
#' and then calculates the average expression for two specified conditions ("PAH" and "CTRL").
#' It computes the log2 fold change between these conditions, performs a Wilcoxon rank-sum test
#' (Mann-Whitney U test) to compare the expression distributions between the conditions for each gene,
#' and calculates adjusted p-values using various methods (BH, FDR, Bonferroni).
#'
#' @param data A data frame or tibble. It must contain at least the following columns:
#'   \itemize{
#'     \item `gene`: Identifier for the gene (or feature).
#'     \item `average_expression`: Numeric values representing expression levels (e.g., averaged per sample or cell).
#'     \item `disease`: A factor or character vector indicating the condition for each expression value,
#'            expected to contain at least "PAH" and "CTRL".
#'   }
#'
#' @return A tibble summarizing the results for each gene. The output tibble includes:
#'   \itemize{
#'     \item `gene`: The gene identifier.
#'     \item `Average_PAH`: Mean of `average_expression` for the "PAH" group.
#'     \item `Average_CTRL`: Mean of `average_expression` for the "CTRL" group.
#'     \item `log2FoldChange`: Log2 ratio of `Average_PAH` / `Average_CTRL`.
#'     \item `p_value`: Raw p-value from the Wilcoxon rank-sum test comparing "PAH" vs "CTRL".
#'            `NA` if the test could not be performed (e.g., insufficient data in one group) or resulted in an error.
#'     \item `adjusted_p_value.BH`: Benjamini-Hochberg adjusted p-value.
#'     \item `adjusted_p_value.FDR`: Alias for Benjamini-Hochberg adjusted p-value (method = "fdr" is the same as "BH").
#'     \item `adjusted_p_value.Bonferroni`: Bonferroni adjusted p-value.
#'   }
#'   The results are filtered to remove rows where `log2FoldChange` or `p_value` are `NA`, `NaN`, or infinite.
#'
#' @details
#' The function uses `dplyr` for data manipulation and `stats::wilcox.test` for the statistical test.
#' It assumes the `average_expression` column represents values suitable for direct comparison between
#' the "PAH" and "CTRL" groups within each `gene`.
#' The Wilcoxon test is performed using `exact = FALSE` for potentially faster computation with larger datasets,
#' enabling the use of normal approximation with continuity correction.
#' Error handling is included for the Wilcoxon test; if an error occurs (e.g., due to insufficient data points
#' or identical values in both groups), the `p_value` for that gene will be `NA`.
#' Note: The function name `Test.wilcox.onAverage` uses dots. Consider renaming using underscores
#' (`test_wilcox_on_average`) or camelCase (`testWilcoxOnAverage`) for better style consistency in packages.
#'
#' @importFrom dplyr group_by summarise mutate filter %>%
#' @importFrom stats wilcox.test p.adjust
#' @importFrom rlang .data
#'
#' @export
#'
#' @examples
#' # --- Example Setup ---
#' # Create a dummy tibble mimicking scRNA-seq average expression data
#' set.seed(123)
#' n_genes <- 10
#' n_samples_pah <- 5
#' n_samples_ctrl <- 5
#' gene_names <- paste0("Gene", 1:n_genes)
#'
#' dummy_data <- tibble::tibble(
#'   gene = rep(gene_names, each = n_samples_pah + n_samples_ctrl),
#'   disease = rep(c(rep("PAH", n_samples_pah), rep("CTRL", n_samples_ctrl)), times = n_genes),
#'   # Simulate some difference for the first few genes
#'   average_expression = rnorm(n_genes * (n_samples_pah + n_samples_ctrl),
#'                              mean = ifelse(gene %in% paste0("Gene", 1:3) & disease == "PAH", 5, 3),
#'                              sd = 1.5)
#' )
#' # Introduce some NAs or edge cases
#' dummy_data$average_expression[1] <- NA # NA value
#' dummy_data$average_expression[dummy_data$gene == "Gene5" & dummy_data$disease == "PAH"] <- 2 # Identical values
#' dummy_data$average_expression[dummy_data$gene == "Gene5" & dummy_data$disease == "CTRL"] <- 2 # Identical values
#' dummy_data <- dplyr::filter(dummy_data, !(gene == "Gene6" & disease == "CTRL")) # Missing group
#'
#' print("Input Dummy Data (head):")
#' print(head(dummy_data))
#'
#' # --- Run the function ---
#' results_table <- test_wilcox_on_average(dummy_data)
#'
#' # --- View Results ---
#' print("Output Results Table:")
#' print(results_table)
#'
#' # Check results for a specific gene where a difference was simulated
#' print("Results for Gene1:")
#' print(dplyr::filter(results_table, gene == "Gene1"))
#'
#' # Check results for Gene5 (where Wilcox test might yield NA p-value due to identical values)
#' print("Results for Gene5:")
#' print(dplyr::filter(results_table, gene == "Gene5"))
#'
#' # Check results for Gene6 (where one group was missing) - should be filtered out
#' print("Results for Gene6 (should be absent):")
#' print(dplyr::filter(results_table, gene == "Gene6"))
#'
test_wilcox_on_average <- function(data) {

  # Input validation (basic checks)
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame or tibble.")
  }
  required_cols <- c("gene", "average_expression", "disease")
  if (!all(required_cols %in% names(data))) {
    stop(sprintf("'data' must contain columns: %s", paste(required_cols, collapse = ", ")))
  }
  if (!("PAH" %in% data$disease && "CTRL" %in% data$disease)) {
    warning("Input data does not contain both 'PAH' and 'CTRL' values in the 'disease' column. Results may be empty or unexpected.", call. = FALSE)
  }

  # Perform calculations using dplyr and tidyr
  # Use .data pronoun for non-standard evaluation safety in packages
  result <- data %>%
    dplyr::group_by(.data$gene) %>%
    dplyr::summarise(
      # Calculate mean expression for each condition, handling potential NAs
      Average_PAH = mean(.data$average_expression[.data$disease == "PAH"], na.rm = TRUE),
      Average_CTRL = mean(.data$average_expression[.data$disease == "CTRL"], na.rm = TRUE),

      # Calculate log2 fold change, handle division by zero or negative means
      # (log2(0/x) = -Inf, log2(x/0) = Inf, log2(0/0) = NaN)
      log2FoldChange = dplyr::case_when(
        Average_PAH > 0 & Average_CTRL > 0 ~ log2(Average_PAH / Average_CTRL),
        Average_PAH > 0 & (Average_CTRL == 0 | is.na(Average_CTRL)) ~ Inf, # Treat 0 or NA denominator as Inf fold change
        (Average_PAH == 0 | is.na(Average_PAH)) & Average_CTRL > 0 ~ -Inf, # Treat 0 or NA numerator as -Inf fold change
        TRUE ~ NaN # Handles 0/0, NA/NA, NA/0, 0/NA etc.
      ),


      # Perform Wilcoxon test
      p_value = {
        # Extract values for each group, removing NAs explicitly
        pah_values <- .data$average_expression[.data$disease == "PAH" & !is.na(.data$average_expression)]
        ctrl_values <- .data$average_expression[.data$disease == "CTRL" & !is.na(.data$average_expression)]

        # Perform test only if both groups have at least one non-NA value
        if (length(pah_values) > 0 && length(ctrl_values) > 0) {
          tryCatch(
            {
              # Check for zero variance within groups (wilcox.test issues warning)
              # Check if all values are identical between groups (wilcox.test gives p=1 or error depending on version/settings)
              if (length(unique(c(pah_values, ctrl_values))) <= 1) {
                # If all values across both groups are identical, p-value is arguably 1 (no difference)
                # or NA (test is ill-defined). Let's return NA for safety/clarity.
                NA_real_
              } else {
                stats::wilcox.test(
                  pah_values,
                  ctrl_values,
                  exact = FALSE # Use approximation for speed
                )$p.value
              }
            },
            # Catch errors during the test (e.g., other edge cases)
            error = function(e) {
              # Optional: Issue a warning specific to the gene
              # warning(sprintf("Wilcox test failed for gene '%s': %s", unique(.data$gene), e$message), call. = FALSE)
              NA_real_
            }
          )
        } else {
          # If one or both groups have no non-NA data, p-value is NA
          NA_real_
        }
      }, # End of p_value calculation block
      .groups = 'drop' # Drop grouping for subsequent operations
    ) %>%
    # Adjust p-values across all genes
    dplyr::mutate(
      # Use p.adjust only on non-NA p-values
      # Note: p.adjust handles NAs gracefully, returning NA for NA inputs
      adjusted_p_value.BH = stats::p.adjust(.data$p_value, method = "BH"),
      adjusted_p_value.FDR = stats::p.adjust(.data$p_value, method = "fdr"), # FDR is alias for BH
      adjusted_p_value.Bonferroni = stats::p.adjust(.data$p_value, method = "bonferroni")
    ) %>%
    # Filter out rows with NA/NaN/Inf results that make interpretation difficult
    dplyr::filter(
      !is.na(.data$log2FoldChange) &
        !is.na(.data$p_value) &
        !is.nan(.data$log2FoldChange) &
        is.finite(.data$log2FoldChange) # Filters out Inf and -Inf log2FC
      # Keep rows with NA p-values if log2FC is finite, as fold change might still be informative
      # Or, add !is.na(p_value) if only complete cases are desired.
      # Current filter keeps rows with valid, finite log2FC and valid (non-NA) p-value.
    )

  return(result)
}

#' @rdname test_wilcox_on_average
#' @export
Test.wilcox.onAverage <- test_wilcox_on_average

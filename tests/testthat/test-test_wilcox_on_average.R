test_that("test_wilcox_on_average works correctly", {
  
  # Create test data
  set.seed(123)
  n_genes <- 5
  n_samples_pah <- 4
  n_samples_ctrl <- 4
  gene_names <- paste0("Gene", 1:n_genes)
  
  test_data <- data.frame(
    gene = rep(gene_names, each = n_samples_pah + n_samples_ctrl),
    disease = rep(c(rep("PAH", n_samples_pah), rep("CTRL", n_samples_ctrl)), times = n_genes),
    average_expression = rnorm(n_genes * (n_samples_pah + n_samples_ctrl), mean = 3, sd = 1.5)
  )
  
  # Make Gene1 have higher expression in PAH
  test_data$average_expression[test_data$gene == "Gene1" & test_data$disease == "PAH"] <- 
    test_data$average_expression[test_data$gene == "Gene1" & test_data$disease == "PAH"] + 2
  
  # Test the function
  results <- test_wilcox_on_average(test_data)
  
  # Basic checks
  expect_s3_class(results, "data.frame")
  expect_true(nrow(results) > 0)
  expect_true(nrow(results) <= n_genes)
  
  # Check required columns
  required_cols <- c("gene", "Average_PAH", "Average_CTRL", "log2FoldChange", 
                     "p_value", "adjusted_p_value.BH", "adjusted_p_value.FDR", 
                     "adjusted_p_value.Bonferroni")
  expect_true(all(required_cols %in% colnames(results)))
  
  # Check that p-values are between 0 and 1 (excluding NAs)
  expect_true(all(results$p_value >= 0 & results$p_value <= 1, na.rm = TRUE))
  expect_true(all(results$adjusted_p_value.BH >= 0 & results$adjusted_p_value.BH <= 1, na.rm = TRUE))
  
  # Check that fold changes are numeric and finite
  expect_true(all(is.finite(results$log2FoldChange)))
  expect_true(all(is.numeric(results$log2FoldChange)))
  
  # Check that Gene1 has positive log2FC (higher in PAH)
  gene1_result <- results[results$gene == "Gene1", ]
  if(nrow(gene1_result) > 0) {
    expect_true(gene1_result$log2FoldChange > 0)
  }
})

test_that("test_wilcox_on_average handles edge cases", {
  
  # Test with missing required columns
  bad_data <- data.frame(x = 1:5, y = 6:10)
  expect_error(test_wilcox_on_average(bad_data))
  
  # Test with non-data.frame input
  expect_error(test_wilcox_on_average(list(a = 1, b = 2)))
  
  # Test with missing disease groups
  partial_data <- data.frame(
    gene = c("Gene1", "Gene1"),
    disease = c("PAH", "PAH"),
    average_expression = c(1, 2)
  )
  expect_warning(test_wilcox_on_average(partial_data))
})

test_that("backward compatibility alias works", {
  
  # Create simple test data
  test_data <- data.frame(
    gene = c("Gene1", "Gene1", "Gene1", "Gene1"),
    disease = c("PAH", "PAH", "CTRL", "CTRL"),
    average_expression = c(5, 6, 3, 4)
  )
  
  # Test that old function name works
  results_new <- test_wilcox_on_average(test_data)
  results_old <- Test.wilcox.onAverage(test_data)
  
  # Should produce identical results
  expect_equal(results_new, results_old)
})

test_that("test_wilcox_on_average handles NA values correctly", {
  
  # Create test data with NAs
  test_data <- data.frame(
    gene = c("Gene1", "Gene1", "Gene1", "Gene1", "Gene2", "Gene2", "Gene2", "Gene2"),
    disease = c("PAH", "PAH", "CTRL", "CTRL", "PAH", "PAH", "CTRL", "CTRL"),
    average_expression = c(5, NA, 3, 4, 7, 8, 5, 6)
  )
  
  results <- test_wilcox_on_average(test_data)
  
  # Should still return results
  expect_s3_class(results, "data.frame")
  expect_true(nrow(results) > 0)
  
  # Check that averages are computed correctly with NA handling
  expect_true(all(is.finite(results$Average_PAH) | is.na(results$Average_PAH)))
  expect_true(all(is.finite(results$Average_CTRL) | is.na(results$Average_CTRL)))
})
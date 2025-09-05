#' Process Tibbles and Calculate Averages
#'
#' This function processes a list of tibbles, calculates averages using a provided
#' test function, and saves the results. It provides progress tracking and verbose
#' output options.
#'
#' @param tibble_list A list of tibbles to process
#' @param user_suffix Character string to append to variable names (default: "APS")
#' @param save_name Base name for the saved results file (default: "Results_list")
#' @param save_path Path where to save the results (optional)
#' @param test_wilcox_fn Function to calculate averages (required)
#' @param verbose Logical; whether to print detailed messages (default: TRUE)
#' @param show_progress Logical; whether to show progress bars (default: TRUE)
#'
#' @return Invisibly returns a list containing the processed results
#'
#' @import progress
#' @importFrom stats setNames
#'
#' @examples
#' \dontrun{
#' results <- process_tibbles_and_calculate(
#'     tibble_list = my_tibble_list,
#'     test_wilcox_fn = my_test_function
#' )
#' }
#'
#' @export
calculate_wilcox_on_tibble_list <- function(tibble_list,
                                          user_suffix = "APS",
                                          save_name = "Results_list",
                                          save_path = NULL,
                                          test_wilcox_fn = NULL,
                                          verbose = TRUE,
                                          show_progress = TRUE) {

  # Input validation
  if (!is.list(tibble_list)) {
    stop("Input must be a list of tibbles")
  }

  if (is.null(test_wilcox_fn)) {
    stop("Please provide the test_wilcox function")
  }

  # Helper function for verbose output
  verbose_msg <- function(...) {
    if (verbose) message(...)
  }

  # Initialize progress bars if enabled
  if (show_progress) {
    pb_overall <- progress::progress_bar$new(
      format = "Overall Progress [:bar] :percent eta: :eta",
      total = 3,
      clear = FALSE,
      width = 80
    )

    pb_tibbles <- progress::progress_bar$new(
      format = "Processing Tibbles [:bar] :percent eta: :eta",
      total = length(tibble_list),
      clear = FALSE,
      width = 80
    )
  }

  # Create environment for storing variables
  assign_env <- parent.frame()

  # Process each tibble and store in environment
  verbose_msg("\nStep 1: Processing tibbles and assigning to environment...")
  assigned_names <- list()

  for (name in names(tibble_list)) {
    base_name <- paste0(name, ".", user_suffix)
    assign(base_name, tibble_list[[name]], envir = assign_env)
    assigned_names$base <- c(assigned_names$base, base_name)

    if (verbose) {
      verbose_msg(sprintf("  - Assigned: %s", base_name))
    }

    if (show_progress) pb_tibbles$tick()
  }

  if (show_progress) pb_overall$tick()

  # Calculate averages
  verbose_msg("\nStep 2: Calculating averages...")
  if (show_progress) {
    pb_averages <- progress::progress_bar$new(
      format = "Calculating Averages [:bar] :percent eta: :eta",
      total = length(tibble_list),
      clear = FALSE,
      width = 80
    )
  }

  for (name in names(tibble_list)) {
    base_name <- paste0(name, ".", user_suffix)
    avg_name <- paste0(base_name, ".Average")

    avg_result <- test_wilcox_fn(tibble_list[[name]])
    assign(avg_name, avg_result, envir = assign_env)
    assigned_names$avg <- c(assigned_names$avg, avg_name)

    if (verbose) {
      verbose_msg(sprintf("  - Calculated average: %s", avg_name))
    }

    if (show_progress) pb_averages$tick()
  }

  if (show_progress) pb_overall$tick()

  # Create results list
  verbose_msg("\nStep 3: Creating and saving results list...")
  results_list <- list()

  if (show_progress) {
    pb_list <- progress::progress_bar$new(
      format = "Creating Results List [:bar] :percent eta: :eta",
      total = length(assigned_names$avg),
      clear = FALSE,
      width = 80
    )
  }

  for (avg_name in assigned_names$avg) {
    results_list[[avg_name]] <- get(avg_name, envir = assign_env)

    if (verbose) {
      verbose_msg(sprintf("  - Added to results list: %s", avg_name))
    }

    if (show_progress) pb_list$tick()
  }

  # Save results if path provided
  if (!is.null(save_path)) {
    full_path <- file.path(save_path, paste0(save_name, ".Rds"))
    saveRDS(results_list, file = full_path, compress = TRUE)
    verbose_msg(sprintf("\nResults saved to: %s", full_path))
  }

  if (show_progress) pb_overall$tick()

  # Print summary if verbose
  if (verbose) {
    verbose_msg("\nSummary:")
    verbose_msg(sprintf("Number of tibbles processed: %d", length(tibble_list)))
    verbose_msg(sprintf("Variables created in environment: %d", length(unlist(assigned_names))))
    verbose_msg("\nCreated variables:")
    verbose_msg(paste("  ", unlist(assigned_names), collapse = "\n"))

    if (!is.null(save_path)) {
      verbose_msg(sprintf("\nResults saved to: %s", full_path))
    }
  }

  # Return results list invisibly
  invisible(results_list)
}

#' Process a List of Tibbles, Apply a Function, and Optionally Save Results
#'
#' @description
#' This function iterates through a named list of tibbles (or data frames). For each item,
#' it assigns the original tibble and the result of applying a user-provided function
#' (`test_wilcox_fn`) to the calling environment. It collects the results from
#' `test_wilcox_fn` into a list, which is returned invisibly and can optionally be
#' saved to an RDS file. By default, it uses the exported `Test.wilcox.onAverage` function
#' from this package for the analysis step.
#'
#' @details
#' **Important Note on Side Effects:** This function assigns variables directly into the
#' calling environment (`parent.frame()`) using `assign()`. This includes the original
#' tibbles (with `user_suffix` appended to their names) and the results of the
#' `test_wilcox_fn` (with `.Average` further appended). Assigning to the parent
#' environment is generally discouraged in R package development. Consider modifying
#' the function to return all necessary outputs within a structured list.
#'
#' The function performs these main steps:
#' 1. Validates inputs.
#' 2. Iterates through `tibble_list`, assigning each tibble to `parent.frame()`
#'    with the name `list_element_name.user_suffix`.
#' 3. Iterates again, applying `test_wilcox_fn` (defaults to `Test.wilcox.onAverage`)
#'    to each tibble. Assigns the result to `parent.frame()` with the name
#'    `list_element_name.user_suffix.Average`. Includes basic error handling.
#' 4. Collects all successfully computed results into a list.
#' 5. If `save_path` is provided, saves the results list as an RDS file.
#' 6. Optionally displays progress bars and verbose messages.
#'
#' **Default Function:** The default function `Test.wilcox.onAverage` (exported by this
#' package) expects input tibbles/data frames with columns: `gene`, `average_expression`,
#' and `disease` (containing exactly two levels, e.g., "PAH" and "CTRL").
#'
#' @param tibble_list A named list where each element is a tibble or data frame.
#'        See details for requirements if using the default `test_wilcox_fn`.
#' @param user_suffix A character string suffix. Defaults to "APS".
#' @param save_name A character string for the saved RDS file name. Defaults to "Results_list".
#' @param save_path An optional character string directory path for saving. Defaults to NULL.
#' @param test_wilcox_fn A function to apply to each tibble. Defaults to
#'        `Test.wilcox.onAverage`. Must accept a single tibble/data frame as input.
#' @param verbose Logical. If `TRUE` (default), print progress messages.
#' @param show_progress Logical. If `TRUE` (default), display progress bars. Requires 'progress'.
#'
#' @return Invisibly returns a named list of successfully computed results.
#'
#' @importFrom progress progress_bar
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # (Examples remain the same as in the previous artifact, they demonstrate
#' # calling with the default and with a custom function)
#'
#' # --- Example Setup ---
#' set.seed(123)
#' make_dummy_tibble <- function(name_prefix, n_genes = 5, n_pah = 3, n_ctrl = 3) {
#'   gene_names <- paste0(name_prefix, "_Gene", 1:n_genes)
#'   tibble::tibble(
#'     gene = rep(gene_names, each = n_pah + n_ctrl),
#'     disease = rep(c(rep("PAH", n_pah), rep("CTRL", n_ctrl)), times = n_genes),
#'     average_expression = rnorm(n_genes * (n_pah + n_ctrl), mean = 5, sd = 1.5) +
#'                          ifelse(disease == "PAH", rnorm(n_genes * (n_pah + n_ctrl), 1, 0.5), 0)
#'   )
#' }
#' tibble1 <- make_dummy_tibble("SetA")
#' tibble2 <- make_dummy_tibble("SetB")
#' my_analysis_list <- list(Group1 = tibble1, Group2 = tibble2)
#'
#' # --- Usage with Default Function ---
#' rm(list = ls(pattern = "\\.DefaultTest|\\.DefaultTest\\.Average"))
#' results_default <- Process_tibbles_Calculate_average_wilcox(
#'   tibble_list = my_analysis_list,
#'   user_suffix = "DefaultTest",
#'   verbose = TRUE, show_progress = FALSE
#' )
#' print(results_default)
#' print(ls(pattern = "\\.DefaultTest"))
#' print(head(Group1.DefaultTest.Average))
#' rm(list = ls(pattern = "\\.DefaultTest|\\.DefaultTest\\.Average"))
#'
#' # --- Usage with Custom Function ---
#' simple_row_counter <- function(tb) { nrow(tb) }
#' rm(list = ls(pattern = "\\.CustomTest|\\.CustomTest\\.Average"))
#' results_custom <- Process_tibbles_Calculate_average_wilcox(
#'   tibble_list = my_analysis_list,
#'   user_suffix = "CustomTest",
#'   test_wilcox_fn = simple_row_counter,
#'   verbose = TRUE, show_progress = FALSE
#' )
#' print(results_custom)
#' print(ls(pattern = "\\.CustomTest"))
#' print(Group1.CustomTest.Average)
#' rm(list = ls(pattern = "\\.CustomTest|\\.CustomTest\\.Average"))
#'
#' # --- Usage with Saving (using default function) ---
#' temp_dir <- tempdir()
#' results_saved <- Process_tibbles_Calculate_average_wilcox(
#'   tibble_list = my_analysis_list,
#'   user_suffix = "SaveRun",
#'   save_name = "DefaultWilcoxResults",
#'   save_path = temp_dir,
#'   verbose = FALSE
#' )
#' saved_file <- file.path(temp_dir, "DefaultWilcoxResults.Rds")
#' print(list.files(temp_dir, pattern = "DefaultWilcoxResults\\.Rds"))
#' saved_data <- readRDS(saved_file)
#' print(saved_data)
#' rm(list = ls(pattern = "\\.SaveRun"))
#' file.remove(saved_file)
#' }
#'
Process_tibbles_Calculate_average_wilcox <- function(tibble_list,
                                                     user_suffix = "APS",
                                                     save_name = "Results_list",
                                                     save_path = NULL,
                                                     test_wilcox_fn = Test.wilcox.onAverage, # Set default to the exported function
                                                     verbose = TRUE,
                                                     show_progress = TRUE) {

  # --- Input Validation ---
  # (Validation code remains the same as in the previous artifact)
  if (!is.list(tibble_list)) {
    stop("Input 'tibble_list' must be a list.")
  }
  if (is.null(names(tibble_list)) || any(names(tibble_list) == "")) {
    stop("Input 'tibble_list' must be a named list with non-empty names.")
  }
  if (!all(sapply(tibble_list, function(x) inherits(x, "data.frame")))) {
    stop("All elements in 'tibble_list' must be tibbles or data frames.")
  }
  if (!is.character(user_suffix) || length(user_suffix) != 1 || nchar(user_suffix) == 0) {
    stop("'user_suffix' must be a single, non-empty character string.")
  }
  if (!is.character(save_name) || length(save_name) != 1 || nchar(save_name) == 0) {
    stop("'save_name' must be a single, non-empty character string.")
  }
  if (!is.null(save_path)) {
    if (!is.character(save_path) || length(save_path) != 1) {
      stop("'save_path' must be a single character string or NULL.")
    }
    if (!dir.exists(save_path)) {
      stop("Specified 'save_path' directory does not exist: ", save_path)
    }
  }
  if (!is.function(test_wilcox_fn)) {
    stop("'test_wilcox_fn' must be a function.")
  }
  if (!is.logical(verbose) || length(verbose) != 1) {
    stop("'verbose' must be a single logical value (TRUE or FALSE).")
  }
  if (!is.logical(show_progress) || length(show_progress) != 1) {
    stop("'show_progress' must be a single logical value (TRUE or FALSE).")
  }

  # --- Helper function for verbose output ---
  verbose_msg <- function(...) {
    if (verbose) message(...)
  }

  # --- Initialize progress bars if enabled ---
  # (Code remains the same as in the previous artifact)
  if (show_progress) {
    if (!requireNamespace("progress", quietly = TRUE)) {
      warning("Package 'progress' needed for progress bars but is not installed. Disabling progress bars.", call. = FALSE)
      show_progress <- FALSE
    } else {
      pb_overall <- progress::progress_bar$new(
        format = "Overall Progress [:bar] :percent eta: :eta",
        total = 3,
        clear = FALSE, width = 60
      )
    }
  }
  if (!show_progress) {
    pb_overall <- list(tick = function(...) {})
    pb_assign <- list(tick = function(...) {})
    pb_calc <- list(tick = function(...) {})
  }


  # --- Environment for storing variables ---
  assign_env <- parent.frame()
  verbose_msg("Assigning variables to environment: ", environmentName(assign_env))

  # --- Process each tibble and assign to environment ---
  # (Code remains the same as in the previous artifact)
  verbose_msg("\nStep 1: Assigning original tibbles to environment...")
  assigned_names <- list(base = character(), avg = character())

  if (show_progress) {
    pb_assign <- progress::progress_bar$new(
      format = "Assigning Tibbles [:bar] :percent eta: :eta",
      total = length(tibble_list), clear = FALSE, width = 60
    )
  }

  for (name in names(tibble_list)) {
    base_name <- paste0(name, ".", user_suffix)
    tryCatch({
      assign(base_name, tibble_list[[name]], envir = assign_env)
      assigned_names$base <- c(assigned_names$base, base_name)
      if (verbose) {
        verbose_msg(sprintf("  - Assigned: %s", base_name))
      }
    }, error = function(e) {
      warning(sprintf("Failed to assign base variable '%s': %s", base_name, e$message), call. = FALSE)
    })
    if (show_progress) pb_assign$tick()
  }

  if (show_progress) pb_overall$tick()


  # --- Apply test_wilcox_fn and assign results ---
  # (Code remains the same as in the previous artifact)
  verbose_msg("\nStep 2: Applying 'test_wilcox_fn' and assigning results...")
  results_list_internal <- list() # Store results internally first

  if (show_progress) {
    pb_calc <- progress::progress_bar$new(
      format = "Calculating Results [:bar] :percent eta: :eta",
      total = length(tibble_list), clear = FALSE, width = 60
    )
  }

  for (name in names(tibble_list)) {
    base_name <- paste0(name, ".", user_suffix) # Reconstruct base name
    avg_name <- paste0(base_name, ".Average")

    avg_result <- tryCatch({
      test_wilcox_fn(tibble_list[[name]])
    }, error = function(e) {
      warning(sprintf("Error applying 'test_wilcox_fn' to tibble '%s': %s", name, e$message), call. = FALSE)
      NULL
    })

    if (!is.null(avg_result)) {
      tryCatch({
        assign(avg_name, avg_result, envir = assign_env)
        assigned_names$avg <- c(assigned_names$avg, avg_name)
        results_list_internal[[avg_name]] <- avg_result
        if (verbose) {
          verbose_msg(sprintf("  - Calculated and assigned result: %s", avg_name))
        }
      }, error = function(e) {
        warning(sprintf("Failed to assign result variable '%s': %s", avg_name, e$message), call. = FALSE)
        results_list_internal[[avg_name]] <- NULL
      })
    } else {
      if (verbose) {
        verbose_msg(sprintf("  - Failed to calculate result for: %s (Skipping assignment)", name))
      }
    }
    if (show_progress) pb_calc$tick()
  }

  if (show_progress) pb_overall$tick()

  # --- Create final results list and optionally save ---
  # (Code remains the same as in the previous artifact)
  verbose_msg("\nStep 3: Finalizing results list and optionally saving...")
  final_results_list <- results_list_internal[!sapply(results_list_internal, is.null)]

  full_path <- NULL
  if (!is.null(save_path)) {
    full_path <- file.path(save_path, paste0(save_name, ".Rds"))
    tryCatch({
      saveRDS(final_results_list, file = full_path, compress = TRUE)
      verbose_msg(sprintf("\nResults list saved successfully to: %s", full_path))
    }, error = function(e) {
      warning(sprintf("Failed to save results list to '%s': %s", full_path, e$message), call. = FALSE)
      full_path <- NULL
    })
  } else {
    verbose_msg("\nResults list not saved (save_path was NULL).")
  }

  if (show_progress) pb_overall$tick()

  # --- Print summary if verbose ---
  # (Code remains the same as in the previous artifact)
  if (verbose) {
    verbose_msg("\n--- Processing Summary ---")
    verbose_msg(sprintf("Number of input tibbles: %d", length(tibble_list)))
    verbose_msg(sprintf("Variables assigned to environment '%s': %d (Tibbles: %d, Results: %d)",
                        environmentName(assign_env),
                        length(assigned_names$base) + length(assigned_names$avg),
                        length(assigned_names$base),
                        length(assigned_names$avg)))
    if (length(assigned_names$base) + length(assigned_names$avg) > 0) {
      verbose_msg("\nAssigned variable names:")
      if(length(assigned_names$base) > 0) verbose_msg(paste("  Tibbles:", paste(assigned_names$base, collapse = ", ")))
      if(length(assigned_names$avg) > 0) verbose_msg(paste("  Results:", paste(assigned_names$avg, collapse = ", ")))
    }
    verbose_msg(sprintf("\nNumber of results successfully computed: %d", length(final_results_list)))

    if (!is.null(save_path) && !is.null(full_path)) {
      verbose_msg(sprintf("Results list saved to: %s", full_path))
    } else if (!is.null(save_path) && is.null(full_path)) {
      verbose_msg("Attempted to save results list, but failed.")
    }
    verbose_msg("--- End of Summary ---")
  }

  # --- Return results list invisibly ---
  invisible(final_results_list)
}

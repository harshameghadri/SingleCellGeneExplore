#' Perform Comprehensive Enrichment Analysis and Pathway Visualization
#'
#' This function conducts a series of enrichment analyses (Enrichr, GO, KEGG, Reactome)
#' for up- and down-regulated genes from differential expression results for multiple
#' cell types. It also integrates Pathview for KEGG pathway visualization and
#' generates various plots (barplots, cnetplots).
#'
#' @param de_results_list A named list of data frames. Each data frame should
#'   contain differential expression results for a specific cell type and must
#'   include columns: "gene" (gene symbol), "p_val_adj" (adjusted p-value),
#'   and "avg_log2FC" (average log2 fold change).
#' @param top_n An integer specifying the maximum number of enriched terms or
#'   pathways to display in plots (e.g., barplots, cnetplots). Default: 20.
#' @param pval_cutoff A numeric value for the adjusted p-value cutoff to
#'   determine statistical significance for enrichment results display and gene selection.
#'   Default: 0.05.
#' @param up_log2fc_threshold A numeric value for the minimum positive log2 fold
#'   change to consider a gene as up-regulated. Default: 0.25.
#' @param down_log2fc_threshold A numeric value for the maximum negative log2 fold
#'   change to consider a gene as down-regulated. Default: -0.25.
#' @param max_genes_for_enrichment An optional numeric value. If provided and positive,
#'   it caps the number of most significant up-regulated genes and, separately,
#'   down-regulated genes (sorted by p-value then absolute fold change) submitted
#'   for enrichment. If NULL, Inf, or non-positive, all genes passing
#'   `pval_cutoff` and log2FC thresholds are used. Default: NULL.
#' @param enrichr_dbs A character vector of Enrichr database names to query.
#'   Default: c("GO_Biological_Process_2023", "GO_Cellular_Component_2023", ...).
#' @param outdir A string specifying the main output directory path where all
#'   results will be saved. Default: "enrichment_results".
#' @param kegg_n An integer specifying the maximum number of top KEGG pathways
#'   (ordered by p-value) to visualize using Pathview. Default: 10.
#' @param use_full_fc_for_pathview A logical value. If TRUE, Pathview visualizations
#'   will use fold change data from all significant genes in the input DE table for
#'   coloring, rather than just the subset currently being analyzed (e.g., up-regulated only).
#'   Default: TRUE.
#' @param pathview_colors A character vector of three colors for Pathview:
#'   low expression, mid-expression (neutral), and high expression.
#'   Default: c("blue", "gray", "red").
#' @param save_all_results A logical value. If TRUE, all enrichment tables
#'   (e.g., from clusterProfiler) are saved with a lenient p-value cutoff (1.0),
#'   allowing inspection of all terms. Plots will still use `pval_cutoff`.
#'   If FALSE, tables are also filtered by `pval_cutoff`. Default: TRUE.
#' @param plot_width A numeric value for the width of individual saved plots.
#'   Default: 8.
#' @param plot_height A numeric value for the height of individual saved plots.
#'   Default: 6.
#' @param plot_units A string specifying the units for `plot_width` and
#'   `plot_height` (e.g., "in", "cm", "mm"). Default: "in".
#' @param plot_dpi A numeric value for the resolution (DPI) of individual saved plots.
#'   Default: 300.
#' @param plot_formats A character vector of formats for saving individual plots
#'   (e.g., "png", "pdf", "svg"). Default: c("png").
#' @param pdf_width A numeric value for the width of summary multi-page PDF plots.
#'   Default: 11.
#' @param pdf_height A numeric value for the height of summary multi-page PDF plots.
#'   Default: 8.5.
#' @param cnet_color_edge A logical value. If TRUE, edges in cnetplots will be
#'   colored based on fold change. Default: TRUE.
#' @param parallel.per.celltype A logical value. If TRUE, processing for each
#'   cell type in `de_results_list` will be attempted in parallel using the
#'   `future` framework. Default: FALSE.
#' @param n_cores An optional integer specifying the number of cores to use for
#'   parallel processing. If NULL, `future::availableCores() - 1` will be used.
#'   Only relevant if `parallel.per.celltype` is TRUE. Default: NULL.
#'
#' @return Invisibly returns a list named by cell types, where each element
#'   contains sub-lists for "UpRegulated_CP_results" and "DownRegulated_CP_results"
#'   (the results from `run_clusterprofiler_and_pathview`).
#'
#' @importFrom dplyr filter arrange distinct pull mutate group_by summarise inner_join select bind_rows
#' @importFrom tibble tibble
#' @importFrom enrichR listEnrichrDbs enrichr plotEnrich setEnrichrSite
#' @importFrom clusterProfiler bitr enrichGO enrichKEGG enrichPathway setReadable
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom ReactomePA annFUN.org
#' @importFrom enrichplot barplot cnetplot
#' @importFrom ggplot2 ggplot ggsave ggtitle theme element_text
#' @importFrom pathview pathview
#' @importFrom KEGGREST keggGet
#' @importFrom future plan multisession sequential availableCores
#' @importFrom future.apply future_lapply
#'
#' @examples
#' \dontrun{
#' # Create dummy DE results list
#' set.seed(123)
#' gene_names <- paste0("Gene", 1:1000)
#' de_data1 <- data.frame(
#'   gene = sample(gene_names, 500),
#'   avg_log2FC = rnorm(500, 0, 1.5),
#'   p_val_adj = runif(500, 0, 0.2)
#' )
#' de_data1$avg_log2FC[1:20] <- rnorm(20, 2, 0.5); de_data1$p_val_adj[1:20] <- runif(20, 0, 0.04)
#' de_data1$avg_log2FC[21:40] <- rnorm(20, -2, 0.5); de_data1$p_val_adj[21:40] <- runif(20, 0, 0.04)
#'
#' de_data2 <- data.frame(
#'   gene = sample(gene_names, 600),
#'   avg_log2FC = rnorm(600, 0, 1),
#'   p_val_adj = runif(600, 0, 0.3)
#' )
#' de_data2$avg_log2FC[1:15] <- rnorm(15, 1.8, 0.4); de_data2$p_val_adj[1:15] <- runif(15, 0, 0.03)
#' de_data2$avg_log2FC[16:35] <- rnorm(20, -1.5, 0.6); de_data2$p_val_adj[16:35] <- runif(20, 0, 0.045)
#'
#' my_de_list <- list(Epithelial_Type1 = de_data1, Stromal_TypeA = de_data2)
#'
#' # Run the analysis
#' results <- perform_enrichment_dual_pathview_refactored(
#'   de_results_list = my_de_list,
#'   outdir = "my_enrichment_output_roxygen_example",
#'   pval_cutoff = 0.05,
#'   up_log2fc_threshold = 0.5,
#'   down_log2fc_threshold = -0.5,
#'   max_genes_for_enrichment = 150, # Analyze top 150 up/down
#'   kegg_n = 5,
#'   plot_formats = c("png", "pdf"),
#'   parallel.per.celltype = TRUE,
#'   n_cores = 2
#' )
#'
#' # To use all significant genes (default behavior for max_genes_for_enrichment):
#' results_all_sig <- perform_enrichment_dual_pathview_refactored(
#'   de_results_list = my_de_list,
#'   outdir = "my_enrichment_output_roxygen_all_sig",
#'   max_genes_for_enrichment = NULL # or Inf, or omit the parameter
#' )
#' }
#'
perform_enrichment_dual_pathview_refactored <- function(
    de_results_list,
    top_n = 20,
    pval_cutoff = 0.05,
    up_log2fc_threshold = 0.25,
    down_log2fc_threshold = -0.25,
    max_genes_for_enrichment = NULL,
    enrichr_dbs = c(
      "GO_Biological_Process_2023",
      "GO_Cellular_Component_2023",
      "GO_Molecular_Function_2023",
      "WikiPathways_2024_Human", # Note: DB names can change, ensure these are current
      "Reactome_Pathways_2024", # Or specific Reactome e.g., "Reactome_Pathways_Human_2022"
      "KEGG_2021_Human",
      "ChEA_2022",
      "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X",
      "TRANSFAC_and_JASPAR_PWMs",
      "MSigDB_Hallmark_2020"
    ),
    outdir = "enrichment_results",
    kegg_n = 10,
    use_full_fc_for_pathview = TRUE,
    pathview_colors = c("blue", "gray", "red"),
    save_all_results = TRUE,
    plot_width = 8,
    plot_height = 6,
    plot_units = "in",
    plot_dpi = 300,
    plot_formats = c("png"),
    pdf_width = 11,
    pdf_height = 8.5,
    cnet_color_edge = TRUE,
    parallel.per.celltype = FALSE,
    n_cores = NULL
) {
  # Load required packages quietly
  suppressPackageStartupMessages({
    library(dplyr)
    library(tibble)
    library(enrichR)
    library(clusterProfiler)
    library(org.Hs.eg.db) # Ensure this is installed
    library(ReactomePA)
    library(enrichplot)
    library(ggplot2)
    library(pathview)
    if (!requireNamespace("KEGGREST", quietly = TRUE)) {
      message("Installing KEGGREST package as it's not found...")
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install("KEGGREST", update = FALSE)
    }
    library(KEGGREST)
    if (!requireNamespace("future", quietly = TRUE)) {
      message("Installing future package as it's not found...")
      install.packages("future")
    }
    library(future)
    if (!requireNamespace("future.apply", quietly = TRUE)) {
      message("Installing future.apply package as it's not found...")
      install.packages("future.apply")
    }
    library(future.apply)
  })
  
  # --- Input Validation ---
  if (!is.list(de_results_list) || length(de_results_list) == 0) {
    stop("`de_results_list` must be a non-empty list of data frames.")
  }
  if (!all(sapply(de_results_list, is.data.frame))) {
    stop("Each element in `de_results_list` must be a data frame.")
  }
  required_cols <- c("gene", "p_val_adj", "avg_log2FC")
  for (ct_name in names(de_results_list)) {
    if (!all(required_cols %in% colnames(de_results_list[[ct_name]]))) {
      stop(paste0("Data frame for cell type '", ct_name,
                  "' is missing one or more required columns: 'gene', 'p_val_adj', 'avg_log2FC'."))
    }
  }
  if (!is.numeric(top_n) || top_n <= 0) stop("`top_n` must be a positive number.")
  if (!is.numeric(pval_cutoff) || pval_cutoff < 0 || pval_cutoff > 1) stop("`pval_cutoff` must be between 0 and 1.")
  if (!is.numeric(up_log2fc_threshold)) stop("`up_log2fc_threshold` must be numeric.")
  if (!is.numeric(down_log2fc_threshold)) stop("`down_log2fc_threshold` must be numeric.")
  
  if (!is.null(max_genes_for_enrichment)) {
    if (!is.numeric(max_genes_for_enrichment) ||
        (max_genes_for_enrichment <= 0 && !is.infinite(max_genes_for_enrichment))) {
      warning(paste0("`max_genes_for_enrichment` (currently: ", max_genes_for_enrichment, ") must be NULL, Inf, or a positive number. ",
                     "Defaulting to NULL (all significant genes)."))
      max_genes_for_enrichment <- NULL # Default to using all significant if invalid
    }
  }
  
  if (!is.character(enrichr_dbs) || length(enrichr_dbs) == 0) stop("`enrichr_dbs` must be a non-empty character vector.")
  if (!is.character(outdir) || length(outdir) != 1) stop("`outdir` must be a single string.")
  if (!is.numeric(kegg_n) || kegg_n <= 0) stop("`kegg_n` must be a positive number.")
  if (!is.logical(use_full_fc_for_pathview) || length(use_full_fc_for_pathview) != 1) stop("`use_full_fc_for_pathview` must be a single logical value.")
  if (!is.character(pathview_colors) || length(pathview_colors) != 3) stop("`pathview_colors` must be a character vector of length 3.")
  if (!is.logical(save_all_results) || length(save_all_results) != 1) stop("`save_all_results` must be a single logical value.")
  if (!is.numeric(plot_width) || plot_width <= 0) stop("`plot_width` must be a positive number.")
  if (!is.numeric(plot_height) || plot_height <= 0) stop("`plot_height` must be a positive number.")
  if (!plot_units %in% c("in", "cm", "mm")) stop("`plot_units` must be 'in', 'cm', or 'mm'.")
  if (!is.numeric(plot_dpi) || plot_dpi <= 0) stop("`plot_dpi` must be a positive number.")
  if (!is.character(plot_formats) || length(plot_formats) == 0 || !all(sapply(plot_formats, is.character))) stop("`plot_formats` must be a non-empty character vector.")
  if (!is.numeric(pdf_width) || pdf_width <= 0) stop("`pdf_width` must be a positive number.")
  if (!is.numeric(pdf_height) || pdf_height <= 0) stop("`pdf_height` must be a positive number.")
  if (!is.logical(cnet_color_edge) || length(cnet_color_edge) != 1) stop("`cnet_color_edge` must be a single logical value.")
  if (!is.logical(parallel.per.celltype) || length(parallel.per.celltype) != 1) stop("`parallel.per.celltype` must be a single logical value.")
  if (!is.null(n_cores) && (!is.numeric(n_cores) || n_cores <= 0 || n_cores %% 1 != 0)) stop("`n_cores` must be NULL or a positive integer.")
  
  
  # --- Helper Functions ---
  create_dir_if_needed <- function(dir_path) {
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
    }
  }
  
  run_enrichment <- function(fun, ...) {
    tryCatch({
      result <- fun(...)
      if (is.null(result) || !inherits(result, "enrichResult") || nrow(result@result) == 0) {
        return(NULL)
      }
      result
    }, error = function(e) {
      warning(paste0("clusterProfiler enrichment (", deparse(substitute(fun)), ") failed: ", conditionMessage(e)))
      return(NULL)
    })
  }
  
  clean_pathway_name <- function(pathway_name, kegg_id) {
    if (is.null(pathway_name) || pathway_name == "") {
      pathway_name_clean <- paste0("Pathway_", kegg_id)
    } else {
      pathway_name_clean <- gsub("[\\/:*?\"<>| ]+", "_", pathway_name) # Replace problematic characters with underscore
      pathway_name_clean <- gsub("[^[:alnum:]_.-]", "", pathway_name_clean) # Remove any remaining non-alphanumeric (except _, ., -)
      pathway_name_clean <- gsub("_+", "_", pathway_name_clean) # Replace multiple underscores with single
      pathway_name_clean <- gsub("^_|_$", "", pathway_name_clean) # Trim leading/trailing underscores
      pathway_name_clean <- substr(pathway_name_clean, 1, 100) # Truncate if too long
      if (pathway_name_clean == "") { # Fallback if cleaning results in empty string
        pathway_name_clean <- paste0("Pathway_", kegg_id)
      }
    }
    return(pathway_name_clean)
  }
  
  run_and_rename_pathview <- function(kegg_id, kegg_desc, fc_vector, out_suffix_pv, target_dir, local_pathway_name_cache) {
    pathview_generated <- FALSE
    current_wd <- getwd()
    on.exit(setwd(current_wd), add = TRUE) # Ensure WD is reset
    
    tryCatch({
      create_dir_if_needed(target_dir)
      setwd(target_dir) # pathview writes to current WD
      
      fc_limit <- max(abs(fc_vector), na.rm = TRUE)
      if(!is.finite(fc_limit) || fc_limit == 0) fc_limit <- 1 # Default limit if NA or 0
      
      # Actual pathview call
      pv_obj <- pathview::pathview(
        gene.data  = fc_vector,
        pathway.id = kegg_id,
        species    = "hsa",
        limit      = list(gene = fc_limit, cpd = 1), # cpd limit is often 1
        low        = list(gene = pathview_colors[1], cpd = "blue"), # Colors for gene and compound data
        mid        = list(gene = pathview_colors[2], cpd = "gray"),
        high       = list(gene = pathview_colors[3], cpd = "red"),
        node.sum   = "mean", # Method to summarize multiple gene products mapping to the same node
        out.suffix = out_suffix_pv, # Suffix for output files
        kegg.native = TRUE, # Generate KEGG native graph image
        same.layer = TRUE, # Try to plot gene and compound data on the same layer
        sign.pos = NA, # Position of the signature, NA for default
        verbose = FALSE # Suppress verbose messages from pathview
      )
      
      # File renaming logic (pathview output names can be tricky)
      # pathview with kegg.native=TRUE typically creates pathway.id.suffix.png
      # but sometimes just pathway.id.png if suffix is simple.
      # Let's try to find the generated PNG.
      expected_png_kegg_native <- paste0(kegg_id, ".png") # pathview now often creates this directly
      expected_png_with_suffix <- paste0(kegg_id, ".", out_suffix_pv, ".png")
      expected_xml <- paste0(kegg_id, ".xml") # XML file with graph data
      
      png_file_found <- NULL
      if (file.exists(expected_png_kegg_native)) {
        png_file_found <- expected_png_kegg_native
      } else if (file.exists(expected_png_with_suffix)) {
        png_file_found <- expected_png_with_suffix
      } else {
        # Fallback: search for any PNG starting with the KEGG ID in the current (target_dir)
        potential_files <- list.files(pattern = paste0("^", kegg_id, ".*\\.png$"), ignore.case = TRUE)
        if(length(potential_files) > 0) png_file_found <- potential_files[1]
      }
      
      
      if (!is.null(png_file_found)) {
        pathview_generated <- TRUE
        pathway_name <- local_pathway_name_cache[[kegg_id]] # Check cache first
        
        if (is.null(pathway_name)) {
          pathway_info <- tryCatch({
            KEGGREST::keggGet(kegg_id)
          }, error = function(e) {
            warning(paste("KEGGREST::keggGet failed for", kegg_id, ":", conditionMessage(e)))
            NULL
          })
          if (!is.null(pathway_info) && length(pathway_info) > 0 && "NAME" %in% names(pathway_info[[1]])) {
            pathway_name <- pathway_info[[1]]$NAME
            local_pathway_name_cache[[kegg_id]] <- pathway_name # Update cache
          } else {
            pathway_name <- kegg_desc # Fallback to description from enrichResult
          }
        }
        
        pathway_name_clean <- clean_pathway_name(pathway_name, kegg_id)
        new_png_name <- paste0(out_suffix_pv, "_", kegg_id, "_", pathway_name_clean, ".png")
        old_png_path <- png_file_found # This is already in target_dir because we did setwd(target_dir)
        
        tryCatch({
          file.rename(old_png_path, new_png_name)
        }, warning = function(w){
          warning(paste("Could not rename pathview PNG '", old_png_path, "' to '", new_png_name, "': ", w$message))
        }, error = function(e){
          warning(paste("Error renaming pathview PNG '", old_png_path, "' to '", new_png_name, "': ", e$message))
        })
        
        # Rename XML file if it exists
        if (file.exists(expected_xml)) {
          new_xml_name <- paste0(out_suffix_pv, "_", kegg_id, "_", pathway_name_clean, ".xml")
          tryCatch({
            file.rename(expected_xml, new_xml_name)
          }, warning = function(w){}, error = function(e){}) # Silently try to rename XML
        }
        
      } else {
        warning(paste("Pathview PNG file not found for pathway ID:", kegg_id, "with suffix:", out_suffix_pv, "in directory:", getwd()))
      }
    }, error = function(e) {
      warning(paste("Pathview visualization failed for pathway ID", kegg_id, ":", conditionMessage(e)))
    }, finally = {
      setwd(current_wd) # Crucial to reset WD
    })
    return(list(generated = pathview_generated, cache = local_pathway_name_cache))
  }
  
  
  # --- Main Setup ---
  create_dir_if_needed(outdir)
  csv_dir_base <- file.path(outdir, "csv_results") # Base CSV directory
  create_dir_if_needed(csv_dir_base)
  
  
  ## ---- Helper Function: enrichR Analysis ----
  run_enrichr_analysis <- function(gene_list, label, celltype, enrichr_csv_dir_celltype) {
    # Define output directories for this specific celltype and direction
    celltype_specific_outdir <- file.path(outdir, celltype) # e.g., enrichment_results/CellTypeA
    enrichr_main_dir_ct <- file.path(celltype_specific_outdir, "enrichR") # e.g., enrichment_results/CellTypeA/enrichR
    direction_specific_dir <- file.path(enrichr_main_dir_ct, label) # e.g., enrichment_results/CellTypeA/enrichR/UpRegulated
    plots_individual_dir <- file.path(direction_specific_dir, "plots_individual")
    
    create_dir_if_needed(direction_specific_dir)
    create_dir_if_needed(plots_individual_dir)
    
    enrichr_results_list <- list()
    # Try to set Enrichr site (e.g., "Enrichr", "FlyEnrichr", "WormEnrichr")
    # This might be needed if the default site is down or for specific organisms
    tryCatch(enrichR::setEnrichrSite("Enrichr"), error = function(e) message("Could not set Enrichr site, using default."))
    
    available_dbs_df <- tryCatch(enrichR::listEnrichrDbs(), error = function(e) {
      warning("Could not retrieve Enrichr DB list: ", e$message); NULL
    })
    
    if(is.null(available_dbs_df) || !"libraryName" %in% colnames(available_dbs_df)) {
      warning("Skipping Enrichr for ", celltype, " ", label, " as DB list could not be fetched.")
      return()
    }
    available_dbs_vector <- available_dbs_df$libraryName
    
    for (db_name in enrichr_dbs) {
      if (!db_name %in% available_dbs_vector) {
        warning(paste0(db_name, " is not a recognized EnrichR database. Skipping... ",
                       "(Available DBs start with: ", paste(head(available_dbs_vector), collapse=", "), ")"))
        next
      }
      enriched_data <- tryCatch({
        enrichr(gene_list, databases = db_name) # Returns a list with one element named db_name
      }, error = function(e) {
        warning(paste0("enrichR query failed for DB: '", db_name, "' in ", celltype, " - ", label, ": ", e$message))
        NULL
      })
      
      if (!is.null(enriched_data) && length(enriched_data) > 0 && db_name %in% names(enriched_data)) {
        db_table <- enriched_data[[db_name]]
        if (!is.null(db_table) && nrow(db_table) > 0) {
          enrichr_results_list[[db_name]] <- db_table
        } else {
          # message(paste0("No results from Enrichr for DB: '", db_name, "' in ", celltype, " - ", label))
        }
      }
    }
    
    if (length(enrichr_results_list) == 0) {
      message(paste0("No results obtained from any Enrichr databases for ", celltype, " - ", label))
      return()
    }
    
    # Combine results and save to CSV
    combined_enrichr_df <- dplyr::bind_rows(
      lapply(names(enrichr_results_list), function(db_name) {
        df <- enrichr_results_list[[db_name]]
        # Ensure essential columns are present before mutating
        if(all(c("Term", "Overlap", "Adjusted.P.value") %in% colnames(df))) {
          df <- df %>%
            dplyr::mutate(Database = db_name,
                          Count = suppressWarnings(as.numeric(sub("/.*", "", Overlap))))
          return(df)
        } else {
          warning(paste0("Skipping formatting for Enrichr DB '", db_name, "' due to missing essential columns."))
          return(NULL) # Return NULL if columns are missing
        }
      }),
      .id = NULL # We are adding Database column manually
    )
    # Filter out NULLs if any DBs were skipped
    if(inherits(combined_enrichr_df, "list")) { # Should not happen if bind_rows works as expected
      combined_enrichr_df <- dplyr::bind_rows(combined_enrichr_df[!sapply(combined_enrichr_df, is.null)])
    }
    
    
    if (!is.null(combined_enrichr_df) && nrow(combined_enrichr_df) > 0) {
      # Ensure enrichr_csv_dir_celltype exists (it should be created by the main loop)
      create_dir_if_needed(enrichr_csv_dir_celltype)
      outfile_path <- file.path(enrichr_csv_dir_celltype, paste0(celltype, "_enrichr_results_", label, ".csv"))
      tryCatch({
        write.csv(combined_enrichr_df, outfile_path, row.names = FALSE)
        message(paste("Saved combined Enrichr results to:", outfile_path))
      }, error = function(e) {
        warning(paste("Failed to save combined Enrichr results:", e$message))
      })
    }
    
    # Plotting
    plot_count_for_pdf <- 0
    pdf_summary_path <- file.path(direction_specific_dir, paste0(celltype, "_", label, "_enrichr_plots_summary.pdf"))
    pdf_dev_is_open <- FALSE
    tryCatch({
      pdf(pdf_summary_path, width = pdf_width, height = pdf_height)
      pdf_dev_is_open <- TRUE
      
      for (db_name_plot in names(enrichr_results_list)) {
        res_df_plot <- enrichr_results_list[[db_name_plot]]
        # Ensure columns for plotting are present
        if(all(c("Term", "Adjusted.P.value", "Overlap") %in% colnames(res_df_plot))) {
          res_df_plot_ready <- res_df_plot %>%
            dplyr::mutate(Count = suppressWarnings(as.numeric(sub("/.*", "", Overlap)))) %>%
            dplyr::filter(!is.na(Count) & Adjusted.P.value < pval_cutoff) %>% # Filter for significance for plotting
            dplyr::arrange(Adjusted.P.value)
          
          if (nrow(res_df_plot_ready) > 0) {
            num_terms_to_display <- min(top_n, nrow(res_df_plot_ready))
            plot_title_text <- paste(celltype, label, "-", db_name_plot, "(Top", num_terms_to_display, "sig.)")
            
            # Use enrichR::plotEnrich
            p_enrich <- tryCatch({
              enrichR::plotEnrich(res_df_plot_ready, showTerms = num_terms_to_display, numChar = 60,
                                  y = "Count", orderBy = "Adjusted.P.value", title = plot_title_text) +
                ggplot2::theme(plot.title = ggplot2::element_text(size = 10)) # Adjust title size
            }, error = function(e){
              warning(paste0("plotEnrich failed for DB: '", db_name_plot, "' in ", celltype, " - ", label, ": ", e$message))
              NULL
            })
            
            if(!is.null(p_enrich)) {
              try(print(p_enrich), silent = TRUE) # Print to PDF
              plot_count_for_pdf <- plot_count_for_pdf + 1
              # Save individual plots
              plot_file_basename <- paste0(celltype, "_", label, "_enrichr_", gsub("[^[:alnum:]_]", "_", db_name_plot), "_barplot")
              for (fmt_ext in plot_formats) {
                plot_full_filename <- paste0(plot_file_basename, ".", fmt_ext)
                plot_full_filepath <- file.path(plots_individual_dir, plot_full_filename)
                tryCatch({
                  ggplot2::ggsave(filename = plot_full_filepath, plot = p_enrich,
                                  width = plot_width, height = plot_height, units = plot_units, dpi = plot_dpi)
                }, error = function(e){
                  warning(paste0("ggsave failed for Enrichr plot (format '", fmt_ext, "'): '",
                                 plot_full_filename, "' => ", conditionMessage(e)))
                })
              }
            }
          }
        } else {
          warning(paste0("Skipping plotting for Enrichr DB '", db_name_plot, "' due to missing essential columns for plotting."))
        }
      }
    }, error = function(e_pdf) {
      warning(paste0("Error during Enrichr PDF generation for ", celltype, " - ", label, ": ", conditionMessage(e_pdf)))
    }, finally = {
      if(pdf_dev_is_open) dev.off()
    })
    
    if (plot_count_for_pdf > 0) {
      message(paste0("Saved ", plot_count_for_pdf, " types of Enrichr plots to multi-page PDF: ", pdf_summary_path))
      message(paste0("Individual Enrichr plots (if any) saved in: ", plots_individual_dir))
    } else if (file.exists(pdf_summary_path)) {
      try(file.remove(pdf_summary_path), silent = TRUE) # Remove empty PDF
      message(paste0("No significant Enrichr terms to plot for ", celltype, " - ", label))
    }
  }
  
  
  ## ---- Helper Function: clusterProfiler + Pathview Integration ----
  run_clusterprofiler_and_pathview <- function(gene_list_symbols, label, celltype, de_table_full_celltype, cp_csv_dir_celltype) {
    pathway_name_cache_cp <- list() # Initialize cache locally for this run
    
    # Define output directories
    celltype_specific_outdir_cp <- file.path(outdir, celltype)
    cp_main_dir_ct <- file.path(celltype_specific_outdir_cp, "clusterProfiler")
    direction_specific_dir_cp <- file.path(cp_main_dir_ct, label)
    plots_individual_dir_cp <- file.path(direction_specific_dir_cp, "plots_individual")
    
    create_dir_if_needed(direction_specific_dir_cp)
    create_dir_if_needed(plots_individual_dir_cp)
    
    # --- Gene ID Mapping & FC Vector Prep ---
    fc_vector_for_cnet <- numeric(0) # For cnetplot, uses only genes in current list
    if (all(c("gene", "avg_log2FC") %in% colnames(de_table_full_celltype))) {
      fc_df_cnet <- de_table_full_celltype %>%
        dplyr::filter(gene %in% gene_list_symbols) %>%
        dplyr::distinct(gene, .keep_all = TRUE) %>% # Keep first if duplicates by gene name
        dplyr::select(gene, avg_log2FC) %>%
        dplyr::filter(!is.na(avg_log2FC) & is.numeric(avg_log2FC))
      
      if(nrow(fc_df_cnet) > 0) {
        fc_vector_for_cnet <- fc_df_cnet$avg_log2FC
        names(fc_vector_for_cnet) <- fc_df_cnet$gene
      }
    } else {
      warning("Missing 'gene' or 'avg_log2FC' in de_table_full_celltype, FC vector for cnetplot will be empty.")
    }
    
    
    gene_entrez_map <- tryCatch({
      clusterProfiler::bitr(gene_list_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) %>%
        dplyr::distinct(SYMBOL, .keep_all = TRUE) # Ensure one ENTREZID per SYMBOL if multiple mappings
    }, error = function(e) {
      warning(paste0("bitr mapping from SYMBOL to ENTREZID failed for ", celltype, " - ", label, ": ", e$message))
      NULL
    })
    
    if (is.null(gene_entrez_map) || nrow(gene_entrez_map) == 0) {
      warning(paste0("No ENTREZ IDs successfully mapped for ", celltype, " - ", label, ". Skipping clusterProfiler."))
      return(NULL)
    }
    
    mapped_gene_count <- nrow(gene_entrez_map)
    total_input_genes <- length(unique(gene_list_symbols))
    if (total_input_genes > 0) {
      mapping_eff <- (mapped_gene_count / total_input_genes) * 100
      if (mapping_eff < 100) {
        message(sprintf("Mapped %d out of %d input genes (%.1f%%) to ENTREZ IDs for %s - %s.",
                        mapped_gene_count, total_input_genes, mapping_eff, celltype, label))
      }
    } else {
      warning(paste0("Input gene list for clusterProfiler is empty for ", celltype, " - ", label))
      return(NULL)
    }
    
    entrez_ids_for_enrichment <- unique(gene_entrez_map$ENTREZID)
    entrez_ids_for_enrichment <- entrez_ids_for_enrichment[!is.na(entrez_ids_for_enrichment)] # Remove any NAs
    
    if(length(entrez_ids_for_enrichment) == 0) {
      warning(paste0("No valid (non-NA) ENTREZ IDs obtained after mapping for ", celltype, " - ", label, ". Skipping clusterProfiler."))
      return(NULL)
    }
    
    # Prepare fc_for_pathview (can be full DE table or subset)
    fc_vector_for_pathview <- NULL
    if (use_full_fc_for_pathview) {
      if (all(c("gene", "p_val_adj", "avg_log2FC") %in% colnames(de_table_full_celltype))) {
        all_sig_genes_df <- de_table_full_celltype %>%
          dplyr::filter(!is.na(p_val_adj) & p_val_adj < pval_cutoff,
                        !is.na(avg_log2FC) & abs(avg_log2FC) >= min(abs(up_log2fc_threshold), abs(down_log2fc_threshold), na.rm = TRUE)) %>%
          dplyr::distinct(gene, .keep_all = TRUE) %>%
          dplyr::select(gene, avg_log2FC)
        
        if(nrow(all_sig_genes_df) > 0) {
          all_sig_fc_entrez_map <- tryCatch({
            clusterProfiler::bitr(all_sig_genes_df$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) %>%
              dplyr::distinct(SYMBOL, .keep_all = TRUE)
          }, error = function(e){
            warning(paste0("bitr mapping failed for full FC vector for Pathview (", celltype, " - ", label, "): ", e$message))
            NULL
          })
          
          if(!is.null(all_sig_fc_entrez_map) && nrow(all_sig_fc_entrez_map) > 0) {
            all_sig_fc_df_entrez <- dplyr::inner_join(all_sig_genes_df, all_sig_fc_entrez_map, by = c("gene" = "SYMBOL")) %>%
              dplyr::filter(!is.na(ENTREZID) & !is.na(avg_log2FC))
            
            if(nrow(all_sig_fc_df_entrez) > 0) {
              # Handle cases where one ENTREZID might map to multiple gene symbols with different FCs by averaging
              all_sig_fc_agg_df <- all_sig_fc_df_entrez %>%
                dplyr::group_by(ENTREZID) %>%
                dplyr::summarise(avg_log2FC = mean(avg_log2FC, na.rm = TRUE), .groups = 'drop')
              fc_vector_for_pathview <- all_sig_fc_agg_df$avg_log2FC
              names(fc_vector_for_pathview) <- all_sig_fc_agg_df$ENTREZID
              message(sprintf("Pathview will use FC data from %d significant genes (full DE table) for %s - %s.",
                              length(fc_vector_for_pathview), celltype, label))
            }
          }
        }
      } else {
        warning("Missing required columns in de_table_full_celltype for 'use_full_fc_for_pathview=TRUE'.")
      }
    } else { # Use subset FC vector (from current gene_list_symbols)
      if (length(fc_vector_for_cnet) > 0) { # fc_vector_for_cnet has SYMBOL names
        subset_fc_entrez_map <- tryCatch({
          clusterProfiler::bitr(names(fc_vector_for_cnet), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) %>%
            dplyr::distinct(SYMBOL, .keep_all = TRUE)
        }, error = function(e){
          warning(paste0("bitr mapping failed for subset FC vector for Pathview (", celltype, " - ", label, "): ", e$message))
          NULL
        })
        
        if(!is.null(subset_fc_entrez_map) && nrow(subset_fc_entrez_map) > 0) {
          fc_vector_cnet_df <- data.frame(SYMBOL = names(fc_vector_for_cnet), avg_log2FC = fc_vector_for_cnet)
          subset_fc_df_entrez <- dplyr::inner_join(fc_vector_cnet_df, subset_fc_entrez_map, by = "SYMBOL") %>%
            dplyr::filter(!is.na(ENTREZID) & !is.na(avg_log2FC))
          
          if(nrow(subset_fc_df_entrez) > 0){
            subset_fc_agg_df <- subset_fc_df_entrez %>%
              dplyr::group_by(ENTREZID) %>%
              dplyr::summarise(avg_log2FC = mean(avg_log2FC, na.rm = TRUE), .groups = 'drop')
            fc_vector_for_pathview <- subset_fc_agg_df$avg_log2FC
            names(fc_vector_for_pathview) <- subset_fc_agg_df$ENTREZID
            message(sprintf("Pathview will use FC data from %d genes (current list subset) for %s - %s.",
                            length(fc_vector_for_pathview), celltype, label))
          }
        }
      }
    }
    
    if(is.null(fc_vector_for_pathview) || length(fc_vector_for_pathview) == 0 || !is.numeric(fc_vector_for_pathview) || all(is.na(fc_vector_for_pathview))) {
      warning(paste0("FC vector for Pathview is invalid or empty for ", celltype, " - ", label, ". Pathview might be skipped or lack coloring."))
      fc_vector_for_pathview <- NULL # Ensure it's NULL if problematic
    }
    
    
    # --- Enrichment analyses ---
    # For GO, use the main pval_cutoff for filtering results for saving/plotting
    ego_bp_res <- run_enrichment(clusterProfiler::enrichGO, gene = entrez_ids_for_enrichment, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH", pvalueCutoff = pval_cutoff, qvalueCutoff = pval_cutoff)
    ego_cc_res <- run_enrichment(clusterProfiler::enrichGO, gene = entrez_ids_for_enrichment, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "CC", pAdjustMethod = "BH", pvalueCutoff = pval_cutoff, qvalueCutoff = pval_cutoff)
    ego_mf_res <- run_enrichment(clusterProfiler::enrichGO, gene = entrez_ids_for_enrichment, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "MF", pAdjustMethod = "BH", pvalueCutoff = pval_cutoff, qvalueCutoff = pval_cutoff)
    
    # For KEGG and Reactome, use a more lenient cutoff if save_all_results is TRUE for table saving
    # Plots will still be filtered by the main pval_cutoff later.
    table_save_pval_cutoff <- if(save_all_results) 1.0 else pval_cutoff
    
    ekegg_res <- run_enrichment(clusterProfiler::enrichKEGG, gene = entrez_ids_for_enrichment, organism = "hsa", pvalueCutoff = table_save_pval_cutoff, qvalueCutoff = table_save_pval_cutoff)
    react_enrich_res <- run_enrichment(ReactomePA::enrichPathway, gene = entrez_ids_for_enrichment, organism = "human", pvalueCutoff = table_save_pval_cutoff, qvalueCutoff = table_save_pval_cutoff, readable = TRUE)
    
    cp_results_all <- list(ego_bp = ego_bp_res, ego_cc = ego_cc_res, ego_mf = ego_mf_res, ekegg = ekegg_res, react_enrich = react_enrich_res)
    cp_results_all <- cp_results_all[!sapply(cp_results_all, is.null)] # Remove NULLs if any enrichment failed
    
    if(length(cp_results_all) == 0) {
      message(paste0("No clusterProfiler results were generated for ", celltype, " - ", label))
      return(NULL)
    }
    
    # --- Save RDS and CSV results ---
    rds_file_path <- file.path(direction_specific_dir_cp, paste0(celltype, "_clusterProfiler_results_", label, ".rds"))
    tryCatch({
      saveRDS(cp_results_all, file = rds_file_path)
      message(paste("Saved combined clusterProfiler results (RDS) to:", rds_file_path))
    }, error = function(e) {
      warning(paste("Failed to save clusterProfiler RDS object:", e$message))
    })
    
    num_csv_saved <- 0
    create_dir_if_needed(cp_csv_dir_celltype) # Ensure celltype-specific CSV dir exists
    
    for(res_type_name in names(cp_results_all)) {
      enrich_obj_to_save <- cp_results_all[[res_type_name]]
      if (!is.null(enrich_obj_to_save) && inherits(enrich_obj_to_save, "enrichResult") && nrow(enrich_obj_to_save@result) > 0) {
        df_to_save <- NULL
        # Make results readable (gene symbols) for GO and KEGG before saving CSV
        if (res_type_name %in% c("ekegg", "ego_bp", "ego_cc", "ego_mf")) {
          df_to_save <- tryCatch({
            clusterProfiler::setReadable(enrich_obj_to_save, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")@result
          }, error = function(e){
            warning(paste0("setReadable failed for ", res_type_name, " during CSV saving for ", celltype, " - ", label, ". Using original IDs. Error: ", e$message))
            as.data.frame(enrich_obj_to_save@result) # Fallback
          })
        } else if (res_type_name == "react_enrich") { # Reactome is already readable if specified
          df_to_save <- as.data.frame(enrich_obj_to_save@result)
        }
        
        if(!is.null(df_to_save) && nrow(df_to_save) > 0) {
          # Filter by pval_cutoff for CSV if save_all_results is FALSE
          if (!save_all_results && "p.adjust" %in% colnames(df_to_save)) {
            df_to_save <- df_to_save %>% dplyr::filter(p.adjust <= pval_cutoff)
          }
          
          if (nrow(df_to_save) > 0) { # Check again after potential filtering
            csv_file_prefix <- switch(res_type_name,
                                      ego_bp = "GO_BP", ego_cc = "GO_CC", ego_mf = "GO_MF",
                                      ekegg = "KEGG", react_enrich = "Reactome", "Unknown_Enrichment")
            csv_file_name_full <- paste0(celltype, "_", csv_file_prefix, "_", label, ".csv")
            csv_file_path_full <- file.path(cp_csv_dir_celltype, csv_file_name_full) # Save to celltype specific dir
            
            tryCatch({
              write.csv(df_to_save, file = csv_file_path_full, row.names = FALSE)
              num_csv_saved <- num_csv_saved + 1
            }, error = function(e){
              warning(paste0("Failed to write CSV for ", res_type_name, " (", celltype, " - ", label, "): ", e$message))
            })
          }
        }
      }
    }
    if(num_csv_saved > 0) {
      message(paste0("Saved ", num_csv_saved, " clusterProfiler result tables to: ", cp_csv_dir_celltype))
    } else {
      message(paste0("No clusterProfiler result tables met criteria for saving for ", celltype, " - ", label))
    }
    
    # --- Create and Save Plots ---
    plots_for_pdf_list <- list()
    
    # Barplots
    for(res_type_name_plot in names(cp_results_all)){
      enrich_obj_plot <- cp_results_all[[res_type_name_plot]]
      if (!is.null(enrich_obj_plot) && inherits(enrich_obj_plot, "enrichResult") && nrow(enrich_obj_plot@result) > 0) {
        # Filter data for plotting based on pval_cutoff
        plot_data_filtered <- enrich_obj_plot
        if ("p.adjust" %in% colnames(plot_data_filtered@result)) {
          plot_data_filtered@result <- plot_data_filtered@result %>% dplyr::filter(p.adjust <= pval_cutoff)
        }
        
        if(nrow(plot_data_filtered@result) > 0 && "Count" %in% colnames(plot_data_filtered@result) && is.numeric(plot_data_filtered@result$Count)) {
          plot_data_filtered@result <- plot_data_filtered@result %>% dplyr::filter(!is.na(Count)) # Ensure Count is not NA
          if(nrow(plot_data_filtered@result) > 0) {
            num_categories_to_show <- min(top_n, nrow(plot_data_filtered@result))
            if (num_categories_to_show > 0) {
              plot_title_str <- paste(celltype, label, res_type_name_plot, "(Top", num_categories_to_show, "sig.)")
              p_bar_plot <- tryCatch({
                enrichplot::barplot(plot_data_filtered, showCategory = num_categories_to_show, x = "Count") +
                  ggplot2::ggtitle(plot_title_str) +
                  ggplot2::theme(plot.title = ggplot2::element_text(size=10))
              }, error = function(e){
                warning(paste0("Barplot generation failed for '", plot_title_str, "': ", conditionMessage(e)))
                NULL
              })
              if(!is.null(p_bar_plot)) plots_for_pdf_list[[paste0(res_type_name_plot, "_bar")]] <- p_bar_plot
            }
          }
        } else {
          # message(paste0("No valid data (after filtering by pval_cutoff or missing Count) for barplot: ", res_type_name_plot, " in ", celltype, " - ", label))
        }
      }
    }
    
    # Cnetplots (typically for GO-BP and Reactome)
    # fc_vector_for_cnet has SYMBOL names, cnetplot needs ENTREZ if enrich_obj is not readable
    # However, setReadable is called inside generate_cnetplot, so SYMBOL names for fc_vector_for_cnet is fine.
    for(res_type_name_cnet in c("ego_bp", "react_enrich")){ # Add others if desired
      enrich_obj_cnet <- cp_results_all[[res_type_name_cnet]]
      if (!is.null(enrich_obj_cnet) && inherits(enrich_obj_cnet, "enrichResult") && nrow(enrich_obj_cnet@result) > 0 && length(fc_vector_for_cnet) > 0) {
        # Filter data for plotting based on pval_cutoff
        plot_data_cnet_filtered <- enrich_obj_cnet
        if ("p.adjust" %in% colnames(plot_data_cnet_filtered@result)) {
          plot_data_cnet_filtered@result <- plot_data_cnet_filtered@result %>% dplyr::filter(p.adjust <= pval_cutoff)
        }
        
        if(nrow(plot_data_cnet_filtered@result) > 0) {
          enrich_readable_cnet <- tryCatch({
            clusterProfiler::setReadable(plot_data_cnet_filtered, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
          }, error = function(e){
            warning(paste0("setReadable failed for cnetplot (", res_type_name_cnet, ") in ", celltype, " - ", label, ". Error: ", e$message))
            plot_data_cnet_filtered # Fallback to non-readable
          })
          
          num_categories_cnet <- min(top_n, nrow(enrich_readable_cnet@result))
          if(num_categories_cnet > 0) {
            # Normal CNetplot
            cnet_title_normal <- paste(celltype, label, "-", res_type_name_cnet, "CNet (Normal, Top", num_categories_cnet, "sig.)")
            p_cnet_normal_plot <- tryCatch({
              enrichplot::cnetplot(enrich_readable_cnet, showCategory = num_categories_cnet,
                                   foldChange = fc_vector_for_cnet, # fc_vector_for_cnet has SYMBOL names
                                   circular = FALSE, colorEdge = cnet_color_edge,
                                   node_label = "all", # "gene", "category", "all", "none"
                                   cex_label_gene = 0.7, cex_label_category = 0.9) +
                ggplot2::ggtitle(cnet_title_normal) + ggplot2::theme(plot.title = ggplot2::element_text(size=10))
            }, error = function(e) {warning(paste("Normal cnetplot failed for",cnet_title_normal,":",conditionMessage(e))); NULL})
            if(!is.null(p_cnet_normal_plot)) plots_for_pdf_list[[paste0(res_type_name_cnet, "_cnet_normal")]] <- p_cnet_normal_plot
            
            # Circular CNetplot
            cnet_title_circular <- paste(celltype, label, "-", res_type_name_cnet, "CNet (Circular, Top", num_categories_cnet, "sig.)")
            p_cnet_circular_plot <- tryCatch({
              enrichplot::cnetplot(enrich_readable_cnet, showCategory = num_categories_cnet,
                                   foldChange = fc_vector_for_cnet,
                                   circular = TRUE, colorEdge = cnet_color_edge,
                                   node_label = "all",
                                   cex_label_gene = 0.7, cex_label_category = 0.9) +
                ggplot2::ggtitle(cnet_title_circular) + ggplot2::theme(plot.title = ggplot2::element_text(size=10))
            }, error = function(e) {warning(paste("Circular cnetplot failed for",cnet_title_circular,":",conditionMessage(e))); NULL})
            if(!is.null(p_cnet_circular_plot)) plots_for_pdf_list[[paste0(res_type_name_cnet, "_cnet_circular")]] <- p_cnet_circular_plot
          }
        }
      } else {
        # message(paste0("Skipping cnetplot for ", res_type_name_cnet, " in ", celltype, " - ", label, " due to no significant terms or empty FC vector."))
      }
    }
    plots_for_pdf_list <- plots_for_pdf_list[!sapply(plots_for_pdf_list, is.null)] # Remove NULLs
    
    # Save plots to multi-page PDF and individual files
    if (length(plots_for_pdf_list) > 0) {
      pdf_summary_path_cp <- file.path(direction_specific_dir_cp, paste0(celltype, "_", label, "_clusterProfiler_plots_summary.pdf"))
      pdf_dev_cp_is_open <- FALSE
      tryCatch({
        pdf(pdf_summary_path_cp, width = pdf_width, height = pdf_height)
        pdf_dev_cp_is_open <- TRUE
        for(plot_obj_to_print in plots_for_pdf_list) {
          try(print(plot_obj_to_print), silent = TRUE)
        }
      }, error = function(e_pdf_cp) {
        warning(paste0("Failed to create multi-page PDF for clusterProfiler plots (", celltype, " - ", label, "): ", conditionMessage(e_pdf_cp)))
      }, finally = {
        if(pdf_dev_cp_is_open) dev.off()
      })
      
      if (file.exists(pdf_summary_path_cp) && file.info(pdf_summary_path_cp)$size > 0) {
        message(paste0("Saved multi-page clusterProfiler summary PDF to: ", pdf_summary_path_cp))
      } else if (file.exists(pdf_summary_path_cp)) {
        try(file.remove(pdf_summary_path_cp), silent = TRUE) # Remove empty PDF
      }
      
      
      num_individual_plots_cp <- 0
      for (plot_item_name in names(plots_for_pdf_list)) {
        plot_obj_to_save <- plots_for_pdf_list[[plot_item_name]]
        plot_file_basename_cp <- paste0(celltype, "_", label, "_clusterProfiler_", plot_item_name)
        plot_saved_this_type <- FALSE
        for (fmt_ext_cp in plot_formats) {
          plot_full_filename_cp <- paste0(plot_file_basename_cp, ".", fmt_ext_cp)
          plot_full_filepath_cp <- file.path(plots_individual_dir_cp, plot_full_filename_cp)
          tryCatch({
            ggplot2::ggsave(filename = plot_full_filepath_cp, plot = plot_obj_to_save,
                            width = plot_width, height = plot_height, units = plot_units, dpi = plot_dpi)
            if(!plot_saved_this_type) { num_individual_plots_cp <- num_individual_plots_cp + 1; plot_saved_this_type <- TRUE}
          }, error = function(e_ggsave){
            warning(paste0("ggsave failed for clusterProfiler plot (format '", fmt_ext_cp, "'): '",
                           plot_full_filename_cp, "' => ", conditionMessage(e_ggsave)))
          })
        }
      }
      if (num_individual_plots_cp > 0) {
        message(paste0("Saved ", num_individual_plots_cp, " types of individual clusterProfiler plots to: ", plots_individual_dir_cp))
      }
    } else {
      message(paste0("No valid clusterProfiler plots generated for ", celltype, " - ", label))
    }
    
    # --- KEGG Pathview Integration ---
    if (!is.null(ekegg_res) && inherits(ekegg_res, "enrichResult") && nrow(ekegg_res@result) > 0 && !is.null(fc_vector_for_pathview)) {
      pathview_main_dir_ct <- file.path(celltype_specific_outdir_cp, "pathview") # pathview folder at celltype level
      pathview_subdir_directional <- file.path(pathview_main_dir_ct, label) # Direction specific
      create_dir_if_needed(pathview_subdir_directional)
      
      ekegg_results_df_pv <- ekegg_res@result
      # Filter for significance for Pathview visualization
      if ("p.adjust" %in% colnames(ekegg_results_df_pv)) {
        ekegg_results_df_pv_filtered <- ekegg_results_df_pv %>% dplyr::filter(p.adjust <= pval_cutoff)
      } else {
        warning("p.adjust column missing in KEGG results, cannot filter for Pathview based on significance.")
        ekegg_results_df_pv_filtered <- ekegg_results_df_pv # Proceed with unfiltered if p.adjust is missing
      }
      
      
      if(nrow(ekegg_results_df_pv_filtered) > 0) {
        if(!all(c("ID", "Description") %in% colnames(ekegg_results_df_pv_filtered))) {
          warning(paste0("Missing 'ID' or 'Description' columns in KEGG results for Pathview (", celltype, " - ", label, ")."))
        } else {
          ekegg_ordered_for_pv <- ekegg_results_df_pv_filtered %>% dplyr::arrange(p.adjust) # Sort by p.adjust
          top_kegg_pathways_to_plot <- head(ekegg_ordered_for_pv, min(kegg_n, nrow(ekegg_ordered_for_pv)))
          
          if (nrow(top_kegg_pathways_to_plot) > 0) {
            message(paste0("Generating Pathview visualizations for top ", nrow(top_kegg_pathways_to_plot),
                           " significant KEGG pathways for ", celltype, " - ", label, "..."))
            num_pathview_generated <- 0
            pathview_run_statuses <- logical(nrow(top_kegg_pathways_to_plot))
            
            for(idx in 1:nrow(top_kegg_pathways_to_plot)) {
              current_kegg_id <- top_kegg_pathways_to_plot$ID[idx]
              current_kegg_desc <- top_kegg_pathways_to_plot$Description[idx]
              message(paste0("  Processing KEGG pathway for Pathview: ", current_kegg_id, " (", current_kegg_desc, ")"))
              
              pathview_out_suffix <- paste0(celltype, "_", label) # Suffix for pathview output files
              
              pv_run_result <- run_and_rename_pathview(
                kegg_id = current_kegg_id,
                kegg_desc = current_kegg_desc,
                fc_vector = fc_vector_for_pathview, # This is ENTREZ named
                out_suffix_pv = pathview_out_suffix,
                target_dir = pathview_subdir_directional,
                local_pathway_name_cache = pathway_name_cache_cp # Pass and update cache
              )
              pathway_name_cache_cp <- pv_run_result$cache # Update cache from helper
              pathview_run_statuses[idx] <- pv_run_result$generated
            }
            num_pathview_generated <- sum(pathview_run_statuses)
            if(num_pathview_generated > 0) {
              message(paste0("Successfully generated ", num_pathview_generated, " Pathview images in: ", pathview_subdir_directional))
            } else {
              message(paste0("No Pathview images were generated for ", celltype, " - ", label, " (check warnings)."))
            }
          } else {
            message(paste0("No significant KEGG pathways met criteria for Pathview visualization for ", celltype, " - ", label))
          }
        }
      } else {
        message(paste0("No KEGG results passed p.adjust filtering (<= ", pval_cutoff, ") for Pathview for ", celltype, " - ", label))
      }
    } else {
      if(is.null(fc_vector_for_pathview)) {
        message(paste0("Skipping Pathview for ", celltype, " - ", label, " due to invalid or empty fold change vector."))
      } else {
        message(paste0("Skipping Pathview for ", celltype, " - ", label, " as no significant KEGG enrichment results were found."))
      }
    }
    return(cp_results_all) # Return the full list of CP results for this direction
  } # end run_clusterprofiler_and_pathview
  
  
  # --- Setup Parallel Backend ---
  if (parallel.per.celltype) {
    num_workers <- n_cores
    if (is.null(num_workers) || !is.numeric(num_workers) || num_workers <= 0) { # Ensure num_workers is valid
      num_workers <- max(1, future::availableCores(omit = 1)) # Leave one core free
    }
    message(paste0("Setting up parallel processing for cell types using ", num_workers, " workers (via future::multisession)."))
    future::plan(future::multisession, workers = num_workers)
    on.exit(future::plan(future::sequential), add = TRUE) # Important to reset plan
  } else {
    message("Running analysis sequentially for each cell type.")
    future::plan(future::sequential) # Explicitly set sequential plan
  }
  
  # --- Main Processing Loop (Parallel or Sequential) ---
  all_celltype_results <- future.apply::future_lapply(names(de_results_list), function(current_celltype_name) {
    # Ensure necessary packages are loaded within the worker/session
    # This is good practice for future_lapply, though loaded above, environments can differ.
    suppressPackageStartupMessages({
      library(dplyr); library(clusterProfiler); library(enrichR); library(ReactomePA);
      library(enrichplot); library(ggplot2); library(pathview); library(KEGGREST); library(org.Hs.eg.db)
    })
    
    message(paste0("\n================================================="))
    message(paste0("Processing cell type: ", current_celltype_name, " (Worker PID: ", Sys.getpid(), ")"))
    message(paste0("================================================="))
    
    current_de_table <- de_results_list[[current_celltype_name]]
    
    # Basic validation of the current DE table
    if (!all(c("gene", "p_val_adj", "avg_log2FC") %in% colnames(current_de_table))) {
      warning(paste0("Skipping cell type '", current_celltype_name, "' due to missing required columns in its DE table."))
      return(NULL) # Return NULL for this cell type if data is bad
    }
    # Ensure numeric types and filter out NA genes or empty gene names
    current_de_table <- current_de_table %>%
      dplyr::mutate(
        p_val_adj = suppressWarnings(as.numeric(as.character(p_val_adj))), # Handle factors if any
        avg_log2FC = suppressWarnings(as.numeric(as.character(avg_log2FC)))
      ) %>%
      dplyr::filter(!is.na(gene) & gene != "" & !is.na(p_val_adj) & !is.na(avg_log2FC))
    
    if(nrow(current_de_table) == 0) {
      warning(paste0("Skipping cell type '", current_celltype_name, "' as it has no valid data rows after basic cleaning."))
      return(NULL)
    }
    
    # Create celltype-specific output directories
    celltype_outdir_main <- file.path(outdir, current_celltype_name)
    create_dir_if_needed(celltype_outdir_main)
    
    # Define and create celltype-specific CSV directories (within the main csv_dir_base)
    celltype_csv_main_dir <- file.path(csv_dir_base, current_celltype_name)
    enrichr_csv_for_celltype <- file.path(celltype_csv_main_dir, "Enrichr_CSVs")
    cp_csv_for_celltype <- file.path(celltype_csv_main_dir, "ClusterProfiler_CSVs")
    create_dir_if_needed(celltype_csv_main_dir)
    create_dir_if_needed(enrichr_csv_for_celltype)
    create_dir_if_needed(cp_csv_for_celltype)
    
    # Define filtering thresholds (use defaults if parameters are not numeric, though validation should catch this)
    up_fc_thresh_eff <- if(is.numeric(up_log2fc_threshold)) up_log2fc_threshold else 0.25
    down_fc_thresh_eff <- if(is.numeric(down_log2fc_threshold)) down_log2fc_threshold else -0.25
    pval_cutoff_eff <- if(is.numeric(pval_cutoff)) pval_cutoff else 0.05
    
    # --- Gene Selection with max_genes_for_enrichment ---
    # Get all significant up-regulated genes, sorted
    all_sig_up_df <- current_de_table %>%
      dplyr::filter(p_val_adj < pval_cutoff_eff, avg_log2FC > up_fc_thresh_eff) %>%
      dplyr::arrange(p_val_adj, desc(abs(avg_log2FC))) # Sort by p-value, then by magnitude of FC
    
    # Get all significant down-regulated genes, sorted
    all_sig_down_df <- current_de_table %>%
      dplyr::filter(p_val_adj < pval_cutoff_eff, avg_log2FC < down_fc_thresh_eff) %>%
      dplyr::arrange(p_val_adj, desc(abs(avg_log2FC))) # Sort by p-value, then by magnitude of FC
    
    selected_sig_up_genes <- character(0) # Initialize
    if (!is.null(max_genes_for_enrichment) && is.numeric(max_genes_for_enrichment) && max_genes_for_enrichment > 0 && !is.infinite(max_genes_for_enrichment)) {
      if (nrow(all_sig_up_df) > 0) {
        num_to_select_up <- min(max_genes_for_enrichment, nrow(all_sig_up_df))
        selected_sig_up_genes <- head(all_sig_up_df$gene, num_to_select_up)
        message(sprintf("For %s (Up): Selected top %d (out of %d initially significant) genes based on p_val_adj and |log2FC|.",
                        current_celltype_name, length(selected_sig_up_genes), nrow(all_sig_up_df)))
      } else {
        message(sprintf("For %s (Up): No genes passed initial significance filters (p < %.2f, log2FC > %.2f).",
                        current_celltype_name, pval_cutoff_eff, up_fc_thresh_eff))
      }
    } else { # Use all significant genes
      if (nrow(all_sig_up_df) > 0) {
        selected_sig_up_genes <- all_sig_up_df$gene
        message(sprintf("For %s (Up): Using all %d significant genes (p < %.2f, log2FC > %.2f).",
                        current_celltype_name, length(selected_sig_up_genes), pval_cutoff_eff, up_fc_thresh_eff))
      } else {
        message(sprintf("For %s (Up): No genes passed initial significance filters (p < %.2f, log2FC > %.2f).",
                        current_celltype_name, pval_cutoff_eff, up_fc_thresh_eff))
      }
    }
    selected_sig_up_genes <- unique(selected_sig_up_genes) # Ensure uniqueness
    
    selected_sig_down_genes <- character(0) # Initialize
    if (!is.null(max_genes_for_enrichment) && is.numeric(max_genes_for_enrichment) && max_genes_for_enrichment > 0 && !is.infinite(max_genes_for_enrichment)) {
      if (nrow(all_sig_down_df) > 0) {
        num_to_select_down <- min(max_genes_for_enrichment, nrow(all_sig_down_df))
        selected_sig_down_genes <- head(all_sig_down_df$gene, num_to_select_down)
        message(sprintf("For %s (Down): Selected top %d (out of %d initially significant) genes based on p_val_adj and |log2FC|.",
                        current_celltype_name, length(selected_sig_down_genes), nrow(all_sig_down_df)))
      } else {
        message(sprintf("For %s (Down): No genes passed initial significance filters (p < %.2f, log2FC < %.2f).",
                        current_celltype_name, pval_cutoff_eff, down_fc_thresh_eff))
      }
    } else { # Use all significant genes
      if (nrow(all_sig_down_df) > 0) {
        selected_sig_down_genes <- all_sig_down_df$gene
        message(sprintf("For %s (Down): Using all %d significant genes (p < %.2f, log2FC < %.2f).",
                        current_celltype_name, length(selected_sig_down_genes), pval_cutoff_eff, down_fc_thresh_eff))
      } else {
        message(sprintf("For %s (Down): No genes passed initial significance filters (p < %.2f, log2FC < %.2f).",
                        current_celltype_name, pval_cutoff_eff, down_fc_thresh_eff))
      }
    }
    selected_sig_down_genes <- unique(selected_sig_down_genes) # Ensure uniqueness
    
    if (length(selected_sig_up_genes) == 0 && length(selected_sig_down_genes) == 0) {
      warning(paste0("No significant genes selected for enrichment for cell type '", current_celltype_name,
                     "' after applying all filters (p-value, log2FC, max_genes_for_enrichment). Skipping further analysis for this cell type."))
      return(NULL)
    }
    
    # Save the *selected* gene lists to celltype_csv_main_dir (not Enrichr/CP specific subdirs)
    if (length(selected_sig_up_genes) > 0) {
      tryCatch(write.csv(data.frame(gene = selected_sig_up_genes),
                         file.path(celltype_csv_main_dir, paste0(current_celltype_name, "_UpRegulated_Selected_For_Enrichment.csv")), row.names = FALSE),
               error = function(e) warning(paste0("Failed to save selected Up-regulated gene list for ", current_celltype_name, ": ", e$message)))
    }
    if (length(selected_sig_down_genes) > 0) {
      tryCatch(write.csv(data.frame(gene = selected_sig_down_genes),
                         file.path(celltype_csv_main_dir, paste0(current_celltype_name, "_DownRegulated_Selected_For_Enrichment.csv")), row.names = FALSE),
               error = function(e) warning(paste0("Failed to save selected Down-regulated gene list for ", current_celltype_name, ": ", e$message)))
    }
    
    # Run Enrichr Analysis
    message(paste0("\n--- Running Enrichr Analysis for ", current_celltype_name, " ---"))
    if (length(selected_sig_up_genes) > 0) {
      message(paste0("Processing UpRegulated (", length(selected_sig_up_genes), " genes) for Enrichr..."))
      run_enrichr_analysis(selected_sig_up_genes, "UpRegulated", current_celltype_name, enrichr_csv_for_celltype)
    } else {
      message(paste0("Skipping Enrichr for UpRegulated in ", current_celltype_name, " (no genes selected)."))
    }
    if (length(selected_sig_down_genes) > 0) {
      message(paste0("Processing DownRegulated (", length(selected_sig_down_genes), " genes) for Enrichr..."))
      run_enrichr_analysis(selected_sig_down_genes, "DownRegulated", current_celltype_name, enrichr_csv_for_celltype)
    } else {
      message(paste0("Skipping Enrichr for DownRegulated in ", current_celltype_name, " (no genes selected)."))
    }
    
    # Run clusterProfiler and Pathview Analysis
    message(paste0("\n--- Running clusterProfiler/Pathview Analysis for ", current_celltype_name, " ---"))
    cp_res_up <- NULL; cp_res_down <- NULL
    if (length(selected_sig_up_genes) > 0) {
      message(paste0("Processing UpRegulated (", length(selected_sig_up_genes), " genes) for clusterProfiler/Pathview..."))
      cp_res_up <- run_clusterprofiler_and_pathview(selected_sig_up_genes, "UpRegulated", current_celltype_name, current_de_table, cp_csv_for_celltype)
    } else {
      message(paste0("Skipping clusterProfiler/Pathview for UpRegulated in ", current_celltype_name, " (no genes selected)."))
    }
    if (length(selected_sig_down_genes) > 0) {
      message(paste0("Processing DownRegulated (", length(selected_sig_down_genes), " genes) for clusterProfiler/Pathview..."))
      cp_res_down <- run_clusterprofiler_and_pathview(selected_sig_down_genes, "DownRegulated", current_celltype_name, current_de_table, cp_csv_for_celltype)
    } else {
      message(paste0("Skipping clusterProfiler/Pathview for DownRegulated in ", current_celltype_name, " (no genes selected)."))
    }
    
    # Run Bidirectional Pathview Analysis (uses the selected up/down genes)
    pathway_name_cache_bidir_worker <- list() # Local cache for this worker's bidirectional run
    if (length(selected_sig_up_genes) > 0 && length(selected_sig_down_genes) > 0) {
      message(paste0("\n--- Running Bidirectional KEGG Pathview Analysis for ", current_celltype_name, " ---"))
      pathview_main_dir_bidir <- file.path(celltype_outdir_main, "pathview") # Main pathview dir for celltype
      bidir_pathview_output_dir <- file.path(pathview_main_dir_bidir, "Bidirectional")
      create_dir_if_needed(bidir_pathview_output_dir)
      
      combined_genes_for_bidir <- unique(c(selected_sig_up_genes, selected_sig_down_genes))
      message(paste0("Total unique genes for bidirectional Pathview (from selected up/down lists): ", length(combined_genes_for_bidir)))
      
      gene_entrez_map_combined <- tryCatch({
        clusterProfiler::bitr(combined_genes_for_bidir, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) %>%
          dplyr::distinct(SYMBOL, .keep_all = TRUE)
      }, error = function(e){
        warning(paste0("bitr mapping failed for combined genes for bidirectional Pathview (", current_celltype_name, "): ", e$message))
        NULL
      })
      
      if (!is.null(gene_entrez_map_combined) && nrow(gene_entrez_map_combined) > 0) {
        entrez_ids_combined_bidir <- unique(gene_entrez_map_combined$ENTREZID)
        entrez_ids_combined_bidir <- entrez_ids_combined_bidir[!is.na(entrez_ids_combined_bidir)]
        
        if(length(entrez_ids_combined_bidir) > 0) {
          # KEGG enrichment for the combined list
          ekegg_res_combined <- run_enrichment(clusterProfiler::enrichKEGG,
                                               gene = entrez_ids_combined_bidir, organism = "hsa",
                                               pvalueCutoff = if(save_all_results) 1.0 else pval_cutoff_eff, # Lenient for table
                                               qvalueCutoff = if(save_all_results) 1.0 else pval_cutoff_eff)
          
          if (!is.null(ekegg_res_combined) && inherits(ekegg_res_combined, "enrichResult") && nrow(ekegg_res_combined@result) > 0) {
            ekegg_readable_df_combined <- tryCatch({
              clusterProfiler::setReadable(ekegg_res_combined, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")@result
            }, error = function(e){
              warning(paste0("setReadable failed for combined KEGG results (", current_celltype_name, "): ", e$message))
              as.data.frame(ekegg_res_combined@result) # Fallback
            })
            
            # Save bidirectional KEGG CSV to the celltype's ClusterProfiler CSV directory
            create_dir_if_needed(cp_csv_for_celltype) # Ensure it exists
            bidir_csv_path <- file.path(cp_csv_for_celltype, paste0(current_celltype_name, "_KEGG_Bidirectional_SelectedGenes.csv"))
            tryCatch({
              write.csv(ekegg_readable_df_combined, bidir_csv_path, row.names = FALSE)
              message(paste0("Saved Bidirectional KEGG enrichment results (from selected genes) to: ", bidir_csv_path))
            }, error = function(e){
              warning(paste0("Failed to write Bidirectional KEGG CSV for ", current_celltype_name, ": ", e$message))
            })
            
            # Filter for Pathview plotting (use main pval_cutoff)
            ekegg_df_for_pv_bidir <- ekegg_res_combined@result
            if ("p.adjust" %in% colnames(ekegg_df_for_pv_bidir)) {
              ekegg_df_for_pv_bidir_filtered <- ekegg_df_for_pv_bidir %>% dplyr::filter(p.adjust <= pval_cutoff_eff)
            } else {
              ekegg_df_for_pv_bidir_filtered <- ekegg_df_for_pv_bidir # No p.adjust, proceed with caution
            }
            
            
            if(nrow(ekegg_df_for_pv_bidir_filtered) > 0) {
              if(!all(c("ID", "Description") %in% colnames(ekegg_df_for_pv_bidir_filtered))) {
                warning(paste0("Missing 'ID' or 'Description' in Bidirectional KEGG results for Pathview (", current_celltype_name, ")."))
              } else {
                ekegg_ordered_bidir_pv <- ekegg_df_for_pv_bidir_filtered %>% dplyr::arrange(p.adjust)
                top_kegg_bidir_to_plot <- head(ekegg_ordered_bidir_pv, min(kegg_n, nrow(ekegg_ordered_bidir_pv)))
                
                # Prepare FC vector for bidirectional Pathview (using combined_genes_for_bidir and original de_table)
                fc_df_for_bidir_pv <- current_de_table %>%
                  dplyr::filter(gene %in% combined_genes_for_bidir, !is.na(avg_log2FC)) %>% # Filter for the combined gene list
                  dplyr::distinct(gene, .keep_all = TRUE) %>%
                  dplyr::select(gene, avg_log2FC)
                
                fc_vector_bidir_entrez <- NULL
                if(nrow(fc_df_for_bidir_pv) > 0) {
                  fc_map_bidir_pv <- tryCatch({
                    clusterProfiler::bitr(fc_df_for_bidir_pv$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) %>%
                      dplyr::distinct(SYMBOL, .keep_all = TRUE)
                  }, error = function(e){ warning(paste0("bitr mapping failed for bidirectional Pathview FC vector: ", e$message)); NULL })
                  
                  if(!is.null(fc_map_bidir_pv) && nrow(fc_map_bidir_pv) > 0) {
                    fc_data_bidir_entrez <- dplyr::inner_join(fc_df_for_bidir_pv, fc_map_bidir_pv, by = c("gene" = "SYMBOL")) %>%
                      dplyr::filter(!is.na(ENTREZID) & !is.na(avg_log2FC)) %>%
                      dplyr::group_by(ENTREZID) %>%
                      dplyr::summarise(avg_log2FC = mean(avg_log2FC, na.rm = TRUE), .groups = 'drop')
                    if(nrow(fc_data_bidir_entrez) > 0){
                      fc_vector_bidir_entrez <- fc_data_bidir_entrez$avg_log2FC
                      names(fc_vector_bidir_entrez) <- fc_data_bidir_entrez$ENTREZID
                    }
                  }
                }
                
                if (!is.null(fc_vector_bidir_entrez) && length(fc_vector_bidir_entrez) > 0 && is.numeric(fc_vector_bidir_entrez) && !all(is.na(fc_vector_bidir_entrez))) {
                  message(paste0("Generated FC vector for bidirectional Pathview with ", length(fc_vector_bidir_entrez), " ENTREZ IDs for ", current_celltype_name, "."))
                  if (nrow(top_kegg_bidir_to_plot) > 0) {
                    message(paste0("Generating Bidirectional Pathview visualizations for top ", nrow(top_kegg_bidir_to_plot), " KEGG pathways for ", current_celltype_name, "..."))
                    num_pathview_bidir_generated <- 0
                    bidir_pv_statuses <- logical(nrow(top_kegg_bidir_to_plot))
                    
                    for(j_idx in 1:nrow(top_kegg_bidir_to_plot)) {
                      kegg_id_bidir <- top_kegg_bidir_to_plot$ID[j_idx]
                      kegg_desc_bidir <- top_kegg_bidir_to_plot$Description[j_idx]
                      message(paste0("  Processing Bidirectional KEGG for Pathview: ", kegg_id_bidir, " (", kegg_desc_bidir, ")"))
                      pathview_out_suffix_bidir <- paste0(current_celltype_name, "_Bidir")
                      
                      pv_run_res_bidir <- run_and_rename_pathview(
                        kegg_id = kegg_id_bidir, kegg_desc = kegg_desc_bidir,
                        fc_vector = fc_vector_bidir_entrez,
                        out_suffix_pv = pathview_out_suffix_bidir,
                        target_dir = bidir_pathview_output_dir,
                        local_pathway_name_cache = pathway_name_cache_bidir_worker
                      )
                      pathway_name_cache_bidir_worker <- pv_run_res_bidir$cache
                      bidir_pv_statuses[j_idx] <- pv_run_res_bidir$generated
                    }
                    num_pathview_bidir_generated <- sum(bidir_pv_statuses)
                    if(num_pathview_bidir_generated > 0) {
                      message(paste0("Successfully generated ", num_pathview_bidir_generated, " Bidirectional Pathview images in: ", bidir_pathview_output_dir))
                    } else {
                      message(paste0("No Bidirectional Pathview images were generated for ", current_celltype_name, " (check warnings)."))
                    }
                  } else {
                    message(paste0("No Bidirectional KEGG pathways met criteria for Pathview after filtering for ", current_celltype_name, "."))
                  }
                } else {
                  warning(paste0("Bidirectional Pathview FC vector is invalid or empty for ", current_celltype_name, ". Skipping Pathview for bidirectional results."))
                }
              }
            } else {
              message(paste0("No combined KEGG results passed p.adjust filtering for Pathview for ", current_celltype_name, "."))
            }
          } else {
            message(paste0("No combined KEGG enrichment results found for bidirectional Pathview for ", current_celltype_name, "."))
          }
        } else {
          warning(paste0("Could not map combined genes to ENTREZ IDs for bidirectional Pathview for ", current_celltype_name, "."))
        }
      } else {
        warning(paste0("Failed to map combined genes for bidirectional Pathview for ", current_celltype_name, "."))
      }
    } else {
      if(length(selected_sig_up_genes) > 0 || length(selected_sig_down_genes) > 0) { # Only message if at least one direction had genes
        message(paste0("\nSkipping Bidirectional Pathview for ", current_celltype_name, " (requires both Up- and Down-regulated selected genes)."))
      }
    }
    
    # Return results for this cell type (e.g., the clusterProfiler objects)
    # This structure can be adjusted based on what's most useful to return
    list(
      CellTypeName = current_celltype_name,
      UpRegulated_CP_Results = cp_res_up,
      DownRegulated_CP_Results = cp_res_down
      # Enrichr results are saved to files but not typically returned directly in large lists
    )
    
  }, future.seed = TRUE) # future.seed=TRUE for reproducibility if using random numbers within future (not typical here)
  
  # Assign names to the results list based on input de_results_list names
  names(all_celltype_results) <- names(de_results_list)
  
  message(paste0("\n================================================="))
  message(paste0("Enrichment analysis pipeline complete. All results saved in: ", normalizePath(outdir)))
  message(paste0("================================================="))
  
  invisible(all_celltype_results) # Return combined results invisibly
}

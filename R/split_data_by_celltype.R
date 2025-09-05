#' Split Gene Expression Data by Cell Type
#' 
#' @description
#' This function takes a gene expression data frame where column names are formatted as 
#' "CELLTYPE_DISEASE_SUBJECT" and transforms it into a list of tibbles, one for each cell type.
#' Each tibble contains the gene expression data in long format with columns for gene, cell type,
#' disease, subject, and expression values.
#'
#' @param df A data frame or tibble containing gene expression data. Must have a column named "gene"
#'        and other columns named in the format "CELLTYPE_DISEASE_SUBJECT".
#' @return A named list where each element is a tibble containing data for a specific cell type.
#'        Each tibble has columns: gene, cell_type, disease, subject, average_expression.
#'
#' @examples
#' \dontrun{
#' # Example data frame
#' gene_data <- data.frame(
#'   gene = c("GENE1", "GENE2"),
#'   AT2_CTRL_S1 = c(5.2, 3.1),
#'   AT2_PAH_S2 = c(6.7, 4.2),
#'   AT1_CTRL_S1 = c(2.1, 1.5)
#' )
#' 
#' # Split by cell type
#' cell_type_data <- split_data_by_cell_type2(gene_data)
#' 
#' # Access AT2 data
#' at2_data <- cell_type_data[["AT2"]]
#' }
#'
#' @importFrom dplyr select left_join all_of
#' @importFrom tidyr pivot_longer
#' @export
split_data_by_cell_type <- function(df) {
  # Load necessary packages
  require(dplyr)
  require(tidyr)
  
  # Check if input is a data frame
  if (!is.data.frame(df)) {
    stop("Input must be a data frame or tibble")
  }
  
  # Check if 'gene' column exists
  if (!"gene" %in% colnames(df)) {
    stop("Input data frame must contain a 'gene' column")
  }
  
  # Get all column names except 'gene'
  all_cols <- colnames(df)[colnames(df) != "gene"]
  
  # Create a data frame to map column names to their components
  col_info <- data.frame(
    column = all_cols,
    stringsAsFactors = FALSE
  )
  
  # Extract cell type (everything before the first underscore)
  col_info$cell_type <- sapply(col_info$column, function(x) {
    strsplit(x, "_")[[1]][1]
  })
  
  # Extract disease (always the second part after splitting by underscore)
  col_info$disease <- sapply(col_info$column, function(x) {
    parts <- strsplit(x, "_")[[1]]
    if (length(parts) < 2) {
      warning(paste("Column", x, "does not follow expected format (CELLTYPE_DISEASE_SUBJECT)"))
      return(NA)
    }
    parts[2]
  })
  
  # Extract subject (everything after the disease part)
  col_info$subject <- sapply(col_info$column, function(x) {
    parts <- strsplit(x, "_")[[1]]
    if (length(parts) < 3) {
      warning(paste("Column", x, "does not follow expected format (CELLTYPE_DISEASE_SUBJECT)"))
      return(NA)
    }
    # Join all parts after the 2nd one
    paste(parts[3:length(parts)], collapse = "_")
  })
  
  # Create a list to store results by cell type
  cell_type_data <- list()
  
  # Get unique cell types
  cell_types <- unique(col_info$cell_type)
  
  # For each cell type, create a tibble with the data
  for (ct in cell_types) {
    # Get columns for current cell type
    ct_cols <- col_info$column[col_info$cell_type == ct]
    ct_info <- col_info[col_info$cell_type == ct, ]
    
    # Select gene and relevant columns for this cell type
    ct_df <- df %>% 
      dplyr::select(gene, dplyr::all_of(ct_cols))
    
    # Convert to long format
    ct_long <- ct_df %>% 
      tidyr::pivot_longer(
        cols = -gene, 
        names_to = "column", 
        values_to = "average_expression"
      )
    
    # Join with the column info to get cell_type, disease, and subject
    ct_long <- ct_long %>% 
      dplyr::left_join(ct_info, by = "column") %>% 
      dplyr::select(gene, cell_type, disease, subject, average_expression)
    
    # Add to the list
    cell_type_data[[ct]] <- ct_long
  }
  
  # Add summary attribute
  attr(cell_type_data, "summary") <- list(
    n_cell_types = length(cell_types),
    cell_types = cell_types,
    n_genes = nrow(df),
    diseases = unique(col_info$disease),
    subjects = unique(col_info$subject)
  )
  
  return(cell_type_data)
}
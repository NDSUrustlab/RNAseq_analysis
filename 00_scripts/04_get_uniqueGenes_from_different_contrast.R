library(tidyverse)
library(dplyr)
################################################################################

###################fucntions###################
#function to get unique genes
get_n_save_unique_genes <- function(path, comp1, comp2, sufix, result_dir) {
  # Construct file paths
  file1_path <- paste0(path, "/DE_", comp1, sufix)
  file2_path <- paste0(path, "/DE_", comp2, sufix)

  # Read the CSV files
  file1 <- read.csv(file1_path)
  file2 <- read.csv(file2_path)

  # Check if the files have been read correctly
  if (is.null(file1) || is.null(file2)) {
    stop("One or both of the files could not be read.")
  }

  # Find unique genes in comp1 compared to comp2
  unique_genes_in_comp1 <- anti_join(file1, file2, by = "Gene")

  # Construct the output file path
  output_file_path <- paste0(result_dir, "/unique_DE_in_", comp1, "_compared_to_", comp2, sufix)

  # Write the unique genes to the output file
  write.csv(unique_genes_in_comp1, output_file_path, row.names = FALSE)

  # Return a message indicating successful completion
  return(paste("Unique genes in", comp1, "compared to", comp2, "have been saved to", output_file_path))
}

######################################
# Define the directory path
results_dir_with_log2FC_threshold <- "BLS_DE_Results/01_Res-vs-Sus/With_log2FC_threshold/01_unique_genes"

# Check if the directory exists, and if not, create it
if (!dir.exists(results_dir_with_log2FC_threshold)) {
  dir.create(results_dir_with_log2FC_threshold, recursive = TRUE)
}

# Res vs Sus
# get unique genes in trt compared to mock
get_n_save_unique_genes("BLS_DE_Results/01_Res-vs-Sus/With_log2FC_threshold", "U17BI1_vs_U89BI1", "U17BM1_vs_U89BM1", "_padj-0.05_log2foldChange-1.csv", results_dir_with_log2FC_threshold)
get_n_save_unique_genes("BLS_DE_Results/01_Res-vs-Sus/With_log2FC_threshold", "U17BI3_vs_U89BI3", "U17BM3_vs_U89BM3", "_padj-0.05_log2foldChange-1.csv", results_dir_with_log2FC_threshold)

######################################
# Define the directory path
results_dir_with_log2FC_threshold <- "BLS_DE_Results/02_Trt-vs-Mock/With_log2FC_threshold/01_unique_genes"

# Check if the directory exists, and if not, create it
if (!dir.exists(results_dir_with_log2FC_threshold)) {
  dir.create(results_dir_with_log2FC_threshold, recursive = TRUE)
}

# Res vs Sus
# get unique genes in trt compared to mock
get_n_save_unique_genes("BLS_DE_Results/02_Trt-vs-Mock/With_log2FC_threshold", "U17BI1_vs_U17BM1", "U89BI1_vs_U89BM1", "_padj-0.05_log2foldChange-1.csv", results_dir_with_log2FC_threshold)
get_n_save_unique_genes("BLS_DE_Results/02_Trt-vs-Mock/With_log2FC_threshold", "U17BI3_vs_U17BM3", "U89BI3_vs_U89BM3", "_padj-0.05_log2foldChange-1.csv", results_dir_with_log2FC_threshold)
get_n_save_unique_genes("BLS_DE_Results/02_Trt-vs-Mock/With_log2FC_threshold", "U89BI3_vs_U89BM3", "U17BI3_vs_U17BM3", "_padj-0.05_log2foldChange-1.csv", results_dir_with_log2FC_threshold)

######################################
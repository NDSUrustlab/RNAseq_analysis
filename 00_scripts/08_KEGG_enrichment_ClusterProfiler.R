library(clusterProfiler)
library(pathview)
library(tidyverse)
################################################################################
# create functions for the analysis
perform_kegg_enrichment <- function(ko_file, path, de_file, suffix, results_dir, threshold_padj, threshold_log2FoldChange) {
  # Convert thresholds to numeric
  threshold_padj <- as.numeric(threshold_padj)
  threshold_log2FoldChange <- as.numeric(threshold_log2FoldChange)
  
  # Read in KO data
  ko_data <- read.table(ko_file, header = FALSE, col.names = c("Gene", "KO"), fill = TRUE)
  
  # Get rid of duplicate transcripts in “universe” through gsub and regex
  ko_data[,1] <- sub("\\.[0-9]+$", "", ko_data[,1])
  ko_data <- distinct(ko_data)
  
  # Filter out entries without KO annotations
  ko_data <- ko_data[ko_data$KO != "", ]
  
  # Read in DE genes data
  de_genes <- read.csv(paste0(path, "/unique_DE_in_", de_file, suffix, ".csv"))
  
  # Ensure log2FoldChange and padj are numeric
  de_genes$log2FoldChange <- as.numeric(de_genes$log2FoldChange)
  de_genes$padj <- as.numeric(de_genes$padj)
  
  # Annotate upregulated and downregulated genes
  de_genes <- de_genes %>% mutate(regulation = case_when(
    log2FoldChange >= threshold_log2FoldChange & padj < threshold_padj ~ 'UP',
    log2FoldChange <= -threshold_log2FoldChange & padj < threshold_padj ~ 'DOWN',
    TRUE ~ 'NO'
  ))
  
  # Split the data into up- and down-regulated genes
  upregulated_genes <- de_genes %>% filter(regulation == 'UP')
  downregulated_genes <- de_genes %>% filter(regulation == 'DOWN')
  
  # Map DE genes to KO terms
  upregulated_genes_mapped <- upregulated_genes %>% inner_join(ko_data, by = "Gene")
  downregulated_genes_mapped <- downregulated_genes %>% inner_join(ko_data, by = "Gene")
  
  # Perform KEGG enrichment analysis on up-regulated genes
  if (nrow(upregulated_genes_mapped) > 0) {
    ekegg_up <- enrichKEGG(
      gene = upregulated_genes_mapped$KO,
      organism = "ko",
      keyType = "kegg",
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      minGSSize = 10,
      maxGSSize = 500,
      qvalueCutoff = 0.2,
      use_internal_data = FALSE
    )
  } else {
    warning("No upregulated genes found for KEGG enrichment analysis.")
  }
  
  # Perform KEGG enrichment analysis on down-regulated genes
  if (nrow(downregulated_genes_mapped) > 0) {
    ekegg_down <- enrichKEGG(
      gene = downregulated_genes_mapped$KO,
      organism = "ko",
      keyType = "kegg",
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      minGSSize = 10,
      maxGSSize = 500,
      qvalueCutoff = 0.2,
      use_internal_data = FALSE
    )
  } else {
    warning("No downregulated genes found for KEGG enrichment analysis.")
  }
  
  # Save results
  # Outfile name
  outfile_up <- paste0(results_dir, "/KEGG_enrichment_of_unique_DE_in_", de_file, "_padj-", threshold_padj, "_log2foldChange-", threshold_log2FoldChange, "_UPREGULATED")
  outfile_down <- paste0(results_dir, "/KEGG_enrichment_of_unique_DE_in_", de_file, "_padj-", threshold_padj, "_log2foldChange-", threshold_log2FoldChange, "_DOWNREGULATED")
  
  # Save in csv format
  if (!is.null(ekegg_up)) {
    write.csv(as.data.frame(ekegg_up), paste0(outfile_up, ".csv"), row.names = FALSE)
  }
  if (!is.null(ekegg_down)) {
    write.csv(as.data.frame(ekegg_down), paste0(outfile_down, ".csv"), row.names = FALSE)
  }
  
  # Plotting options
  # Example for up-regulated genes
  if (!is.null(ekegg_up)) {
    dotplot_up <- dotplot(ekegg_up, showCategory = 10)
    barplot_up <- barplot(ekegg_up, showCategory = 10)
    ggsave(paste0(outfile_up, "_dotplot.png"), plot = dotplot_up, dpi = 300, width = 10, height = 10)
    ggsave(paste0(outfile_up, "_barplot.png"), plot = barplot_up, dpi = 300, width = 10, height = 10)
  }
  
  # Example for down-regulated genes
  if (!is.null(ekegg_down)) {
    dotplot_down <- dotplot(ekegg_down, showCategory = 10)
    barplot_down <- barplot(ekegg_down, showCategory = 10)
    ggsave(paste0(outfile_down, "_dotplot.png"), plot = dotplot_down, dpi = 300, width = 10, height = 10)
    ggsave(paste0(outfile_down, "_barplot.png"), plot = barplot_down, dpi = 300, width = 10, height = 10)
  }
}
################################################################################
#Res v Sus
# Define the directory path
results_dir <- "BLS_Enrichment_Results/01_Res-vs-Sus/02_ORA_KEGG_enrichment"

# Check if the directory exists, and if not, create it
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

path <- "BLS_DE_Results/01_Res-vs-Sus/With_log2FC_threshold/01_unique_genes"
ko_file = "../../12_KEGG_analysis/02_with_renamed_genes/Aumb_BRAKER_RNA_seq_Iso_seq_genes_KO-terms.txt"

perform_kegg_enrichment(ko_file, path, "U17BI1_vs_U89BI1_compared_to_U17BM1_vs_U89BM1", "_padj-0.05_log2foldChange-1", results_dir, 0.05, 1)
perform_kegg_enrichment(ko_file, path, "U17BI3_vs_U89BI3_compared_to_U17BM3_vs_U89BM3", "_padj-0.05_log2foldChange-1", results_dir, 0.05, 1)
################################################################################
#Trt vs Mock
# Define the directory path
results_dir <- "BLS_Enrichment_Results/02_Trt-vs-Mock/02_ORA_KEGG_enrichment"

# Check if the directory exists, and if not, create it
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

path <- "BLS_DE_Results/02_Trt-vs-Mock/With_log2FC_threshold/01_unique_genes"
ko_file = "../../12_KEGG_analysis/02_with_renamed_genes/Aumb_BRAKER_RNA_seq_Iso_seq_genes_KO-terms.txt"
perform_kegg_enrichment(ko_file, path, "U17BI1_vs_U17BM1_compared_to_U89BI1_vs_U89BM1", "_padj-0.05_log2foldChange-1", results_dir, 0.05, 1)
perform_kegg_enrichment(ko_file, path, "U17BI3_vs_U17BM3_compared_to_U89BI3_vs_U89BM3", "_padj-0.05_log2foldChange-1", results_dir, 0.05, 1)
perform_kegg_enrichment(ko_file, path, "U89BI3_vs_U89BM3_compared_to_U17BI3_vs_U17BM3", "_padj-0.05_log2foldChange-1", results_dir, 0.05, 1)
################################################################################
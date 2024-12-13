### download the libraries
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install("EnhancedVolcano")

################################################################################

## load the libraries

library(DESeq2)
library(tidyverse)
library(EnhancedVolcano)
library(VennDiagram)
library(UpSetR)
library(pheatmap)


#################STEP: I#################
## preparing counting data

desired_order <- c("U17BI1R1", "U17BI1R2", "U17BI1R3", "U17BM1R1", "U17BM1R2", "U17BM1R3", "U17BI3R1", "U17BI3R2", "U17BI3R3", "U17BM3R1", "U17BM3R2", "U17BM3R3", "U89BI1R1", "U89BI1R2", "U89BI1R3", "U89BM1R1", "U89BM1R2", "U89BM1R3", "U89BI3R1", "U89BI3R2", "U89BI3R3", "U89BM3R1", "U89BM3R2", "U89BM3R3")

## reading the counts data
counts_data <- read.csv('bls_feature_count_for_DE.csv', sep="\t", row.names = 1)
counts_data <- counts_data %>% dplyr::select(all_of(desired_order))
head(counts_data)
dim(counts_data)

# reading meta data from sampleInfo file
colData <- read.csv('group_sampleInfo.csv', sep="\t", row.names = 1)
# Convert all columns to factors
colData$group <- as.factor(colData$group)
head(colData)

## making sure everything is right
# making sure the row names in colData matches to column names in counts_data
all(colnames(counts_data) %in% rownames(colData))

# are they in the same order?
all(colnames(counts_data) == rownames(colData))

#################STEP: II#################
## constructing a DESeqDataSet object

dds <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = colData,
                              design = ~ group)

dds

# pre-filtering: removing rows with low gene counts
# keeping rows that have at least 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds

#################STEP: III#################
## Running DESeq

dds <- DESeq(dds)
resultsNames(dds)
################################################################################
####################################################################
# Function to perform DESeq2 comparison and save DE genes to CSV
DE_comparison_and_save_csv_n_volcano_plot <- function(dds, group1, group2, results_dir, threshold_padj, threshold_log2FoldChange) {
  # Create the contrast and perform the analysis
  res <- results(dds, contrast=c("group", group1, group2))
  
  # create a df
  res_df <- as.data.frame(res)
  
  # Filter DE genes based on adjusted p-value and log2foldchange
  filtered_res <- res_df %>%
    filter(abs(log2FoldChange) > threshold_log2FoldChange & padj < threshold_padj) 
  
  # Add gene names as a column
  filtered_res$Gene <- rownames(filtered_res)
  
  # make the gene name column as the first column
  filtered_res <- filtered_res[, c("Gene", setdiff(names(filtered_res), "Gene"))]
  
  # Generate a filename based on group names, using the correct file path separator
  filename <- paste0(results_dir, "/", "DE_", gsub("_", "-", group1), "_vs_", gsub("_", "-", group2), "_padj-", threshold_padj, "_log2foldChange-", threshold_log2FoldChange)
  
  # Save to CSV file
  write.csv(filtered_res, file = paste0(filename, ".csv"), row.names = FALSE)
  
  ## now create a volcano plot
  # Convert to data frame with extra information for volcano plot
  res_df <- as.data.frame(res) %>%
    mutate(
      regulation = case_when(
        padj < threshold_padj & log2FoldChange > threshold_log2FoldChange ~ "Upregulated",
        padj < threshold_padj & log2FoldChange < -threshold_log2FoldChange ~ "Downregulated",
        TRUE ~ "Not Significant"
      )
    )
  
  # Calculate counts for upregulated and downregulated genes
  upregulated_count <- res_df %>%
    filter(regulation == "Upregulated") %>%
    nrow()
  
  downregulated_count <- res_df %>%
    filter(regulation == "Downregulated") %>%
    nrow()
  
  # Create custom key-value pairs for up and down regulated expression by log2 fold-change
  keyvals <- ifelse(res$log2FoldChange < -threshold_log2FoldChange & res$padj < threshold_padj, 'blue',
                    ifelse(res$log2FoldChange > threshold_log2FoldChange & res$padj < threshold_padj, 'red', 'grey'))
  
  # Ensure all NA values in keyvals are set to 'grey'
  keyvals[is.na(keyvals)] <- 'grey'
  
  # Create a vector for names matching the unique colors in keyvals
  names(keyvals) <- keyvals
  names(keyvals)[keyvals == 'red'] <- paste('Upregulated:', upregulated_count)
  names(keyvals)[keyvals == 'blue'] <- paste('Downregulated:', downregulated_count)
  names(keyvals)[keyvals == 'grey'] <- 'Not significant'
  
  # Combine the caption for the main plot
  main_caption <- paste0(
    'Log2 fold change cutoff: ', threshold_log2FoldChange, '\n',
    'Adjusted p-value cutoff: ', threshold_padj
  )
  
  # main title
  main_tile <- paste0(group1, " vs ", group2)
  
  # Create the Enhanced Volcano plot
  volcano_plot <- EnhancedVolcano(res,
                                  lab = NA,
                                  x = 'log2FoldChange',
                                  y = 'padj',
                                  title = main_tile,
                                  subtitle = NULL,
                                  xlab = bquote(~Log[2]~ 'fold change'),
                                  ylab = bquote(~-Log[10]~ 'adjP-value'),
                                  pCutoff = threshold_padj,
                                  FCcutoff = threshold_log2FoldChange,
                                  pointSize = 2,
                                  labSize = 4.5,
                                  colAlpha = 4/5,
                                  colCustom = keyvals,
                                  caption = main_caption,
                                  legendPosition = 'bottom',
                                  axisLabSize = 22,         # Increase axis label font size
                                  titleLabSize = 24,        # Increase title font size
                                  subtitleLabSize = 14,     # Increase subtitle font size
                                  captionLabSize = 18,      # Increase caption font size
                                  legendLabSize = 20)       # Increase legend text size
  
  # Save plot as PNG
  png(paste0(filename, ".png"), width = 10, height = 10, units = "in", res = 300)
  print(volcano_plot)
  dev.off()
  
}

# Perform DE comparisons and get DE genes
DE_comparison <- function(dds, group1, group2, threshold_padj, threshold_log2FoldChange) {
  # Create the contrast and perform the analysis
  res <- results(dds, contrast=c("group", group1, group2))
  
  # create a df
  res_df <- as.data.frame(res)
  
  # Filter DE genes based on adjusted p-value and log2foldchange
  filtered_res <- res_df %>%
    filter(abs(log2FoldChange) > threshold_log2FoldChange & padj < threshold_padj)
  
  # Return the list of DE genes
  return(filtered_res)
}

# Create a Upset Plot
create_upset_plot <- function(all_genes, results_dir, file_name = "UpSetPlot.png", text_scale = 1.5, point_size = 4, 
                              sets_x_label = "Number of DEGs", mainbar_y_label = "Shared DEGs") {
  
  # Convert to UpSetR data format
  upset_data <- fromList(all_genes)
  
  # Create UpSet plot
  upset_plot <- upset(
    upset_data,
    order.by = c("freq", "degree"),
    #main.bar.color = "#56B4E9",  # Light blue
    #sets.bar.color = "#E69F00",  # Reddish-orange
    text.scale = text_scale,
    point.size = point_size,
    sets.x.label = sets_x_label,
    mainbar.y.label = mainbar_y_label
  )
  
  # Save the UpSet plot
  png(paste0(results_dir, "/", file_name), width = 10, height = 6, units = "in", res = 1000)
  print(upset_plot)
  dev.off()
}

#################STEP: IV#################
# set the padj and log2foldchange thresholds
threshold_padj <- 0.05
threshold_log2FoldChange <- 1

#################STEP: V#################
# Resistant vs Susceptible (4 Comparison)

# Define the directory path
results_dir <- "BLS_DE_Results/01_Res-vs-Sus/With_log2FC_threshold"

# Check if the directory exists, and if not, create it
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

## With log2FC > 1
DE_comparison_and_save_csv_n_volcano_plot(dds, "U17BI1", "U89BI1", results_dir, threshold_padj, threshold_log2FoldChange)
DE_comparison_and_save_csv_n_volcano_plot(dds, "U17BI3", "U89BI3", results_dir, threshold_padj, threshold_log2FoldChange)
DE_comparison_and_save_csv_n_volcano_plot(dds, "U17BM1", "U89BM1", results_dir, threshold_padj, threshold_log2FoldChange)
DE_comparison_and_save_csv_n_volcano_plot(dds, "U17BM3", "U89BM3", results_dir, threshold_padj, threshold_log2FoldChange)

# Perform DE comparisons and get DE genes for upset Plot
DE_genes_U17BI1_vs_U89BI1 <- DE_comparison(dds, "U17BI1", "U89BI1", threshold_padj, threshold_log2FoldChange)
DE_genes_U17BI3_vs_U89BI3 <- DE_comparison(dds, "U17BI3", "U89BI3", threshold_padj, threshold_log2FoldChange)
DE_genes_U17BM1_vs_U89BM1 <- DE_comparison(dds, "U17BM1", "U89BM1", threshold_padj, threshold_log2FoldChange)
DE_genes_U17BM3_vs_U89BM3 <- DE_comparison(dds, "U17BM3", "U89BM3", threshold_padj, threshold_log2FoldChange)

# to make upset plot
# Combine all DE genes into a list
all_genes <- list(
  "U17BI1 vs U89BI1" = rownames(DE_genes_U17BI1_vs_U89BI1),
  "U17BI3 vs U89BI3" = rownames(DE_genes_U17BI3_vs_U89BI3),
  "U17BM1 vs U89BM1" = rownames(DE_genes_U17BM1_vs_U89BM1),
  "U17BM3 vs U89BM3" = rownames(DE_genes_U17BM3_vs_U89BM3)
)

create_upset_plot(all_genes, results_dir, file_name = "UpSetPlot_Res_VS_Sus.png")

####################################
#Treated vs Untreated (4 Comparisons)

# Define the directory path
results_dir <- "BLS_DE_Results/02_Trt-vs-Mock/With_log2FC_threshold"

# Check if the directory exists, and if not, create it
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

DE_comparison_and_save_csv_n_volcano_plot(dds, "U89BI1", "U89BM1", results_dir, threshold_padj, threshold_log2FoldChange)
DE_comparison_and_save_csv_n_volcano_plot(dds, "U89BI3", "U89BM3", results_dir, threshold_padj, threshold_log2FoldChange)
DE_comparison_and_save_csv_n_volcano_plot(dds, "U17BI1", "U17BM1", results_dir, threshold_padj, threshold_log2FoldChange)
DE_comparison_and_save_csv_n_volcano_plot(dds, "U17BI3", "U17BM3", results_dir, threshold_padj, threshold_log2FoldChange)


# Perform DE comparisons and get DE genes for upset Plot
DE_genes_U89BI1_vs_U89BM1 <- DE_comparison(dds, "U89BI1", "U89BM1", threshold_padj, threshold_log2FoldChange)
DE_genes_U89BI3_vs_U89BM3 <- DE_comparison(dds, "U89BI3", "U89BM3", threshold_padj, threshold_log2FoldChange)
DE_genes_U17BI1_vs_U17BM1 <- DE_comparison(dds, "U17BI1", "U17BM1", threshold_padj, threshold_log2FoldChange)
DE_genes_U17BI3_vs_U17BM3 <- DE_comparison(dds, "U17BI3", "U17BM3", threshold_padj, threshold_log2FoldChange)

# to make upset plot
# Combine all DE genes into a list
all_genes <- list(
  "U89BI1 vs U89BM1" = rownames(DE_genes_U89BI1_vs_U89BM1),
  "U89BI3 vs U89BM3" = rownames(DE_genes_U89BI3_vs_U89BM3),
  "U17BI1 vs U17BM1" = rownames(DE_genes_U17BI1_vs_U17BM1),
  "U17BI3 vs U17BM3" = rownames(DE_genes_U17BI3_vs_U17BM3)
)

create_upset_plot(all_genes, results_dir, file_name = "UpSetPlot_Trt_VS_Mock.png")

####################################
#Time points (4 Comparisons)

# Define the directory path
results_dir <- "BLS_DE_Results/03_Time1d-vs-3d/With_log2FC_threshold"

# Check if the directory exists, and if not, create it
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

DE_comparison_and_save_csv_n_volcano_plot(dds, "U17BI1", "U17BI3", results_dir, threshold_padj, threshold_log2FoldChange)
DE_comparison_and_save_csv_n_volcano_plot(dds, "U17BM1", "U17BM3", results_dir, threshold_padj, threshold_log2FoldChange)
DE_comparison_and_save_csv_n_volcano_plot(dds, "U89BI1", "U89BI3", results_dir, threshold_padj, threshold_log2FoldChange)
DE_comparison_and_save_csv_n_volcano_plot(dds, "U89BM1", "U89BM3", results_dir, threshold_padj, threshold_log2FoldChange)

# Perform DE comparisons and get DE genes for upset Plot
DE_genes_U17BI1_vs_U17BI3 <- DE_comparison(dds, "U17BI1", "U17BI3", threshold_padj, threshold_log2FoldChange)
DE_genes_U17BM1_vs_U17BM3 <- DE_comparison(dds, "U17BM1", "U17BM3", threshold_padj, threshold_log2FoldChange)
DE_genes_U89BI1_vs_U89BI3 <- DE_comparison(dds, "U89BI1", "U89BI3", threshold_padj, threshold_log2FoldChange)
DE_genes_U89BM1_vs_U89BM3 <- DE_comparison(dds, "U89BM1", "U89BM3", threshold_padj, threshold_log2FoldChange)

# to make upset plot
# Combine all DE genes into a list
all_genes <- list(
  "U17BI1 vs U17BI3" = rownames(DE_genes_U17BI1_vs_U17BI3),
  "U17BM1 vs U17BM3" = rownames(DE_genes_U17BM1_vs_U17BM3),
  "U89BI1 vs U89BI3" = rownames(DE_genes_U89BI1_vs_U89BI3),
  "U89BM1 vs U89BM3" = rownames(DE_genes_U89BM1_vs_U89BM3)
)

create_upset_plot(all_genes, results_dir, file_name = "UpSetPlot_1d_VS_3d.png")

####################################

# Function to perform DESeq2 comparison and count DE genes
DE_comparison_and_count <- function(dds, group1, group2, threshold_padj, threshold_log2FoldChange) {
  # Create the contrast and perform the analysis
  res <- results(dds, contrast=c("group", group1, group2))
  
  # Create a data frame
  res_df <- as.data.frame(res)
  
  # Filter DE genes based on adjusted p-value and log2foldchange
  upregulated <- res_df %>%
    filter(log2FoldChange > threshold_log2FoldChange & padj < threshold_padj) %>%
    nrow()
  downregulated <- res_df %>%
    filter(log2FoldChange < -threshold_log2FoldChange & padj < threshold_padj) %>%
    nrow()
  
  # Return the counts
  return(c(upregulated, downregulated))
}

# Define thresholds
threshold_padj <- 0.05
threshold_log2FoldChange <- 1

# List of comparisons
comparisons <- list(
  c("U17BI1", "U89BI1"),
  c("U17BI3", "U89BI3"),
  c("U17BM1", "U89BM1"),
  c("U17BM3", "U89BM3"),
  c("U17BI1", "U17BM1"),
  c("U17BI3", "U17BM3"),
  c("U89BI1", "U89BM1"),
  c("U89BI3", "U89BM3"),
  c("U17BI1", "U17BI3"),
  c("U17BM1", "U17BM3"),
  c("U89BI1", "U89BI3"),
  c("U89BM1", "U89BM3")
)

# Initialize a data frame to store the results
results_df <- data.frame(
  contrast_name = character(),
  up = integer(),
  down = integer(),
  stringsAsFactors = FALSE
)

# Run DESeq2 comparisons and store results
for (comp in comparisons) {
  contrast_name <- paste(comp[1], "vs", comp[2])
  counts <- DE_comparison_and_count(dds, comp[1], comp[2], threshold_padj, threshold_log2FoldChange)
  results_df <- rbind(results_df, data.frame(contrast_name = contrast_name, up = counts[1], down = counts[2]))
}

head(results_df)

# With log2FC threshold
# Save the results to a CSV file
write.csv(results_df, "BLS_DE_Results/DE_genes_counts_padj-0.05_log2FC-1.csv", row.names = FALSE)

# with log2FC threshold
# Convert contrast_name to factor to maintain order in plotting
results_df$contrast_name <- factor(results_df$contrast_name, levels = results_df$contrast_name)

plot <- ggplot(results_df, aes(x = contrast_name)) +
  geom_bar(aes(y = up, fill = "Upregulated"), stat = "identity", position = "dodge") +
  geom_bar(aes(y = -down, fill = "Downregulated"), stat = "identity", position = "dodge") +
  geom_text(aes(y = up, label = up), vjust = ifelse(results_df$up >= 0, 1.5, 1), size = 3, position = position_dodge(width = 0.9), color = "black") +
  geom_text(aes(y = -down, label = down), vjust = ifelse(results_df$down >= 0, -0.5, 1), size = 3, position = position_dodge(width = 0.9), color = "white") +
  scale_y_continuous(labels = function(x) abs(x), breaks = c(-6000, -4000, -2000, 0, 2000, 4000, 6000)) +
  labs(x = "", y = "Number of genes", fill = "") +
  #ggtitle("Counts of Upregulated and Downregulated Genes by Contrast") +
  theme_minimal() +
  scale_fill_manual(values = c("Upregulated" = "#E87538", "Downregulated" = "#004D40")) +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text = element_text(color = "black", size = 14),
        legend.text = element_text(size = 14),
        legend.position = "bottom") 
  #annotate("text", x = length(unique(results_df$contrast_name))/2, y = -5000, 
           #label = paste("adj. p-value <", threshold_padj, "\n log2FoldChange >", threshold_log2FoldChange), 
           #size = 6, hjust = 1.5, vjust = 0.5)

print(plot)

# Save the plot
ggsave("BLS_DE_Results/DE_count_grouped_bar_plot_with_log2FC_threshold-1.png", plot, width = 9, height = 7, units = "in", dpi = 300)

############################################################################################################
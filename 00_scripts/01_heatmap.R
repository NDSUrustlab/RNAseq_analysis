## load the libraries

library(DESeq2)
library(tidyverse)
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

####################################################################
# Extract normalized counts
normalized_counts <- counts(dds, normalized = TRUE)

# Select differentially expressed genes
res <- results(dds)
sig_genes <- rownames(subset(res, padj < 0.05 & abs(log2FoldChange) > 1))
sig_counts <- normalized_counts[sig_genes, ]

# Log transformation
log_counts <- log2(sig_counts + 1)

# Create sample annotations
annotation_col <- as.data.frame(colData)
# Add two spaces to the second column values
annotation_col$group <- paste0(annotation_col$group, "  ")

################################################################################
# Define a custom color palette using RColorBrewer
custom_colors <- colorRampPalette(brewer.pal(11, "RdBu"))(100)

# Generate the heatmap with custom colors
p <- pheatmap(log_counts, 
              scale = "row", 
              cluster_rows=TRUE,
              cluster_cols=TRUE,
              clustering_distance_rows = "euclidean", 
              clustering_distance_cols = "euclidean", 
              clustering_method = "complete", 
              annotation_col = annotation_col,
              show_rownames = FALSE,  # Hide row names (gene names)
              show_colnames = TRUE,
              color = custom_colors,  # Use custom colors
              fontsize_col = 12,      # Adjust the font size for column labels
              width = 12,             # Adjust the width of the plot
              height = 15)            # Adjust the height of the plot

# Save the plot with 300 DPI resolution
png("BLS_DE_Results/gene-expression_pheatmap_plot.png", width = 10, height = 6, units = "in", res = 300)
print(p)  # Print the pheatmap plot
dev.off()  # Close the PNG device
################################################################################







# Generate the heatmap with multiple annotations
pheatmap(log_counts, 
         scale = "row", 
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean", 
         clustering_method = "complete", 
         annotation_col = annotation_col,
         show_rownames = FALSE,  # Hide row names (gene names)
         show_colnames = TRUE,
         fontsize_col = 8,       # Adjust the font size for column labels
         width = 10,             # Adjust the width of the plot
         height = 15)            # Adjust the height of the plot




library(RColorBrewer)

# Define a custom color palette
custom_colors <- colorRampPalette(c("blue", "white", "red"))(100)

# Generate the heatmap with custom colors
pheatmap(log_counts, 
         scale = "row", 
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean", 
         clustering_method = "complete", 
         annotation_col = annotation_col,
         show_rownames = FALSE,  # Hide row names (gene names)
         show_colnames = TRUE,
         color = custom_colors,  # Use custom colors
         fontsize_col = 8,       # Adjust the font size for column labels
         width = 10,             # Adjust the width of the plot
         height = 15)            # Adjust the height of the plot





# Install and load the viridis package if necessary
if (!requireNamespace("viridis", quietly = TRUE))
  install.packages("viridis")
library(viridis)

# Generate the heatmap with viridis colors
pheatmap(log_counts, 
         scale = "row", 
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean", 
         clustering_method = "complete", 
         annotation_col = annotation_col,
         show_rownames = FALSE,  # Hide row names (gene names)
         show_colnames = TRUE,
         color = viridis,   # Use viridis colors
         fontsize_col = 8,       # Adjust the font size for column labels
         width = 10,             # Adjust the width of the plot
         height = 15)            # Adjust the height of the plot




# Define a custom color palette using RColorBrewer
custom_colors <- colorRampPalette(brewer.pal(11, "RdYlBu"))(100)

# Generate the heatmap with custom colors
pheatmap(log_counts, 
         scale = "row", 
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean", 
         clustering_method = "complete", 
         annotation_col = annotation_col,
         show_rownames = FALSE,  # Hide row names (gene names)
         show_colnames = TRUE,
         color = custom_colors,  # Use custom colors
         fontsize_col = 8,       # Adjust the font size for column labels
         width = 10,             # Adjust the width of the plot
         height = 15)            # Adjust the height of the plot








pheatmap(assay(vst(dds))[top_genes,], cluster_rows=TRUE, show_rownames=FALSE,
         cluster_cols=TRUE, annotation_col=as.data.frame(colData(dds)))



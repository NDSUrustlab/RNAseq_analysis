## load the libraries

library(DESeq2)
library(tidyverse)

##################################
## Preparing the data
#reading the metadata
sample_info <- read.csv('sample_info.csv', row.names = 1)
# Convert all columns to factors
sample_info[, c("genotype", "time", "condition")] <- lapply(sample_info[, c("genotype", "time", "condition")], as.factor)
head(sample_info)

#reading the counts data
counts <- read.csv('bls_feature_count_for_DE.csv', sep="\t", row.names = 1)
head(counts)

#order of the samples
desired_order <- c("U17BI1R1", "U17BI1R2", "U17BI1R3", "U17BM1R1", "U17BM1R2", "U17BM1R3", "U17BI3R1", "U17BI3R2", "U17BI3R3", "U17BM3R1", "U17BM3R2", "U17BM3R3", "U89BI1R1", "U89BI1R2", "U89BI1R3", "U89BM1R1", "U89BM1R2", "U89BM1R3", "U89BI3R1", "U89BI3R2", "U89BI3R3", "U89BM3R1", "U89BM3R2", "U89BM3R3")

#reordering the column names to match to the sampleinfo rownames
counts <- counts %>% dplyr::select(desired_order)

#checking if colnames in counts is similar with rownames in sample_info and if the names are in order
all(colnames(counts) %in% rownames(sample_info))
all(colnames(counts) == rownames(sample_info))

##################################
## Create a DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = sample_info,
                              design = ~ genotype + time + condition + genotype:time + genotype:condition + condition:time + genotype:condition:time)
dds

# pre-filtering: removing rows with low gene counts
# keeping rows that have at least 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds

##################################
## Run DESeq
dds$genotype <- relevel(dds$genotype, ref = "Susceptible")
dds$time <- relevel(dds$time, ref = "3d")
# dds$condition <- relevel(dds$condition, ref = "treated")

dds <- DESeq(dds)
resultsNames(dds)

##################################
# Define the directory path
# Check if the directory exists, and if not, create it
if (!dir.exists("BLS_DE_Results")) {
  dir.create("BLS_DE_Results", recursive = TRUE)
}
##################################

##################################
# PCA Plot
vsd <- vst(dds)
plotPCA(vsd, intgroup=c("genotype", "condition", "time"))

pcaData <- plotPCA(vsd, intgroup = c("genotype", "time", "condition"), returnData = TRUE)
explained_variance <- round(100 * attr(pcaData, "percentVar"))

PCAplot <- ggplot(pcaData, 
                  aes(PC1, PC2, color = condition, shape = time)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", explained_variance[1], "% variance")) +
  ylab(paste0("PC2: ", explained_variance[2], "% variance")) + 
  coord_fixed(ratio = 1) +
  #ggtitle("PCA Score Plot - Differentiating PI 554389 (Resistant) and 554417 (Susceptible) Over Time") +
  theme(
    text = element_text(size = 18),  # Adjust the overall text size
    axis.text = element_text(size = 18),  # Adjust axis text size
    legend.title = element_text(size = 20, face = "bold"),  # Adjust legend title size and style
    legend.text = element_text(size = 14),  # Adjust legend text size
    legend.position = "bottom",  # Change legend position to the bottom
    legend.box = "horizontal",  # Display legend items horizontally
    panel.background = element_blank(),     # Remove panel background
    panel.border = element_rect(color = "#000000", fill = NA)  # Add border around the plot
  )+
  stat_ellipse(aes(group=genotype), alpha = 0.15, geom = "polygon", linetype =1, type = "norm", color="black") +
  annotate(
    "text",
    x = 25, 
    y = 70,
    label = "PI 554389",
    size = 8
  ) +
  annotate(
    "text",
    x = -42, 
    y = 60,
    label = "PI 554417",
    size = 8
  ) +
  scale_shape_manual(name = "DAI",  # Change legend title to DAI
                     breaks = c("1d", "3d"),  # Specify breaks for the legend
                     values = c("1d" = 5, "3d" = 1),  # Assign shapes to corresponding levels
                     labels = c("1d", "3d")) +  # Specify labels for legend
  guides(
    color = guide_legend(title = "Condition"),  # Customize color legend title
    shape = guide_legend(title = "DAI")  # Customize shape legend title
  )

print(PCAplot)
ggsave("BLS_DE_Results/PCAplot.png", PCAplot, dpi=1000, width=7, height=10)

#################################
# Sample Distance Heatmap
vsd_hp <- vst(dds, blind = FALSE)

# Calculate sample-to-sample distances
sampleDists <- dist(t(assay(vsd_hp)))

# Create a matrix of the distances
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(dds)
colnames(sampleDistMatrix) <- colnames(dds)

# Generate the heatmap
HeatMap <- pheatmap(sampleDistMatrix,
             clustering_distance_rows = sampleDists,
             clustering_distance_cols = sampleDists,
             display_numbers = FALSE, # optionally show distances
             fontsize = 14  # Adjust fontsize for axis text
             #main = "Sample Distance Heatmap"
             )

ggsave("BLS_DE_Results/Heatmap_sample_distance.png", HeatMap, dpi=1000, width=12, height=10)

#################################
# Plot dispersion estimates

# Extract the mean normalized counts
means <- rowMeans(counts(dds, normalized = TRUE))

# Extract dispersion estimates
gene_disp <- mcols(dds)$dispGeneEst
fitted_disp <- mcols(dds)$dispFit
final_disp <- mcols(dds)$dispersion

# Create a data frame for plotting
disp_data <- data.frame(
  mean = means,
  gene_disp = gene_disp,
  fitted_disp = fitted_disp,
  final_disp = final_disp
)

plotDispEstsEnhanced <- function(disp_data, main = "Dispersion Estimates") {
  ggplot(disp_data, aes(x = mean)) +
    geom_point(aes(y = gene_disp), color = "black", size = 1, alpha = 0.5) +
    geom_point(aes(y = final_disp), color = "blue", size = 1, alpha = 0.5) +
    geom_line(aes(y = fitted_disp), color = "red", size = 1) +
    scale_x_log10() +
    scale_y_log10() +
    labs(
      #title = main,
      x = "Mean of normalized counts",
      y = "Dispersion"
    ) +
    theme_minimal() +
    theme(
      #plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      axis.text = element_text(size = 10)
    )
}

# Use the enhanced plotting function
p <- plotDispEstsEnhanced(disp_data)

ggsave("BLS_DE_Results/Dispersion_estimate_plot.png", p, dpi=1000, width=6, height=6)
#################################
library(clusterProfiler)
library(tidyverse)
library(topGO)
library(Rgraphviz)


##install.packages("/mmfs1/scratch/jatinder.singh/RNAseq_aumb_parents/16_DE_n_enrichment/02_bls_DE/org.Aumbellulata.eg.db", repos = NULL)
# Uninstall the org.Aumbellulata.eg.db package
#remove.packages("org.Aumbellulata.eg.db")
library(org.Aumbellulata.eg.db)

################################################################################
# create functions for the analysis

perform_go_enrichment <- function(path, file, suffix, results_dir, threshold_padj, threshold_log2FoldChange) {
  # Convert thresholds to numeric
  threshold_padj <- as.numeric(threshold_padj)
  threshold_log2FoldChange <- as.numeric(threshold_log2FoldChange)
  
  # Read in the data
  de_genes <- read.csv(paste0(path, "/unique_DE_in_", file, suffix, ".csv"))
  
  # Debugging step: Check the structure of de_genes
  print(head(de_genes))
  print(str(de_genes))
  
  # Ensure log2FoldChange and padj are numeric
  de_genes$log2FoldChange <- as.numeric(de_genes$log2FoldChange)
  de_genes$padj <- as.numeric(de_genes$padj)
  
  # Annotate upregulated and downregulated genes
  de_genes <- de_genes %>% mutate(regulation = case_when(
    log2FoldChange >= threshold_log2FoldChange & padj < threshold_padj ~ 'UP',
    log2FoldChange <= -threshold_log2FoldChange & padj < threshold_padj ~ 'DOWN',
    TRUE ~ 'NO'
  ))
  
  # Debugging step: Check the regulation column
  print(table(de_genes$regulation))
  
  # Split the data into up- and down-regulated genes
  upregulated_genes <- de_genes %>% filter(regulation == 'UP')
  downregulated_genes <- de_genes %>% filter(regulation == 'DOWN')
  
  # Debugging step: Ensure there are genes in each category
  print(paste("Number of upregulated genes:", nrow(upregulated_genes)))
  print(paste("Number of downregulated genes:", nrow(downregulated_genes)))
  
  # Perform GO enrichment analysis on up-regulated genes
  if (nrow(upregulated_genes) > 0) {
    up_gene_list <- upregulated_genes$Gene
    ego_up <- enrichGO(
      gene = up_gene_list,
      OrgDb = org.Aumbellulata.eg.db,  # Replace with the correct organism database
      keyType = "GID",
      ont = "all",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2,
      readable = FALSE
    )
  } else {
    ego_up <- NULL
    warning("No upregulated genes found for GO enrichment analysis.")
  }
  
  # Perform GO enrichment analysis on down-regulated genes
  if (nrow(downregulated_genes) > 0) {
    down_gene_list <- downregulated_genes$Gene
    ego_down <- enrichGO(
      gene = down_gene_list,
      OrgDb = org.Aumbellulata.eg.db,  # Replace with the correct organism database
      keyType = "GID",
      ont = "all",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2,
      readable = FALSE
    )
  } else {
    ego_down <- NULL
    warning("No downregulated genes found for GO enrichment analysis.")
  }
  
  # Outfile name
  outfile_up <- paste0(results_dir, "/GO_enrichment_of_unique_DE_in_", file, "_padj-", threshold_padj, "_log2foldChange-", threshold_log2FoldChange, "_UPREGULATED")
  outfile_down <- paste0(results_dir, "/GO_enrichment_of_unique_DE_in_", file, "_padj-", threshold_padj, "_log2foldChange-", threshold_log2FoldChange, "_DOWNREGULATED")
  
  # Save in csv format
  if (!is.null(ego_up)) {
    write.csv(as.data.frame(ego_up), paste0(outfile_up, ".csv"), row.names = FALSE)
  }
  if (!is.null(ego_down)) {
    write.csv(as.data.frame(ego_down), paste0(outfile_down, ".csv"), row.names = FALSE)
  }
  
  # Plotting options
  # Example for up-regulated genes
  if (!is.null(ego_up)) {
    dotplot_up <- dotplot(ego_up, showCategory = 10)
    barplot_up <- barplot(ego_up, showCategory = 10)
    ggsave(paste0(outfile_up, "_dotplot.png"), plot = dotplot_up, dpi = 300, width = 10, height = 10)
    ggsave(paste0(outfile_up, "_barplot.png"), plot = barplot_up, dpi = 300, width = 10, height = 10)
  }
  
  # Example for down-regulated genes
  if (!is.null(ego_down)) {
    dotplot_down <- dotplot(ego_down, showCategory = 10)
    barplot_down <- barplot(ego_down, showCategory = 10)
    ggsave(paste0(outfile_down, "_dotplot.png"), plot = dotplot_down, dpi = 300, width = 10, height = 10)
    ggsave(paste0(outfile_down, "_barplot.png"), plot = barplot_down, dpi = 300, width = 10, height = 10)
  }
}

################################################################################
#Res v Sus
# Define the directory path
results_dir <- "BLS_Enrichment_Results/01_Res-vs-Sus/01_ORA_GO_enrichment"

# Check if the directory exists, and if not, create it
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

path <- "BLS_DE_Results/01_Res-vs-Sus/With_log2FC_threshold/01_unique_genes"

perform_go_enrichment(path, "U17BI1_vs_U89BI1_compared_to_U17BM1_vs_U89BM1", "_padj-0.05_log2foldChange-1", results_dir, 0.05, 1)
perform_go_enrichment(path, "U17BI3_vs_U89BI3_compared_to_U17BM3_vs_U89BM3", "_padj-0.05_log2foldChange-1", results_dir, 0.05, 1)

################################################################################
#Trt vs Mock
# Define the directory path
results_dir <- "BLS_Enrichment_Results/02_Trt-vs-Mock/01_ORA_GO_enrichment"

# Check if the directory exists, and if not, create it
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

path <- "BLS_DE_Results/02_Trt-vs-Mock/With_log2FC_threshold/01_unique_genes"

perform_go_enrichment(path, "U17BI1_vs_U17BM1_compared_to_U89BI1_vs_U89BM1", "_padj-0.05_log2foldChange-1", results_dir, 0.05, 1)
perform_go_enrichment(path, "U17BI3_vs_U17BM3_compared_to_U89BI3_vs_U89BM3", "_padj-0.05_log2foldChange-1", results_dir, 0.05, 1)
perform_go_enrichment(path, "U89BI3_vs_U89BM3_compared_to_U17BI3_vs_U17BM3", "_padj-0.05_log2foldChange-1", results_dir, 0.05, 1)

################################################################################

# run enrichment analysis on common genes of unique trt 1d (7036 genes) and unique PI55417 at 1d (7979)

file1 <- read.csv("BLS_DE_Results/01_Res-vs-Sus/With_log2FC_threshold/01_unique_genes/unique_DE_in_U17BI1_vs_U89BI1_compared_to_U17BM1_vs_U89BM1_padj-0.05_log2foldChange-1.csv", row.names = 1)
nrow(file1)
file2 <- read.csv("BLS_DE_Results/02_Trt-vs-Mock/With_log2FC_threshold/01_unique_genes/unique_DE_in_U17BI1_vs_U17BM1_compared_to_U89BI1_vs_U89BM1_padj-0.05_log2foldChange-1.csv", row.names = 1)
nrow(file2)

cg <- Reduce(intersect, list(rownames(file1),
                                      rownames(file2)))
length(cg)


cg_from_file1 <- file1[cg, ]
cg_from_file1$Gene <- rownames(cg_from_file1)
nrow(cg_from_file1)
head(cg_from_file1)

cg_from_file2 <- file2[cg, ]
cg_from_file2$Gene <- rownames(cg_from_file2)
nrow(cg_from_file2)
head(cg_from_file2)


#Trt vs Mock
# Define the directory path
results_dir <- "BLS_Enrichment_Results/03_common_genes_3791/01_ORA_GO_enrichment"

# Check if the directory exists, and if not, create it
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

perform_go_enrichment(cg_from_file1, "cg_from_file1", results_dir, 0.05, 1)
perform_go_enrichment(cg_from_file2, "cg_from_file2", results_dir, 0.05, 1)

perform_go_enrichment <- function(de_genes, file_name, results_dir, threshold_padj, threshold_log2FoldChange) {
  # Convert thresholds to numeric
  threshold_padj <- as.numeric(threshold_padj)
  threshold_log2FoldChange <- as.numeric(threshold_log2FoldChange)
  
  # Debugging step: Check the structure of de_genes
  print(head(de_genes))
  print(str(de_genes))
  
  # Ensure log2FoldChange and padj are numeric
  de_genes$log2FoldChange <- as.numeric(de_genes$log2FoldChange)
  de_genes$padj <- as.numeric(de_genes$padj)
  
  # Annotate upregulated and downregulated genes
  de_genes <- de_genes %>% mutate(regulation = case_when(
    log2FoldChange >= threshold_log2FoldChange & padj < threshold_padj ~ 'UP',
    log2FoldChange <= -threshold_log2FoldChange & padj < threshold_padj ~ 'DOWN',
    TRUE ~ 'NO'
  ))
  
  # Debugging step: Check the regulation column
  print(table(de_genes$regulation))
  
  # Split the data into up- and down-regulated genes
  upregulated_genes <- de_genes %>% filter(regulation == 'UP')
  downregulated_genes <- de_genes %>% filter(regulation == 'DOWN')
  
  # Debugging step: Ensure there are genes in each category
  print(paste("Number of upregulated genes:", nrow(upregulated_genes)))
  print(paste("Number of downregulated genes:", nrow(downregulated_genes)))
  
  # Perform GO enrichment analysis on up-regulated genes
  if (nrow(upregulated_genes) > 0) {
    up_gene_list <- upregulated_genes$Gene
    ego_up <- enrichGO(
      gene = up_gene_list,
      OrgDb = org.Aumbell.eg.db,  # Replace with the correct organism database
      keyType = "GID",
      ont = "all",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2,
      readable = FALSE
    )
  } else {
    ego_up <- NULL
    warning("No upregulated genes found for GO enrichment analysis.")
  }
  
  # Perform GO enrichment analysis on down-regulated genes
  if (nrow(downregulated_genes) > 0) {
    down_gene_list <- downregulated_genes$Gene
    ego_down <- enrichGO(
      gene = down_gene_list,
      OrgDb = org.Aumbell.eg.db,  # Replace with the correct organism database
      keyType = "GID",
      ont = "all",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2,
      readable = FALSE
    )
  } else {
    ego_down <- NULL
    warning("No downregulated genes found for GO enrichment analysis.")
  }
  
  # Outfile name
  outfile_up <- paste0(results_dir, "/GO_enrichment_of_unique_DE_in_", file_name, "_padj-", threshold_padj, "_log2foldChange-", threshold_log2FoldChange, "_UPREGULATED")
  outfile_down <- paste0(results_dir, "/GO_enrichment_of_unique_DE_in_", file_name, "_padj-", threshold_padj, "_log2foldChange-", threshold_log2FoldChange, "_DOWNREGULATED")
  
  # Save in csv format
  if (!is.null(ego_up)) {
    write.csv(as.data.frame(ego_up), paste0(outfile_up, ".csv"), row.names = FALSE)
  }
  if (!is.null(ego_down)) {
    write.csv(as.data.frame(ego_down), paste0(outfile_down, ".csv"), row.names = FALSE)
  }
  
  # Plotting options
  # Example for up-regulated genes
  if (!is.null(ego_up)) {
    dotplot_up <- dotplot(ego_up, showCategory = 10)
    barplot_up <- barplot(ego_up, showCategory = 10)
    ggsave(paste0(outfile_up, "_dotplot.png"), plot = dotplot_up, dpi = 300, width = 10, height = 10)
    ggsave(paste0(outfile_up, "_barplot.png"), plot = barplot_up, dpi = 300, width = 10, height = 10)
  }
  
  # Example for down-regulated genes
  if (!is.null(ego_down)) {
    dotplot_down <- dotplot(ego_down, showCategory = 10)
    barplot_down <- barplot(ego_down, showCategory = 10)
    ggsave(paste0(outfile_down, "_dotplot.png"), plot = dotplot_down, dpi = 300, width = 10, height = 10)
    ggsave(paste0(outfile_down, "_barplot.png"), plot = barplot_down, dpi = 300, width = 10, height = 10)
  }
}







common_genes <- inner_join(file1, file2, by="Gene")
nrow(common_genes)

# Define the directory path
results_dir <- "BLS_Enrichment_Results/03_Common_genes_3791/01_ORA_GO_enrichment"

# Check if the directory exists, and if not, create it
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

path <- "BLS_DE_Results/01_Res-vs-Sus/With_log2FC_threshold/01_unique_genes"

perform_go_enrichment(path, "U17BI1_vs_U89BI1_compared_to_U17BM1_vs_U89BM1", "_padj-0.05_log2foldChange-1", results_dir, 0.05, 1)



################################################################################
file <- file %>% mutate(regulation = case_when(
  log2FoldChange > 1 & padj < 0.05 ~ 'UP',
  log2FoldChange < -1 & padj < 0.05 ~ 'DOWN',
  padj >= 0.05 ~ 'NO'
))

# Split the data into up- and down-regulated genes
upregulated_genes <- file %>% filter(regulation == 'UP')
downregulated_genes <- file %>% filter(regulation == 'DOWN')
nrow(upregulated_genes)
nrow(downregulated_genes)
threshold_log2FoldChange
####################
#checking

tmp <- read.csv("BLS_DE_Results/01_Res-vs-Sus/Without_log2FC_threshold/01_unique_genes/unique_DE_in_U17BI1_vs_U89BI1_compared_to_U17BM1_vs_U89BM1_padj-0.05_Without_log2foldChange-threshold.csv")
nrow(tmp)
head(tmp)

# Assuming `tmp` is your data frame and it contains a column `log2FoldChange`
tmp <- tmp %>% 
  mutate(regulation = case_when(
    log2FoldChange >= 1 ~ "upregulated",
    log2FoldChange <= -1 ~ "downregulated",
    TRUE ~ NA_character_
  ))

# Convert regulation to factor to ensure correct handling in count
tmp$regulation <- as.factor(tmp$regulation)

table(tmp$regulation)
# Count the number of each category
regulation_counts <- tmp %>% 
  count(regulation)

print(regulation_counts)

###################
res <- results(dds, contrast=c("group", "U17BI1", "U89BI1"))

# create a df
res_df <- as.data.frame(res)
unfiltered_res <- res_df %>%
  filter(padj < 0.05)

nrow(unfiltered_res)

# Filter DE genes based on adjusted p-value and log2foldchange
filtered_res <- res_df %>%
  filter(abs(log2FoldChange) >= 1 & padj < 0.05)

nrow(filtered_res)

####
res_1 <- results(dds, contrast=c("group", "U17BM1", "U89BM1"))


# create a df
res_df_1 <- as.data.frame(res_1)
unfiltered_res_1 <- res_df_1 %>%
  filter(padj < 0.05)

nrow(unfiltered_res_1)

# Filter DE genes based on adjusted p-value and log2foldchange
filtered_res_1 <- res_df_1 %>%
  filter(abs(log2FoldChange) >= 1 & padj < 0.05)

nrow(filtered_res_1)
####
# Convert row names to a column
filtered_res <- filtered_res %>% rownames_to_column(var = "row_name")
filtered_res <- filtered_res %>% 
  mutate(regulation = ifelse(log2FoldChange >= 1, "upregulated", ifelse(log2FoldChange <= -1, "downregulated", NA)))

filtered_res_1 <- filtered_res_1 %>% rownames_to_column(var = "row_name")
filtered_res_1 <- filtered_res_1 %>% 
  mutate(regulation = ifelse(log2FoldChange >= 1, "upregulated", ifelse(log2FoldChange <= -1, "downregulated", NA)))

###
nrow(anti_join(filtered_res, filtered_res_1, by = "row_name"))


# Convert row names to a column
unfiltered_res <- unfiltered_res %>% rownames_to_column(var = "row_name")
unfiltered_res <- unfiltered_res %>% 
  mutate(regulation = ifelse(log2FoldChange >= 1, "upregulated", ifelse(log2FoldChange <= -1, "downregulated", NA)))

unfiltered_res_1 <- unfiltered_res_1 %>% rownames_to_column(var = "row_name")
unfiltered_res_1 <- unfiltered_res_1 %>% 
  mutate(regulation = ifelse(log2FoldChange >= 1, "upregulated", ifelse(log2FoldChange <= -1, "downregulated", NA)))

###
nrow(anti_join(unfiltered_res, unfiltered_res_1, by = "row_name"))
table(filtered_res$regulation)
table(filtered_res_1$regulation)
table(unfiltered_res$regulation)
table(unfiltered_res_1$regulation)


###
comp1 <- anti_join(filtered_res, filtered_res_1, by = "row_name")
table(comp1$regulation)

comp2 <- anti_join(unfiltered_res, unfiltered_res_1, by = "row_name")
table(comp2$regulation)

####################

de_genes <- read.csv("BLS_DE_Results/01_Res-vs-Sus/With_log2FC_threshold/01_unique_genes/unique_DE_in_U17BI1_vs_U89BI1_compared_to_U17BM1_vs_U89BM1_padj-0.05_log2foldChange-1.csv")
gene_list <- de_genes$Gene

de_genes_17VS89_trt_1 <- de_genes_17VS89_trt_1 %>% 
  mutate(regulation = ifelse(log2FoldChange >= 1, "upregulated", ifelse(log2FoldChange <= -1, "downregulated", NA)))



# Perform GO enrichment analysis
ego <- enrichGO(
  gene = gene_list,
  OrgDb = org.Aumbell.eg.db,
  keyType = "GID",  # Replace with the appropriate key type if different
  ont = "all",       # Ontology to use: "BP" (Biological Process), "MF" (Molecular Function), or "CC" (Cellular Component)
  pAdjustMethod = "BH", # p-value adjustment method
  pvalueCutoff = 0.05,  # p-value cutoff
  qvalueCutoff = 0.2,   # q-value cutoff
  readable = FALSE       
)

head(ego)

# Calculate term similarity matrix
ego_sim <- pairwise_termsim(ego)

# Plotting options
# 1. Dotplot
dotplot(ego, showCategory = 20) + ggtitle("Dotplot for GO Enrichment Analysis")

# 2. Barplot
barplot(ego, showCategory = 20) + ggtitle("Barplot for GO Enrichment Analysis")

# 3. Enrichment Map
emapplot(ego_sim) + ggtitle("Enrichment Map for GO Enrichment Analysis")

# 4. GO Graph
goplot(ego) + ggtitle("GO Graph for GO Enrichment Analysis")

# 5. Network plot
cnetplot(ego, showCategory = 5) + ggtitle("Network Plot for GO Enrichment Analysis")


################################################################################
all_genes_ontologies <- read.csv("gene_ontologies_all_genes_v2.csv")

term2gene <- all_genes_ontologies %>%
  dplyr::select(go, Gene)

# Mutate to add rank based on log2FoldChange and arrange by log2FoldChange
de_genes_df <- de_genes %>%
  arrange(desc(log2FoldChange)) %>%
  mutate(rank = rank(log2FoldChange, ties.method = "random")) %>%
  arrange(desc(rank))

d <- de_genes_df %>%
  dplyr::select("Gene", "rank")


geneList <- d[, 2]
geneList <- sort(geneList, decreasing = TRUE)
names(geneList) <- as.character(d[, 1])


# Perform GSEA
gsea_result <- GSEA(geneList = geneList,
                    exponent = 1,
                    minGSSize = 10,
                    maxGSSize = 500,
                    eps = 1e-10,
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "BH",
                    TERM2GENE = term2gene,
                    TERM2NAME = NA,
                    verbose = TRUE,
                    seed = FALSE,
                    by = "fgsea")
dotplot(gsea_result) + ggtitle("GSEA")

gsea_result <- gseGO(
  geneList = geneList,
  OrgDb = org.Aumbell.eg.db, # Replace with your OrgDb package if different
  keyType = "GID", # Replace with the appropriate key type if different
  ont = "ALL",          # Ontology: "BP", "MF", "CC", or "ALL"
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  verbose = TRUE
)

# Visualize the results
dotplot(gsea_result, showCategory = 10)
barplot(gsea_result, showCategory = 20)
emapplot(gsea_result)
ridgeplot(gsea_result)
gseaplot2(gsea_result, geneSetID = "GO:0008150") # Replace with your gene set ID

# Inspect the results
head(gsea_result)
sig_gsea <- subset(gsea_result, p.adjust < 0.05)


###############################################
# Read in data ===================================================
list.files(in_path)
df <- read.csv(paste0(in_path, 'severevshealthy_degresults.csv'), row.names = 1)
# Annotate according to differential expression
df <- df %>% mutate(diffexpressed = case_when(
  log2fc > 0 & padj < 0.05 ~ 'UP',
  log2fc < 0 & padj < 0.05 ~ 'DOWN',
  padj > 0.05 ~ 'NO'
))








##################################
# Annotate according to differential expression
df <- de_genes %>% mutate(diffexpressed = case_when(
  log2FoldChange > 0 & padj < 0.05 ~ 'UP',
  log2FoldChange < 0 & padj < 0.05 ~ 'DOWN',
  padj > 0.05 ~ 'NO'
))

# Split the data into up- and down-regulated genes
upregulated_genes <- df %>% filter(diffexpressed == 'UP')
downregulated_genes <- df %>% filter(diffexpressed == 'DOWN')

# Save the lists if needed
write.csv(upregulated_genes, file = paste0(in_path, 'upregulated_genes.csv'))
write.csv(downregulated_genes, file = paste0(in_path, 'downregulated_genes.csv'))

# Perform GO enrichment analysis on up-regulated genes
up_gene_list <- upregulated_genes$Gene
ego_up <- enrichGO(
  gene = up_gene_list,
  OrgDb = org.Aumbell.eg.db,
  keyType = "GID",
  ont = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = FALSE
)

# Perform GO enrichment analysis on down-regulated genes
down_gene_list <- downregulated_genes$Gene
ego_down <- enrichGO(
  gene = down_gene_list,
  OrgDb = org.Aumbell.eg.db,
  keyType = "GID",
  ont = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = FALSE
)

# Plotting options
# Example for up-regulated genes
dotplot(ego_up, showCategory = 20) + ggtitle("Dotplot for GO Enrichment Analysis of Up-regulated Genes")
barplot(ego_up, showCategory = 20) + ggtitle("Barplot for GO Enrichment Analysis of Up-regulated Genes")

# Example for down-regulated genes
dotplot(ego_down, showCategory = 10) + ggtitle("Dotplot for GO Enrichment Analysis of Down-regulated Genes")

# Combine results into a single data frame
ego_up@result$Regulation <- "UP"
ego_down@result$Regulation <- "DOWN"
combined_results <- bind_rows(ego_up@result, ego_down@result)

# Convert combined results back to enrichResult object for plotting
combined_ego <- new("enrichResult", result = combined_results, pvalueCutoff = 0.05, pAdjustMethod = "BH", ont = "BP", organism = "your_organism", gene = c(upregulated_genes, downregulated_genes), keytype = "GID")

# Plot combined results
dotplot(combined_ego, showCategory = 20, split = "Regulation") + facet_grid(~Regulation) + ggtitle("GO Enrichment Analysis of Up- and Down-regulated Genes")


###########################
# Prepare the geneList for GSEA
# Annotate according to differential expression
df <- de_genes

# Split the data into up- and down-regulated genes
upregulated_genes <- df %>% filter(log2FoldChange > 0 & padj < 0.05) %>% pull(Gene)
downregulated_genes <- df %>% filter(log2FoldChange < 0 & padj < 0.05) %>% pull(Gene)

# Perform GSEA for up-regulated genes
geneList_up <- df[df$Gene %in% upregulated_genes, "log2fc"]
names(geneList_up) <- df[df$Gene %in% upregulated_genes, "Gene"]

gsea_result_up <- gseGO(
  geneList = geneList_up,
  OrgDb = org.Aumbell.eg.db,
  keyType = "GID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  verbose = TRUE
)

# Perform GSEA for down-regulated genes
geneList_down <- df[df$Gene %in% downregulated_genes, "log2fc"]
names(geneList_down) <- df[df$Gene %in% downregulated_genes, "Gene"]

gsea_result_down <- gseGO(
  geneList = geneList_down,
  OrgDb = org.Aumbell.eg.db,
  keyType = "GID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  verbose = TRUE
)

# View GSEA results
head(gsea_result_up)
head(gsea_result_down)

# Plotting GSEA results
# Example for up-regulated genes
dotplot(gsea_result_up, showCategory = 20) + ggtitle("GSEA for Up-regulated Genes")

# Example for down-regulated genes
dotplot(gsea_result_down, showCategory = 20) + ggtitle("GSEA for Down-regulated Genes")

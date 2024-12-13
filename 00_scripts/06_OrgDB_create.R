library(tidyverse)
library(GO.db)
library(dplyr)
library(AnnotationDbi)
library(AnnotationForge)


################################################################################
# creating OrgDb

# import the GO terms and gene ids parsed from interproscan output file
genes_go_annotations <- read.csv("gene_ontologies_all_genes.csv")

#filter unique GO terms
uniq_goTerms <- unique(genes_go_annotations$go)

#get the functional terms and description of the go terms in Ae. umbellulata BRAKER3 annotation
GO_terms_description <- AnnotationDbi::select(GO.db, keys = uniq_goTerms, columns = c("GOID", "TERM", "ONTOLOGY", "DEFINITION"))

# Convert to data frame
GO_terms_description_df <- data.frame(
  go_id = GO_terms_description$GOID,
  term = GO_terms_description$TERM,
  description = GO_terms_description$DEFINITION,
  ontology = GO_terms_description$ONTOLOGY
)

head(GO_terms_description_df)

# merge the df by go terms to include gene names
merged_data <- merge(genes_go_annotations, GO_terms_description_df, by.x = "go", by.y = "go_id")

umb_anno <- merged_data %>%
  dplyr::select(Gene, go, ontology, term, description) %>%
  dplyr::rename(GID = Gene, GENENAME = term, ONTOLOGY = ontology, GO = go, DEFINITION = description)

uSYM <- umb_anno %>%
  dplyr::select(GID) %>%
  distinct() %>%
  mutate(Alias = GID)

uGOFile <- umb_anno %>%
  dplyr::select(GID, GO) %>%
  distinct() %>%
  mutate(EVIDENCE = "IEA")

uTermsFile <- umb_anno %>%
  dplyr::select(GID, GENENAME) %>%
  distinct()  %>%
  na.omit()

uOntologyFile <- umb_anno %>%
  dplyr::select(GID, ONTOLOGY) %>%
  distinct()  %>%
  na.omit()

uDescriptionFile <- umb_anno %>%
  dplyr::select(GID, DEFINITION) %>%
  distinct()  %>%
  na.omit()

# Create OrgDb package
orgDb <- AnnotationForge::makeOrgPackage(
  gene_info = uSYM,                 # Merged data containing gene info and GO annotations
  go = uGOFile,
  #terms = uTermsFile,
  #ontology = uOntologyFile,
  #description = uDescriptionFile,
  version = "1.0",                         # Version number
  author = "Jatinder Singh",                    # Author's name
  maintainer = "Jatinder Singh <jatinder.singh@ndsu.edu>",    # Maintainer's email
  tax_id = "4491",
  genus = "Aegilops",
  species = "umbellulata",
  outputDir = ".",   # Directory to save the OrgDb package
  goTable="go"
)
################################################################################

library(tidyverse)
library(data.table)
library(splitstackshape)
interpro <- read.delim("../../13_interproscan_on_braker3_isoseqRNAseq_again/00_finally/02_all_appl/aumb_brakerRNAref_IsoSeqcomp_functional_renamed_astrix_removed_all.tsv", header = F, stringsAsFactors = F, sep = "\t", row.names = NULL)
head(interpro)
new_colnames <- c("protein_id",
                  "seqMD5",
                  "seqLength",
                  "analysis",
                  "signature_accession",
                  "signature_description",
                  "start",
                  "stop",
                  "score",
                  "status",
                  "date",
                  "interpro_accession",
                  "interpro_description",
                  "go",
                  "pathways")
colnames(interpro) <- new_colnames

interpro_go <- 
  interpro %>% 
  dplyr::select(protein_id, go) %>% 
  dplyr::filter(go != "-") %>% 
  dplyr::filter(go != "") %>%
  mutate(protein_id = gsub("\\.\\d+$", "", protein_id))
tail(interpro_go, n = 10)

##split go terms by | for each gene
##Convert the data frame to a data table
interpro_go_split <- setDT(interpro_go)
interpro_go_split <- interpro_go_split[, .(go = unlist(strsplit(as.character(go), "\\|"))), by = protein_id]
dedup_go <- splitted_interpro_go %>% distinct()
tail(dedup_go)
dedup_go

## another method which also split the Go terms by |
## Use this one if you want (Both do the same thing, this one give some warming)
## splitted_interpro_go <- cSplit(indt = interpro_go, splitCols = "go", sep = "|", direction = "long")
## dedup_go  <- splitted_interpro_go %>% distinct()
## dedup_go

# change column name
dedup_go <- dedup_go %>% 
  dplyr::rename("Gene" = "protein_id")
  ##mutate(Gene = gsub("\\.\\d+$", "", Gene))
tail(dedup_go)
dedup_go
write.csv(x = dedup_go, 
          file = "gene_ontologies_all_genes.csv", 
          row.names = F)
################################################################################
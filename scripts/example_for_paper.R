library(dplyr)
library(readr)
library(tidyr)

setwd("~/Desktop/MTA/CMC_project/rmellea/")

genes <- read_tsv("rmellea_genes.tsv", col_names = T, cols(.default = col_guess(), gene_id = col_character()))
isoforms <- read_tsv("rmellea_isoforms.tsv", col_names = T, cols(.default = col_guess(), gene_id = col_character()))





for(x in 1:length(genes$gene_id)) {
  genes[x,11] <- max(genes[x,2:6]) / min(genes[x,2:6])
}

IDs <- filter(genes, V11 < 1.4) %>% filter(isoforms > 1) %>% select(gene_id)


isoselect <- as_tibble()
for (x in IDs$gene_id) {
  iso_temp <- filter(isoforms, gene_id == x)
  isoselect <- rbind(isoselect, iso_temp)
}

isoselect <- filter(isoselect, `ISOFORM-devreg-s-stricto` == T) 



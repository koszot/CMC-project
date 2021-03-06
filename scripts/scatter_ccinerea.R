library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)

setwd("~/Desktop/MTA/CMC_project/ccinerea/")

genes <- read_tsv("ccinerea_genes.tsv", col_names = T, cols(.default = col_guess(), gene_id = col_character()))
isoforms <- read_tsv("ccinerea_isoforms.tsv", col_names = T, cols(.default = col_guess(), gene_id = col_character()))

genes_AS <- genes %>% filter(isoforms > 1)
isoforms_AS <- isoforms %>% filter(`GENE-isoforms` > 1)

for(x in 1:length(genes_AS$gene_id)) {
  genes_AS$gene_FC[x] <- max(genes_AS[x,2:10]) / min(genes_AS[x,2:10])
}

genes_AS_FC <- genes_AS %>% select(gene_id, gene_FC)

for(x in 1:length(isoforms_AS$gene_id)) {
  isoforms_AS$isoform_FC[x] <- max(isoforms_AS[x,3:11]) / min(isoforms_AS[x,3:11])
  if (isoforms_AS$`GENE-devreg-s-stricto`[x] == T & isoforms_AS$`ISOFORM-devreg-s-stricto`[x] == T) {
    isoforms_AS$devreg_factor[x] <- "Gene Devreg / Isoform Devreg"
  } else if (isoforms_AS$`GENE-devreg-s-stricto`[x] == T & isoforms_AS$`ISOFORM-devreg-s-stricto`[x] == F) {
    isoforms_AS$devreg_factor[x] <- "Gene Devreg / Isoform NOT Devreg"
  } else if (isoforms_AS$`GENE-devreg-s-stricto`[x] == F & isoforms_AS$`ISOFORM-devreg-s-stricto`[x] == F) {
    isoforms_AS$devreg_factor[x] <- "Gene NOT Devreg / Isoform NOT Devreg"
  } else if (isoforms_AS$`GENE-devreg-s-stricto`[x] == F & isoforms_AS$`ISOFORM-devreg-s-stricto`[x] == T) {
    isoforms_AS$devreg_factor[x] <- "Gene NOT Devreg / Isoform Devreg"
  }
}

isoforms_AS_FC <- isoforms_AS %>% select(gene_id, transcript_id, isoform_FC, devreg_factor)

full_FC <- left_join(isoforms_AS_FC, genes_AS_FC)

full_FC_2 <- full_FC %>% filter(gene_FC < 1000) %>% filter(isoform_FC < 1000)

full_FC_3 <- full_FC %>% filter(gene_FC < 10) %>% filter(isoform_FC < 10)

e <- ggplot(full_FC, aes(isoform_FC, gene_FC))
e + geom_point(aes(colour = factor(devreg_factor))) +
  theme(legend.position="right" ) +
  labs(x="Isoform's Fold Change", y="Gene's Fold Change")

e <- ggplot(full_FC_2, aes(isoform_FC, gene_FC))
e + geom_point(aes(colour = factor(devreg_factor))) +
  theme(legend.position="right" ) +
  coord_fixed() + guides(col = guide_legend(nrow = 4)) +
  labs(x="Isoform's Fold Change", y="Gene's Fold Change")

e <- ggplot(full_FC_3, aes(isoform_FC, gene_FC))
e + geom_point(aes(colour = factor(devreg_factor))) +
  theme(legend.position="right" ) +
  coord_fixed() + guides(col = guide_legend(nrow = 4)) +
  labs(x="Isoform's Fold Change", y="Gene's Fold Change")

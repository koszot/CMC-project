library(dplyr)
library(readr)
library(tidyr)

setwd("~/Desktop/MTA/CMC_project/ccinerea/")

genes <- read_tsv("ccinerea_genes.tsv", col_names = T, cols(.default = col_guess(), gene_id = col_character()))
devreg <- genes %>% filter(`FB-devreg` == T) %>% select(gene_id)
init <- genes %>% filter(`FB-init` == T) %>% select(gene_id)
as <- genes %>% filter(isoforms > 1) %>% select(gene_id)
not_devreg <- genes %>% filter(`FB-devreg` == F) %>% dplyr::select(gene_id)
not_init <- genes %>% filter(`FB-init` == F) %>% dplyr::select(gene_id)
isoforms <- read_tsv("ccinerea_isoforms.tsv", col_names = T, cols(.default = col_guess(), gene_id = col_character()))

library(ASpli)

##### AS GENES #####

TxDb <- makeTxDbFromGFF(file="ccinerea_AS_genes.gtf", format="gtf")
features <- binGenome(TxDb)
AS_count <- read_lines("ASpli_binFeatures.log")
AS_count <- as_tibble(AS_count[-c(1,2,3,4,5,6,7,8)]) %>%
  separate(value, c("l1", "l2", "l3", "l4", "l5"), sep = "\t") %>%
  dplyr::select(l2, l4) %>%
  separate(l2, c("event", "count"), sep = "=") %>%
  separate(l4, c("event2", "count2"), sep = "=") 
AS_count$count <- as.integer(AS_count$count)
AS_count$count2 <- as.integer(AS_count$count2)
AS_count[1,2] <- AS_count[1,2] + AS_count[7,4]
AS_count[2,2] <- AS_count[2,2] + AS_count[8,4]
AS_count[3,2] <- AS_count[3,2] + AS_count[9,4]
AS_count[4,2] <- AS_count[4,2] + AS_count[10,4]
AS_count <- AS_count[-(5:10),-(3:4)]
AS_AS_count <- AS_count

##### FB DEVREG GENES #####

TxDb <- makeTxDbFromGFF(file="ccinerea_FB_DEVREG_genes.gtf", format="gtf")
features <- binGenome(TxDb)
AS_count <- read_lines("ASpli_binFeatures.log")
AS_count <- as_tibble(AS_count[-c(1,2,3,4,5,6,7,8)]) %>%
  separate(value, c("l1", "l2", "l3", "l4", "l5"), sep = "\t") %>%
  dplyr::select(l2, l4) %>%
  separate(l2, c("event", "count"), sep = "=") %>%
  separate(l4, c("event2", "count2"), sep = "=") 
AS_count$count <- as.integer(AS_count$count)
AS_count$count2 <- as.integer(AS_count$count2)
AS_count[1,2] <- AS_count[1,2] + AS_count[7,4]
AS_count[2,2] <- AS_count[2,2] + AS_count[8,4]
AS_count[3,2] <- AS_count[3,2] + AS_count[9,4]
AS_count[4,2] <- AS_count[4,2] + AS_count[10,4]
AS_count <- AS_count[-(5:10),-(3:4)]
FB_devreg_AS_count <- AS_count

##### FB INIT GENES #####

TxDb <- makeTxDbFromGFF(file="ccinerea_FB_INIT_genes.gtf", format="gtf")
features <- binGenome(TxDb)
AS_count <- read_lines("ASpli_binFeatures.log")
AS_count <- as_tibble(AS_count[-c(1,2,3,4,5,6,7,8)]) %>%
  separate(value, c("l1", "l2", "l3", "l4", "l5"), sep = "\t") %>%
  dplyr::select(l2, l4) %>%
  separate(l2, c("event", "count"), sep = "=") %>%
  separate(l4, c("event2", "count2"), sep = "=") 
AS_count$count <- as.integer(AS_count$count)
AS_count$count2 <- as.integer(AS_count$count2)
AS_count[1,2] <- AS_count[1,2] + AS_count[7,4]
AS_count[2,2] <- AS_count[2,2] + AS_count[8,4]
AS_count[3,2] <- AS_count[3,2] + AS_count[9,4]
AS_count[4,2] <- AS_count[4,2] + AS_count[10,4]
AS_count <- AS_count[-(5:10),-(3:4)]
FB_init_AS_count <- AS_count

##### DEVREG GENES  /w or w/o DEVREG TRANSCRIPTS #####

FB_devreg_genes <- left_join(devreg, isoforms) %>% filter(`GENE-isoforms` > 1)
FB_devreg_genes_iso_devreg <- 0
FB_devreg_genes_iso_not_devreg <- 0
for (x in 1:length(unique(FB_devreg_genes$gene_id))) {
  temp <- FB_devreg_genes[FB_devreg_genes$gene_id %in% unique(FB_devreg_genes$gene_id)[1],]
  if (sum(temp$`ISOFORM-FB-devreg`) > 0) {
    FB_devreg_genes_iso_devreg <- FB_devreg_genes_iso_devreg + 1
  } else {
    FB_devreg_genes_iso_not_devreg <- FB_devreg_genes_iso_not_devreg + 1
  }
}

FB_init_genes <- left_join(init, isoforms) %>% filter(`GENE-isoforms` > 1)
FB_init_genes_iso_devreg <- 0
FB_init_genes_iso_not_devreg <- 0
for (x in 1:length(unique(FB_init_genes$gene_id))) {
  temp <- FB_init_genes[FB_init_genes$gene_id %in% unique(FB_init_genes$gene_id)[x],]
  if (sum(temp$`ISOFORM-FB-init`) > 0) {
    FB_init_genes_iso_devreg <- FB_init_genes_iso_devreg + 1
  } else {
    FB_init_genes_iso_not_devreg <- FB_init_genes_iso_not_devreg + 1
  }
}

##### NOT DEVREG GENES  /w or w/o DEVREG TRANSCRIPTS #####

not_FB_devreg_genes <- left_join(not_devreg, isoforms) %>% filter(`GENE-isoforms` > 1)
not_FB_devreg_genes_iso_devreg <- 0
not_FB_devreg_genes_iso_not_devreg <- 0
for (x in 1:length(unique(not_FB_devreg_genes$gene_id))) {
  temp <- not_FB_devreg_genes[not_FB_devreg_genes$gene_id %in% unique(not_FB_devreg_genes$gene_id)[x],]
  if (sum(temp$`ISOFORM-FB-devreg`) > 0) {
    not_FB_devreg_genes_iso_devreg <- not_FB_devreg_genes_iso_devreg + 1
  } else {
    not_FB_devreg_genes_iso_not_devreg <- not_FB_devreg_genes_iso_not_devreg + 1
  }
}

not_FB_init_genes <- left_join(not_init, isoforms) %>% filter(`GENE-isoforms` > 1)
not_FB_init_genes_iso_devreg <- 0
not_FB_init_genes_iso_not_devreg <- 0
for (x in 1:length(unique(not_FB_init_genes$gene_id))) {
  temp <- not_FB_init_genes[not_FB_init_genes$gene_id %in% unique(not_FB_init_genes$gene_id)[x],]
  if (sum(temp$`ISOFORM-FB-init`) > 0) {
    not_FB_init_genes_iso_devreg <- not_FB_init_genes_iso_devreg + 1
  } else {
    not_FB_init_genes_iso_not_devreg <- not_FB_init_genes_iso_not_devreg + 1
  }
}

##### STATS #####

stats <- as_tibble()

stats[1,1] <- "Overall"
stats[1,2] <- nrow(genes)
stats[2,1] <- "FB-devreg"
stats[2,2] <- nrow(devreg)
stats[3,1] <- "FB-init"
stats[3,2] <- nrow(init)
stats[4,1] <- "Alternatively Spliced"
stats[4,2] <- nrow(as)
stats[5,1] <- "Alternatively Spliced --- ES bins"
stats[5,2] <- AS_AS_count$count[1]
stats[6,1] <- "Alternatively Spliced --- IR bins"
stats[6,2] <- AS_AS_count$count[2]
stats[7,1] <- "Alternatively Spliced --- Alt 5' SS bins"
stats[7,2] <- AS_AS_count$count[3]
stats[8,1] <- "Alternatively Spliced --- Alt 3' SS bins"
stats[8,2] <- AS_AS_count$count[4]
stats[9,1] <- "Alternatively Spliced / FB-devreg"
stats[9,2] <- nrow(dplyr::intersect(devreg, as))
stats[10,1] <- "Alternatively Spliced / FB-devreg -- ES bins"
stats[10,2] <- FB_devreg_AS_count$count[1]
stats[11,1] <- "Alternatively Spliced / FB-devreg -- IR bins"
stats[11,2] <- FB_devreg_AS_count$count[2]
stats[12,1] <- "Alternatively Spliced / FB-devreg -- Alt 5' SS bins"
stats[12,2] <- FB_devreg_AS_count$count[3]
stats[13,1] <- "Alternatively Spliced / FB-devreg -- Alt 3' SS bins"
stats[13,2] <- FB_devreg_AS_count$count[4]
stats[14,1] <- "Alternatively Spliced / FB-devreg / FB-devreg Isoforms"
stats[14,2] <- FB_devreg_genes_iso_devreg
stats[15,1] <- "Alternatively Spliced / FB-devreg / NO FB-devreg Isoforms"
stats[15,2] <- FB_devreg_genes_iso_not_devreg
stats[16,1] <- "Alternatively Spliced / NOT FB-devreg"
stats[16,2] <- length(unique(not_FB_devreg_genes$gene_id))
stats[17,1] <- "Alternatively Spliced / NOT FB-devreg / FB-devreg Isoforms"
stats[17,2] <- not_FB_devreg_genes_iso_devreg
stats[18,1] <- "Alternatively Spliced / NOT FB-devreg / NO FB-devreg Isoforms"
stats[18,2] <- not_FB_devreg_genes_iso_not_devreg
stats[19,1] <- "Alternatively Spliced / FB-init"
stats[19,2] <- nrow(dplyr::intersect(init, as))
stats[20,1] <- "Alternatively Spliced / FB-init --- ES bins"
stats[20,2] <- FB_init_AS_count$count[1]
stats[21,1] <- "Alternatively Spliced / FB-init --- IR bins"
stats[21,2] <- FB_init_AS_count$count[2]
stats[22,1] <- "Alternatively Spliced / FB-init --- Alt 5' SS bins"
stats[22,2] <- FB_init_AS_count$count[3]
stats[23,1] <- "Alternatively Spliced / FB-init --- Alt 3' SS bins"
stats[23,2] <- FB_init_AS_count$count[4]
stats[24,1] <- "Alternatively Spliced / FB-init / FB-init Isoforms"
stats[24,2] <- FB_init_genes_iso_devreg
stats[25,1] <- "Alternatively Spliced / FB-init / NO FB-init Isoforms"
stats[25,2] <- FB_init_genes_iso_not_devreg
stats[26,1] <- "Alternatively Spliced / NOT FB-init"
stats[26,2] <- length(unique(not_FB_init_genes$gene_id))
stats[27,1] <- "Alternatively Spliced / NOT FB-init / FB-init Isoforms"
stats[27,2] <- not_FB_init_genes_iso_devreg
stats[28,1] <- "Alternatively Spliced / NOT FB-init / NO FB-init Isoforms"
stats[28,2] <- not_FB_init_genes_iso_not_devreg

colnames(stats) <- c("Description", "Coprinopsis cinerea")
write_tsv(stats, "ccinerea_stats.tsv")

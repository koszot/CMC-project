library(tidyr)
library(readr)
library(dplyr)
library(stringr)

##### LOAD FILES #####

setwd("~/Desktop/MTA_CMC_project/ltigrinus/")

genes <- read_tsv("ltigrinus_genes.tsv", col_names = T, cols(.default = col_guess(), gene_id = col_character()))
isoforms <- read_tsv("ltigrinus_isoforms.tsv", col_names = T, cols(.default = col_guess(), gene_id = col_character()))
stats <- as_tibble()

annotation <- read_tsv("ltigrinus_corrected_annotation.gtf",  col_names = c("chr", "maker","type", "start", "end", "att1", "strand", "att2", "attributes")) %>%
  separate(attributes, c("transcriptID_label", "transcriptID", "geneID_label", "geneID"), sep = " ")
annotation$geneID <- annotation$geneID %>% 
  str_replace("\"", "") %>%
  str_replace("\";", "")              # transcriptID-k átalakítása

##### BASIC STATISTICS #####

stats[1,1] <- "genes"
stats[1,2] <- nrow(genes)
stats[2,1] <- "FB-devreg genes"
stats[2,2] <- nrow(genes %>% filter(`FB-devreg` == T))
stats[3,1] <- "FB-init genes"
stats[3,2] <- nrow(genes %>% filter(`FB-init` == T))
stats[4,1] <- "genes with AS"
stats[4,2] <- nrow(genes %>% filter(isoforms > 1))
stats[5,1] <- "FB-devreg genes with AS"
stats[5,2] <- nrow(genes %>% filter(isoforms > 1) %>% filter(`FB-devreg` == T))
stats[6,1] <- "FB-init genes with AS"
stats[6,2] <- nrow(genes %>% filter(isoforms > 1) %>% filter(`FB-init` == T))

colnames(stats) <- c("names", "ccinerea")
write_tsv(stats, "ltigrinus_stats.tsv")

##### AS ANNOTATION #####

as <- genes %>% filter(isoforms > 1)
annotation.as <- as_tibble()

for (x in 1:length(as$gene_id)) {
  temp <- annotation[annotation$geneID %in% as$gene_id[x],]
  annotation.as <- rbind(annotation.as, temp)
  pb <- txtProgressBar(min = 1, max = length(as$gene_id), style = 3)        
  setTxtProgressBar(pb, x, title = NULL, label = NULL) 
}

annotation.as$geneID <- annotation.as$geneID %>% 
  str_replace("^", "\"") %>%
  str_replace("$", "\";")              # transcriptID-k átalakítása
merged <- unite(annotation.as, attributes, 9:12, sep = " ")
write.table(merged, file = "ltigrinus_AS_genes.gtf", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

##### FB-DEVREG ANNOTATION #####

fb.devreg <- genes %>% filter(isoforms > 1) %>% filter(`FB-devreg` == T)
annotation.fb.devreg <- as_tibble()

for (x in 1:length(fb.devreg$gene_id)) {
  temp <- annotation[annotation$geneID %in% fb.devreg$gene_id[x],]
  annotation.fb.devreg <- rbind(annotation.fb.devreg, temp)
  pb <- txtProgressBar(min = 1, max = length(fb.devreg$gene_id), style = 3)        
  setTxtProgressBar(pb, x, title = NULL, label = NULL) 
}

annotation.fb.devreg$geneID <- annotation.fb.devreg$geneID %>% 
  str_replace("^", "\"") %>%
  str_replace("$", "\";")              # transcriptID-k átalakítása
merged <- unite(annotation.fb.devreg, attributes, 9:12, sep = " ")
write.table(merged, file = "ltigrinus_FB_DEVREG_genes.gtf", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

##### FB-INIT ANNOTATION #####

fb.init <- genes %>% filter(isoforms > 1) %>% filter(`FB-init` == T)
annotation.fb.init <- as_tibble()

for (x in 1:length(fb.init$gene_id)) {
  temp <- annotation[annotation$geneID %in% fb.init$gene_id[x],]
  annotation.fb.init <- rbind(annotation.fb.init, temp)
  pb <- txtProgressBar(min = 1, max = length(fb.init$gene_id), style = 3)        
  setTxtProgressBar(pb, x, title = NULL, label = NULL) 
}

annotation.fb.init$geneID <- annotation.fb.init$geneID %>% 
  str_replace("^", "\"") %>%
  str_replace("$", "\";")              # transcriptID-k átalakítása
merged <- unite(annotation.fb.init, attributes, 9:12, sep = " ")
write.table(merged, file = "ltigrinus_FB_INIT_genes.gtf", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

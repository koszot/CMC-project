library(tidyr)
library(readr)
library(dplyr)
library(stringr)

##### LOAD FILES #####

setwd("~/Desktop/MTA_CMC_project/aostoyae/")

genes <- read_tsv("aostoyae_genes.tsv")
isoforms <- read_tsv("aostoyae_isoforms.tsv")
stats <- as_tibble()

annotation <- read_tsv("aostoyae_corrected_annotation.gtf",  col_names = c("chr", "maker","type", "start", "end", "att1", "strand", "att2", "attributes")) %>%
  separate(attributes, c("transcriptID_label", "transcriptID", "geneID_label", "geneID"), sep = " ")
annotation$geneID <- annotation$geneID %>% 
  str_replace("\"", "") %>%
  str_replace("\";", "")              # transcriptID-k átalakítása

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
write.table(merged, file = "aostoyae_AS_genes.gtf", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

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
write.table(merged, file = "aostoyae_FB_DEVREG_genes.gtf", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

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
write.table(merged, file = "aostoyae_FB_INIT_genes.gtf", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

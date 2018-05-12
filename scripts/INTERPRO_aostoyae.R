library(dplyr)
library(readr)
library(tidyr)
library(stringr)

setwd("~/Desktop/MTA/CMC_project/aostoyae/")
isoforms <- read_tsv("aostoyae_isoforms.tsv", col_names = T, cols(.default = col_guess(), gene_id = col_character()))
raw_interpro <- read_lines("p3_i2_t47428_Arm_ostoy_v2.prot.tsv")
IPR <- as_tibble(raw_interpro[str_detect(raw_interpro, "\tIPR")], header = F)
colnames(IPR) <- "all"
IPR <- IPR %>%
  separate(all, c("transcript_id", "hash", "length", "type", 
                  "domain_id", "domin_desc", "start", "end", 
                  "E-value", "attribute", "date", "IPR_domain", 
                  "IPR_domain_name", "GO"), sep = "\t")



# Alternatively Spliced / NOT FB-devreg / FB-devreg Isoforms
NOdevreg_as_devreg <- isoforms %>% filter(`GENE-FB-devreg` == F) %>% filter(`GENE-isoforms` > 1) %>% filter(`ISOFORM-FB-devreg` == T) %>% select(gene_id) %>% unique()
# Alternatively Spliced / NOT FB-init / FB-init Isoforms
NOinit_as_init <- isoforms %>% filter(`GENE-FB-init` == F) %>% filter(`GENE-isoforms` > 1) %>% filter(`ISOFORM-FB-init` == T) %>% select(gene_id) %>% unique()

overall <- union(NOdevreg_as_devreg, NOinit_as_init)
intersect(NOdevreg_as_devreg, NOinit_as_init)
filtered_IPR <- as_tibble()

for (x in 1:length(overall$gene_id)) {
  temp <- IPR[IPR$transcript_id %in% overall$gene_id[x],]
  filtered_IPR <- rbind(filtered_IPR, temp)
  pb <- txtProgressBar(min = 1, max = length(overall$gene_id), style = 3)        
  setTxtProgressBar(pb, x, title = NULL, label = NULL) 
}

filtered_IPR <- filtered_IPR[,c(1,12,13,14)]

write_tsv(filtered_IPR, "aostoyae_filtered_IPR.tsv")

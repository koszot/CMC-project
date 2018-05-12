library(dplyr)
library(readr)
library(tidyr)
library(stringr)

##### BEOLVASSUK A SZÜKSÉGES FÁJLOKAT #####

setwd("~/Desktop/MTA/CMC_project/ccinerea/")
isoforms <- read_tsv("ccinerea_isoforms.tsv", col_names = T, cols(.default = col_guess(), gene_id = col_character()))
raw_interpro <- read_lines("Copci_AmutBmut1_GeneCatalog_proteins_20130522.aa.fasta.tsv")
IPR <- as_tibble(raw_interpro[str_detect(raw_interpro, "\tIPR")], header = F)
colnames(IPR) <- "all"
IPR <- IPR %>%
  separate(all, c("transcript_id", "hash", "length", "type", 
                  "domain_id", "domin_desc", "start", "end", 
                  "E-value", "attribute", "date", "IPR_domain", 
                  "IPR_domain_name", "GO"), sep = "\t")

##### MEGCSINÁLJUK A SZÓTÁRFÁJLT AZ isoforms SZÁMÁRA #####

# beolvassuk a headerfájlból készített szótárfájlt
fasta <- read_tsv("../../AS_project/FILES_ccinerea/ENRICHMENT_ccinerea/Copci_AmutBmut1_GeneCatalog_proteins_20130522.aa.fasta.headers", col_names = c("fasta_transcript_id", "gene_id"))
# beolvassuk az eredeti annotációt ami alapján az expressziós analyzis készült
original <- read_tsv("../../AS_project/FILES_ccinerea/GENOME_ccinerea/Copci_AmutBmut1_GeneCatalog_genes_20130522.gff", col_names = c("chr", "maker", "type", "start", "end", "att1", "strand",  "att2", "attributes"))
original <- original %>%
  filter(type == "exon") %>%
  separate(attributes, c("name", "gene_id", "tid", "annotation_transcript_id"), sep = " ") %>%
  select(gene_id, annotation_transcript_id)
original$gene_id <- original$gene_id %>% 
  str_replace("^\"", "") %>%
  str_replace("\";$", "")
original <- original[!duplicated(original), ]
final <- left_join(fasta, original, by = "gene_id")
names(final) <- c("CMC_gene_id", "gene_name", "gene_id")

##### KIKERESSÜK AZ ID-KET #####

# Alternatively Spliced / NOT FB-devreg / FB-devreg Isoforms
NOdevreg_as_devreg <- isoforms %>% filter(`GENE-FB-devreg` == F) %>% filter(`GENE-isoforms` > 1) %>% filter(`ISOFORM-FB-devreg` == T) %>% select(gene_id) %>% unique()
# Alternatively Spliced / NOT FB-init / FB-init Isoforms
NOinit_as_init <- isoforms %>% filter(`GENE-FB-init` == F) %>% filter(`GENE-isoforms` > 1) %>% filter(`ISOFORM-FB-init` == T) %>% select(gene_id) %>% unique()

overall <- union(NOdevreg_as_devreg, NOinit_as_init)
# megoldjuk az transcript/protein ID gondot
overall <- left_join(overall, final, by = "gene_id") %>% select(CMC_gene_id)

intersect(NOdevreg_as_devreg, NOinit_as_init)
filtered_IPR <- as_tibble()

for (x in 1:length(overall$CMC_gene_id)) {
  temp <- IPR[IPR$transcript_id %in% overall$CMC_gene_id[x],]
  filtered_IPR <- rbind(filtered_IPR, temp)
  pb <- txtProgressBar(min = 1, max = length(overall$CMC_gene_id), style = 3)        
  setTxtProgressBar(pb, x, title = NULL, label = NULL) 
}

length(unique(filtered_IPR$transcript_id))
filtered_IPR <- filtered_IPR[,c(1,12,13,14)]

write_tsv(filtered_IPR, "ccinerea_filtered_IPR.tsv")

##### HIBATESZTELÉS (ha itt nincs meg az összes átalakított isoforms ID az IPR összesben akkor baj van) #####

test1 <- unique(IPR$transcript_id)
test2 <- isoforms %>% select(gene_id) %>% unique()
test2 <- left_join(test2, final, by = "gene_id") %>% select(CMC_gene_id)
test2 <- test2$CMC_gene_id

length(test1)
length(intersect(test1, test2))





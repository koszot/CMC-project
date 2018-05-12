library(dplyr)
library(readr)
library(tidyr)

setwd("~/Desktop/MTA/CMC_project/ccinerea/")
isoforms <- read_tsv("ccinerea_isoforms.tsv", col_names = T, cols(.default = col_guess(), gene_id = col_character()))

##### DICTIONARY #####

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

##### STATS #####

# FB-devreg
devreg <- isoforms %>% filter(`GENE-FB-devreg` == T) %>% select(gene_id) %>% unique()
devreg <- left_join(devreg, final, by = "gene_id") %>% select(CMC_gene_id)
devreg$CMC_gene_id <- str_replace(devreg$CMC_gene_id, "^", "Copci_AmutBmut1_")
write_tsv(devreg, "../temp_tables/devreg.tsv")
# FB-init
init <- isoforms %>% filter(`GENE-FB-init` == T) %>% select(gene_id) %>% unique()
init <- left_join(init, final, by = "gene_id") %>% select(CMC_gene_id)
init$CMC_gene_id <- str_replace(init$CMC_gene_id, "^", "Copci_AmutBmut1_")
write_tsv(init, "../temp_tables/init.tsv")
# Alternatively Spliced
as <- isoforms %>% filter(`GENE-isoforms` > 1) %>% select(gene_id) %>% unique()
as <- left_join(as, final, by = "gene_id") %>% select(CMC_gene_id)
as$CMC_gene_id <- str_replace(as$CMC_gene_id, "^", "Copci_AmutBmut1_")
write_tsv(as, "../temp_tables/as.tsv")
# Alternatively Spliced / FB-devreg
devreg_as <- isoforms %>% filter(`GENE-FB-devreg` == T) %>% filter(`GENE-isoforms` > 1) %>% select(gene_id) %>% unique()
devreg_as <- left_join(devreg_as, final, by = "gene_id") %>% select(CMC_gene_id)
devreg_as$CMC_gene_id <- str_replace(devreg_as$CMC_gene_id, "^", "Copci_AmutBmut1_")
write_tsv(devreg_as, "../temp_tables/devreg_as.tsv")
# Alternatively Spliced / FB-devreg / FB-devreg Isoforms
devreg_as_devreg <- isoforms %>% filter(`GENE-FB-devreg` == T) %>% filter(`GENE-isoforms` > 1) %>% filter(`ISOFORM-FB-devreg` == T) %>% select(gene_id) %>% unique()
devreg_as_devreg <- left_join(devreg_as_devreg, final, by = "gene_id") %>% select(CMC_gene_id)
devreg_as_devreg$CMC_gene_id <- str_replace(devreg_as_devreg$CMC_gene_id, "^", "Copci_AmutBmut1_")
write_tsv(devreg_as_devreg, "../temp_tables/devreg_as_devreg.tsv")
# Alternatively Spliced / FB-devreg / NO FB-devreg Isoforms
devreg_as_NOdevreg <- setdiff(devreg_as, devreg_as_devreg)
write_tsv(devreg_as_NOdevreg, "../temp_tables/devreg_as_NOdevreg.tsv")
# Alternatively Spliced / NOT FB-devreg
NOdevreg_as <- isoforms %>% filter(`GENE-FB-devreg` == F) %>% filter(`GENE-isoforms` > 1) %>% select(gene_id) %>% unique()
NOdevreg_as <- left_join(NOdevreg_as, final, by = "gene_id") %>% select(CMC_gene_id)
NOdevreg_as$CMC_gene_id <- str_replace(NOdevreg_as$CMC_gene_id, "^", "Copci_AmutBmut1_")
write_tsv(NOdevreg_as, "../temp_tables/NOdevreg_as.tsv")
# Alternatively Spliced / NOT FB-devreg / FB-devreg Isoforms
NOdevreg_as_devreg <- isoforms %>% filter(`GENE-FB-devreg` == F) %>% filter(`GENE-isoforms` > 1) %>% filter(`ISOFORM-FB-devreg` == T) %>% select(gene_id) %>% unique()
NOdevreg_as_devreg <- left_join(NOdevreg_as_devreg, final, by = "gene_id") %>% select(CMC_gene_id)
NOdevreg_as_devreg$CMC_gene_id <- str_replace(NOdevreg_as_devreg$CMC_gene_id, "^", "Copci_AmutBmut1_")
write_tsv(NOdevreg_as_devreg, "../temp_tables/NOdevreg_as_devreg.tsv")
# Alternatively Spliced / NOT FB-devreg / NO FB-devreg Isoforms
NOdevreg_as_NOdevreg <- setdiff(NOdevreg_as, NOdevreg_as_devreg)
write_tsv(NOdevreg_as_NOdevreg, "../temp_tables/NOdevreg_as_NOdevreg.tsv")
# Alternatively Spliced / FB-init
init_as <- isoforms %>% filter(`GENE-FB-init` == T) %>% filter(`GENE-isoforms` > 1) %>% select(gene_id) %>% unique()
init_as <- left_join(init_as, final, by = "gene_id") %>% select(CMC_gene_id)
init_as$CMC_gene_id <- str_replace(init_as$CMC_gene_id, "^", "Copci_AmutBmut1_")
write_tsv(init_as, "../temp_tables/init_as.tsv")
# Alternatively Spliced / FB-init / FB-init Isoforms
init_as_init <- isoforms %>% filter(`GENE-FB-init` == T) %>% filter(`GENE-isoforms` > 1) %>% filter(`ISOFORM-FB-init` == T) %>% select(gene_id) %>% unique()
init_as_init <- left_join(init_as_init, final, by = "gene_id") %>% select(CMC_gene_id)
init_as_init$CMC_gene_id <- str_replace(init_as_init$CMC_gene_id, "^", "Copci_AmutBmut1_")
write_tsv(init_as_init, "../temp_tables/init_as_init.tsv")
# Alternatively Spliced / FB-init / NO FB-init Isoforms
init_as_NOinit <- setdiff(init_as, init_as_init)
write_tsv(init_as_NOinit, "../temp_tables/init_as_NOinit.tsv")
# Alternatively Spliced / NOT FB-init
NOinit_as <- isoforms %>% filter(`GENE-FB-init` == F) %>% filter(`GENE-isoforms` > 1) %>% select(gene_id) %>% unique()
NOinit_as <- left_join(NOinit_as, final, by = "gene_id") %>% select(CMC_gene_id)
NOinit_as$CMC_gene_id <- str_replace(NOinit_as$CMC_gene_id, "^", "Copci_AmutBmut1_")
write_tsv(NOinit_as, "../temp_tables/NOinit_as.tsv")
# Alternatively Spliced / NOT FB-init / FB-init Isoforms
NOinit_as_init <- isoforms %>% filter(`GENE-FB-init` == F) %>% filter(`GENE-isoforms` > 1) %>% filter(`ISOFORM-FB-init` == T) %>% select(gene_id) %>% unique()
NOinit_as_init <- left_join(NOinit_as_init, final, by = "gene_id") %>% select(CMC_gene_id)
NOinit_as_init$CMC_gene_id <- str_replace(NOinit_as_init$CMC_gene_id, "^", "Copci_AmutBmut1_")
write_tsv(NOinit_as_init, "../temp_tables/NOinit_as_init.tsv")
# Alternatively Spliced / NOT FB-init / NO FB-init Isoforms
NOinit_as_NOinit <- setdiff(NOinit_as, NOinit_as_init)
write_tsv(NOinit_as_NOinit, "../temp_tables/NOinit_as_NOinit.tsv")

library(dplyr)
library(readr)
library(tidyr)

setwd("~/Desktop/MTA/CMC_project/aostoyae/")
isoforms <- read_tsv("aostoyae_isoforms.tsv", col_names = T, cols(.default = col_guess(), gene_id = col_character()))

##### STATS #####

# FB-devreg
devreg <- isoforms %>% filter(`GENE-FB-devreg` == T) %>% select(gene_id) %>% unique()
write_tsv(devreg, "../temp_tables/devreg.tsv")
# FB-init
init <- isoforms %>% filter(`GENE-FB-init` == T) %>% select(gene_id) %>% unique()
write_tsv(init, "../temp_tables/init.tsv")
# Alternatively Spliced
as <- isoforms %>% filter(`GENE-isoforms` > 1) %>% select(gene_id) %>% unique()
write_tsv(as, "../temp_tables/as.tsv")
# Alternatively Spliced / FB-devreg
devreg_as <- isoforms %>% filter(`GENE-FB-devreg` == T) %>% filter(`GENE-isoforms` > 1) %>% select(gene_id) %>% unique()
write_tsv(devreg_as, "../temp_tables/devreg_as.tsv")
# Alternatively Spliced / FB-devreg / FB-devreg Isoforms
devreg_as_devreg <- isoforms %>% filter(`GENE-FB-devreg` == T) %>% filter(`GENE-isoforms` > 1) %>% filter(`ISOFORM-FB-devreg` == T) %>% select(gene_id) %>% unique()
write_tsv(devreg_as_devreg, "../temp_tables/devreg_as_devreg.tsv")
# Alternatively Spliced / FB-devreg / NO FB-devreg Isoforms
devreg_as_NOdevreg <- setdiff(devreg_as, devreg_as_devreg)
write_tsv(devreg_as_NOdevreg, "../temp_tables/devreg_as_NOdevreg.tsv")
# Alternatively Spliced / NOT FB-devreg
NOdevreg_as <- isoforms %>% filter(`GENE-FB-devreg` == F) %>% filter(`GENE-isoforms` > 1) %>% select(gene_id) %>% unique()
write_tsv(NOdevreg_as, "../temp_tables/NOdevreg_as.tsv")
# Alternatively Spliced / NOT FB-devreg / FB-devreg Isoforms
NOdevreg_as_devreg <- isoforms %>% filter(`GENE-FB-devreg` == F) %>% filter(`GENE-isoforms` > 1) %>% filter(`ISOFORM-FB-devreg` == T) %>% select(gene_id) %>% unique()
write_tsv(NOdevreg_as_devreg, "../temp_tables/NOdevreg_as_devreg.tsv")
# Alternatively Spliced / NOT FB-devreg / NO FB-devreg Isoforms
NOdevreg_as_NOdevreg <- setdiff(NOdevreg_as, NOdevreg_as_devreg)
write_tsv(NOdevreg_as_NOdevreg, "../temp_tables/NOdevreg_as_NOdevreg.tsv")
# Alternatively Spliced / FB-init
init_as <- isoforms %>% filter(`GENE-FB-init` == T) %>% filter(`GENE-isoforms` > 1) %>% select(gene_id) %>% unique()
write_tsv(init_as, "../temp_tables/init_as.tsv")
# Alternatively Spliced / FB-init / FB-init Isoforms
init_as_init <- isoforms %>% filter(`GENE-FB-init` == T) %>% filter(`GENE-isoforms` > 1) %>% filter(`ISOFORM-FB-init` == T) %>% select(gene_id) %>% unique()
write_tsv(init_as_init, "../temp_tables/init_as_init.tsv")
# Alternatively Spliced / FB-init / NO FB-init Isoforms
init_as_NOinit <- setdiff(init_as, init_as_init)
write_tsv(init_as_NOinit, "../temp_tables/init_as_NOinit.tsv")
# Alternatively Spliced / NOT FB-init
NOinit_as <- isoforms %>% filter(`GENE-FB-init` == F) %>% filter(`GENE-isoforms` > 1) %>% select(gene_id) %>% unique()
write_tsv(NOinit_as, "../temp_tables/NOinit_as.tsv")
# Alternatively Spliced / NOT FB-init / FB-init Isoforms
NOinit_as_init <- isoforms %>% filter(`GENE-FB-init` == F) %>% filter(`GENE-isoforms` > 1) %>% filter(`ISOFORM-FB-init` == T) %>% select(gene_id) %>% unique()
write_tsv(NOinit_as_init, "../temp_tables/NOinit_as_init.tsv")
# Alternatively Spliced / NOT FB-init / NO FB-init Isoforms
NOinit_as_NOinit <- setdiff(NOinit_as, NOinit_as_init)
write_tsv(NOinit_as_NOinit, "../temp_tables/NOinit_as_NOinit.tsv")



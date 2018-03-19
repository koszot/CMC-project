library(dplyr)
library(readr)
library(tidyr)

setwd("~/Desktop/MTA/CMC_project/scommune/")

genes <- read_tsv("scommune_genes.tsv", col_names = T, cols(.default = col_guess(), gene_id = col_character()))
devreg <- genes %>% filter(`FB-devreg` == T) %>% select(gene_id)
init <- genes %>% filter(`FB-init` == T) %>% select(gene_id)
as <- genes %>% filter(isoforms > 1) %>% select(gene_id)


stats <- tibble()


stats[1,1] <- "AS + DEVREG + INIT"
stats[1,2] <- nrow(intersect(as, intersect(devreg, init)))
stats[2,1] <- "DEVREG + AS"
stats[2,2] <- nrow(intersect(devreg, as)) - nrow(intersect(as, intersect(devreg, init)))
stats[3,1] <- "INIT + AS"
stats[3,2] <- nrow(intersect(init, as)) - nrow(intersect(as, intersect(devreg, init)))
stats[4,1] <- "DEVREG + INIT"
stats[4,2] <- nrow(intersect(devreg, init)) - nrow(intersect(as, intersect(devreg, init)))
stats[5,1] <- "only DEVREG"
stats[5,2] <- nrow(devreg) - stats[1,2] - stats[2,2] - stats[4,2]
stats[6,1] <- "only INIT"
stats[6,2] <- nrow(init) - stats[1,2] - stats[3,2] - stats[4,2]
stats[7,1] <- "only AS"
stats[7,2] <- nrow(as) - stats[1,2] - stats[2,2] - stats[3,2]



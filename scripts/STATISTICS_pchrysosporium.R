library(dplyr)
library(readr)
library(tidyr)

setwd("~/Desktop/MTA_CMC_project/pchrysosporium/")

genes <- read_tsv("pchrysosporium_genes.tsv", col_names = T, cols(.default = col_guess(), gene_id = col_character()))
devreg <- genes %>% filter(`FB-devreg` == T) %>% select(gene_id)
init <- genes %>% filter(`FB-init` == T) %>% select(gene_id)
as <- genes %>% filter(isoforms > 1) %>% select(gene_id)
stats <- as_tibble()

library(ASpli)

##### AS GENES #####

TxDb <- makeTxDbFromGFF(file="pchrysosporium_AS_genes.gtf", format="gtf")
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

TxDb <- makeTxDbFromGFF(file="pchrysosporium_FB_DEVREG_genes.gtf", format="gtf")
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

TxDb <- makeTxDbFromGFF(file="pchrysosporium_FB_INIT_genes.gtf", format="gtf")
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

##### STATS #####

stats[1,1] <- "All Genes"
stats[1,2] <- nrow(genes)
stats[1,3] <- NA

stats[2,1] <- "FB-devreg Genes"
stats[2,2] <- nrow(devreg)
stats[2,3] <- nrow(devreg) / nrow(genes) * 100

stats[3,1] <- "FB-init Genes"
stats[3,2] <- nrow(init)
stats[3,3] <- nrow(init) / nrow(genes) * 100

stats[4,1] <- "AS Genes"
stats[4,2] <- nrow(as)
stats[4,3] <- nrow(as) / nrow(genes) * 100

stats[5,1] <- "AS FB-devreg Genes"
stats[5,2] <- nrow(dplyr::intersect(devreg, as))
stats[5,3] <- nrow(dplyr::intersect(devreg, as)) / nrow(genes) * 100

stats[6,1] <- "AS FB-init Genes"
stats[6,2] <- nrow(dplyr::intersect(init, as))
stats[6,3] <- nrow(dplyr::intersect(init, as)) / nrow(genes) * 100

stats[7,1] <- "FB-deveg and FB-init Genes"
stats[7,2] <- nrow(dplyr::intersect(init, devreg))
stats[7,3] <- nrow(dplyr::intersect(init, devreg)) / nrow(genes) * 100

stats[8,1] <- "AS FB-deveg and AS FB-init Genes"
stats[8,2] <- nrow(dplyr::intersect(init, devreg) %>% dplyr::intersect(as))
stats[8,3] <- nrow(dplyr::intersect(init, devreg) %>% dplyr::intersect(as)) / nrow(genes) * 100

stats[9,1] <- "AS Genes -- ES bins"
stats[9,2] <- AS_AS_count$count[1]
stats[9,3] <- AS_AS_count$count[1] / sum(AS_AS_count$count) * 100
stats[10,1] <- "AS Genes -- IR bins"
stats[10,2] <- AS_AS_count$count[2]
stats[10,3] <- AS_AS_count$count[2] / sum(AS_AS_count$count) * 100
stats[11,1] <- "AS Genes -- Alt 5' SS bins"
stats[11,2] <- AS_AS_count$count[3]
stats[11,3] <- AS_AS_count$count[3] / sum(AS_AS_count$count) * 100
stats[12,1] <- "AS Genes -- Alt 3' SS bins"
stats[12,2] <- AS_AS_count$count[4]
stats[12,3] <- AS_AS_count$count[4] / sum(AS_AS_count$count) * 100

stats[13,1] <- "AS FB-devreg Genes -- ES bins"
stats[13,2] <- FB_devreg_AS_count$count[1]
stats[13,3] <- FB_devreg_AS_count$count[1] / sum(FB_devreg_AS_count$count) * 100
stats[14,1] <- "AS FB-devreg Genes -- IR bins"
stats[14,2] <- FB_devreg_AS_count$count[2]
stats[14,3] <- FB_devreg_AS_count$count[2] / sum(FB_devreg_AS_count$count) * 100
stats[15,1] <- "AS FB-devreg Genes -- Alt 5' SS bins"
stats[15,2] <- FB_devreg_AS_count$count[3]
stats[15,3] <- FB_devreg_AS_count$count[3] / sum(FB_devreg_AS_count$count) * 100
stats[16,1] <- "AS FB-devreg Genes -- Alt 3' SS bins"
stats[16,2] <- FB_devreg_AS_count$count[4]
stats[16,3] <- FB_devreg_AS_count$count[4] / sum(FB_devreg_AS_count$count) * 100

stats[17,1] <- "AS FB-init Genes -- ES bins"
stats[17,2] <- FB_init_AS_count$count[1]
stats[17,3] <- FB_init_AS_count$count[1] / sum(FB_init_AS_count$count) * 100
stats[18,1] <- "AS FB-init Genes -- IR bins"
stats[18,2] <- FB_init_AS_count$count[2]
stats[18,3] <- FB_init_AS_count$count[2] / sum(FB_init_AS_count$count) * 100
stats[19,1] <- "AS FB-init Genes -- Alt 5' SS bins"
stats[19,2] <- FB_init_AS_count$count[3]
stats[19,3] <- FB_init_AS_count$count[3] / sum(FB_init_AS_count$count) * 100
stats[20,1] <- "AS FB-init Genes -- Alt 3' SS bins"
stats[20,2] <- FB_init_AS_count$count[4]
stats[20,3] <- FB_init_AS_count$count[4] / sum(FB_init_AS_count$count) * 100

colnames(stats) <- c("names", "pchrysosporium", "%_pchrysosporium")
write_tsv(stats, "pchrysosporium_stats.tsv")

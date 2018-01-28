library(ASpli)

setwd("~/Desktop/MTA_CMC_project/scommune/")

stats <- read_tsv("scommune_stats.tsv")

##### AS GENES #####

TxDb <- makeTxDbFromGFF(file="scommune_AS_genes.gtf", format="gtf")
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

stats[7,1] <- "genes with AS ES bins"
stats[7,2] <- AS_count[1,2]
stats[8,1] <- "genes with AS IR bins"
stats[8,2] <- AS_count[2,2]
stats[9,1] <- "genes with AS ALT5SS bins"
stats[9,2] <- AS_count[3,2]
stats[10,1] <- "genes with AS ALT3SS bins"
stats[10,2] <- AS_count[4,2]

##### FB DEVREG GENES #####

TxDb <- makeTxDbFromGFF(file="scommune_FB_DEVREG_genes.gtf", format="gtf")
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

stats[11,1] <- "genes with FB-devreg ES bins"
stats[11,2] <- AS_count[1,2]
stats[12,1] <- "genes with FB-devreg IR bins"
stats[12,2] <- AS_count[2,2]
stats[13,1] <- "genes with FB-devreg ALT5SS bins"
stats[13,2] <- AS_count[3,2]
stats[14,1] <- "genes with FB-devreg ALT3SS bins"
stats[14,2] <- AS_count[4,2]

##### FB INIT GENES #####

TxDb <- makeTxDbFromGFF(file="scommune_FB_INIT_genes.gtf", format="gtf")
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

stats[15,1] <- "genes with FB-init ES bins"
stats[15,2] <- AS_count[1,2]
stats[16,1] <- "genes with FB-init IR bins"
stats[16,2] <- AS_count[2,2]
stats[17,1] <- "genes with FB-init ALT5SS bins"
stats[17,2] <- AS_count[3,2]
stats[18,1] <- "genes with FB-init ALT3SS bins"
stats[18,2] <- AS_count[4,2]

##### STATS #####

stats[1,3] <- NA
stats[2,3] <- stats[2,2] / stats[1,2] * 100
stats[3,3] <- stats[3,2] / stats[1,2] * 100
stats[4,3] <- stats[4,2] / stats[1,2] * 100
stats[5,3] <- stats[5,2] / stats[1,2] * 100
stats[6,3] <- stats[6,2] / stats[1,2] * 100
stats[7,3] <- stats[7,2] / sum(stats[7:10,2]) * 100
stats[8,3] <- stats[8,2] / sum(stats[7:10,2]) * 100
stats[9,3] <- stats[9,2] / sum(stats[7:10,2]) * 100
stats[10,3] <- stats[10,2] / sum(stats[7:10,2]) * 100
stats[11,3] <- stats[11,2] / sum(stats[11:14,2]) * 100
stats[12,3] <- stats[12,2] / sum(stats[11:14,2]) * 100
stats[13,3] <- stats[13,2] / sum(stats[11:14,2]) * 100
stats[14,3] <- stats[14,2] / sum(stats[11:14,2]) * 100
stats[15,3] <- stats[15,2] / sum(stats[15:18,2]) * 100
stats[16,3] <- stats[16,2] / sum(stats[15:18,2]) * 100
stats[17,3] <- stats[17,2] / sum(stats[15:18,2]) * 100
stats[18,3] <- stats[18,2] / sum(stats[15:18,2]) * 100

colnames(stats) <- c("names", "scommune", "%_scommune")
write_tsv(stats, "scommune_stats.tsv")

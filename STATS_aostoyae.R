library(tidyr)
library(readr)
library(dplyr)
library(stringr)

setwd("~/Desktop/alternative_splicing/aostoyae/")

isoforms <- read_tsv("aostoyae_CUFFdiff/isoforms.fpkm_tracking")     # betöltjük az isoform FPKM táblát

isoforms <- isoforms %>%
  select(transcript_id = tracking_id, gene_id, RMA_FPKM, VM_FPKM, 
         P1_FPKM, P2_C_FPKM, P2_S_FPKM, YFB_C_FPKM, 
         YFB_S_FPKM, FB_C_FPKM, FB_L_FPKM, FB_S_FPKM)            # kiválasztjuk a minket érdeklő oszlopokat

isoforms <- isoforms[,c(2,1,3:12)]
isoforms <- isoforms %>%
  arrange(gene_id, transcript_id)

##### ISOFORM FPKM RANKING #####

isoforms.rank <- as_tibble()

for (x in 1:length(unique((isoforms$gene_id)))) {
  temp <- isoforms[isoforms$gene_id %in% unique(isoforms$gene_id)[x],]
  temp$RMA_FPKM_Rank <- match(temp$RMA_FPKM,sort(temp$RMA_FPKM, decreasing = T))
  temp$VM_FPKM_Rank <- match(temp$VM_FPKM,sort(temp$VM_FPKM, decreasing = T))
  temp$P1_FPKM_Rank <- match(temp$P1_FPKM,sort(temp$P1_FPKM, decreasing = T))
  temp$P2_C_FPKM_Rank <- match(temp$P2_C_FPKM,sort(temp$P2_C_FPKM, decreasing = T))
  temp$P2_S_FPKM_Rank <- match(temp$P2_S_FPKM,sort(temp$P2_S_FPKM, decreasing = T))
  temp$YFB_C_FPKM_Rank <- match(temp$YFB_C_FPKM,sort(temp$YFB_C_FPKM, decreasing = T))
  temp$YFB_S_FPKM_Rank <- match(temp$YFB_S_FPKM,sort(temp$YFB_S_FPKM, decreasing = T))
  temp$FB_C_FPKM_Rank <- match(temp$FB_C_FPKM,sort(temp$FB_C_FPKM, decreasing = T))
  temp$FB_L_FPKM_Rank <- match(temp$FB_L_FPKM,sort(temp$FB_L_FPKM, decreasing = T))
  temp$FB_S_FPKM_Rank <- match(temp$FB_S_FPKM,sort(temp$FB_S_FPKM, decreasing = T))
  for(y in 1:length(temp$gene_id)) {
    temp$Rank_Sum[y] <- sum(temp$RMA_FPKM_Rank[y],
                         temp$VM_FPKM_Rank[y],
                         temp$P1_FPKM_Rank[y],
                         temp$P2_C_FPKM_Rank[y],
                         temp$P2_S_FPKM_Rank[y],
                         temp$YFB_C_FPKM_Rank[y],
                         temp$YFB_S_FPKM_Rank[y],
                         temp$FB_C_FPKM_Rank[y],
                         temp$FB_L_FPKM_Rank[y],
                         temp$FB_S_FPKM_Rank[y])
    temp$Overall_FPKM[y] <- sum(temp$RMA_FPKM[y],
                            temp$VM_FPKM[y],
                            temp$P1_FPKM[y],
                            temp$P2_C_FPKM[y],
                            temp$P2_S_FPKM[y],
                            temp$YFB_C_FPKM[y],
                            temp$YFB_S_FPKM[y],
                            temp$FB_C_FPKM[y],
                            temp$FB_L_FPKM[y],
                            temp$FB_S_FPKM[y])
  }
  temp <- arrange(temp, Rank_Sum, desc(Overall_FPKM))
  isoforms.rank <- rbind(isoforms.rank, temp)
  pb <- txtProgressBar(min = 1, max = length(unique((isoforms$gene_id))), style = 3)        
  setTxtProgressBar(pb, x, title = NULL, label = NULL) 
}

##### AS EVENTS #####

library(ASpli)

annotation <- read_tsv("~/Desktop/alternative_splicing/aostoyae/aostoyae_genome/aostoyae_AS_annotation.gtf",  col_names = c("chr", "maker","type", "start", "end", "att1", "strand", "att2", "attributes")) %>%
  separate(attributes, c("transcriptID_label", "transcriptID",
                         "geneID_label", "geneID"), sep = " ")

annotation$transcriptID <- annotation$transcriptID %>% 
  str_replace("\"", "") %>%
  str_replace("\";", "")              # transcriptID-k átalakítása

setwd("~/Desktop/AS_project/")

isoforms.rank.AS <- as_tibble()

for (x in 1:length(unique((isoforms.rank$gene_id)))) {
  temp <- isoforms.rank[isoforms.rank$gene_id %in% unique(isoforms.rank$gene_id)[x],]
  if (nrow(temp) == 1) next
  temp$ES_bins[1] <- 0
  temp$IR_bins[1] <- 0
  temp$ALT5SS_bins[1] <- 0
  temp$ALT3SS_bins[1] <- 0
  for (y in 2:length(temp$gene_id)) {
    temp_annotation <- rbind(annotation[annotation$transcriptID %in% temp$transcript_id[1],], annotation[annotation$transcriptID %in% temp$transcript_id[y],])
    temp_annotation$transcriptID <- temp_annotation$transcriptID %>% str_replace("^", "\"") %>% str_replace("$", "\";")
    temp_annotation <- unite(temp_annotation , attributes, 9:12, sep = " ")
    write.table(temp_annotation, file = "./temp_annotation.gtf", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
    TxDb <- makeTxDbFromGFF(file="./temp_annotation.gtf", format="gtf")
    features <- binGenome(TxDb)
    AS_count <- read_lines("./ASpli_binFeatures.log")
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
    temp$ES_bins[y] <- AS_count[1,2]
    temp$IR_bins[y] <- AS_count[2,2]
    temp$ALT5SS_bins[y] <- AS_count[3,2]
    temp$ALT3SS_bins[y] <- AS_count[4,2]
  }
  isoforms.rank.AS <- rbind(isoforms.rank.AS, temp)
  pb <- txtProgressBar(min = 1, max = length(unique((isoforms.rank$gene_id))), style = 3)        
  setTxtProgressBar(pb, x, title = NULL, label = NULL) 
}












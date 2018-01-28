library(tidyr)
library(readr)
library(dplyr)
library(stringr)

##### LOAD FILES #####

setwd("~/Desktop/MTA_CMC_project/rmellea/")

isoforms <- read_tsv("isoforms.fpkm_tracking", col_names = T, cols(.default = col_guess(), tracking_id = col_character(), gene_id = col_character()))   
isoforms <- isoforms %>%
  select(transcript_id = tracking_id, gene_id, VM_FPKM, P_FPKM, YFB_FPKM, FB_K_FPKM, FB_T_FPKM)            # kiválasztjuk a minket érdeklő oszlopokat
isoforms <- isoforms[,c(2,1,3:7)]
isoforms <- isoforms %>%
  arrange(gene_id, transcript_id)

annotation <- read_tsv("rmellea_AS_annotation.gtf",  col_names = c("chr", "maker","type", "start", "end", "att1", "strand", "att2", "attributes")) %>%
  separate(attributes, c("transcriptID_label", "transcriptID", "geneID_label", "geneID"), sep = " ")
annotation$transcriptID <- annotation$transcriptID %>% 
  str_replace("\"", "") %>%
  str_replace("\";", "")              # transcriptID-k átalakítása

genes <- read_tsv("genes.fpkm_tracking", col_names = T, cols(.default = col_guess(), tracking_id = col_character(), gene_id = col_character()))     # betöltjük az isoform FPKM táblát
genes <- genes %>% select(gene_id, VM_FPKM, P_FPKM, YFB_FPKM, FB_K_FPKM, FB_T_FPKM)            # kiválasztjuk a minket érdeklő oszlopokat

##### GENES - 4 FPKM P-FB_C #####

for (x in 1:length(genes$gene_id)) {
  if (max(genes[x,3:6]) >= 4) {
    genes$">= 4 FPKM"[x] <- T
  } else {
    genes$">= 4 FPKM"[x] <- F
  }
  pb <- txtProgressBar(min = 1, max = length(genes$gene_id), style = 3)        
  setTxtProgressBar(pb, x, title = NULL, label = NULL) 
}

##### GENES - DEVREG_T #####

for (x in 1:length(genes$gene_id)) {
  if (is.na((max(genes$P_FPKM[x], genes$YFB_FPKM[x], genes$FB_T_FPKM[x])/min(genes$P_FPKM[x], genes$YFB_FPKM[x], genes$FB_T_FPKM[x])))) {
    genes$DEVREG_T[x] <- NA 
  } else if ((max(genes$P_FPKM[x], genes$YFB_FPKM[x], genes$FB_T_FPKM[x])/min(genes$P_FPKM[x], genes$YFB_FPKM[x], genes$FB_T_FPKM[x])) >= 4) {
    genes$DEVREG_T[x] <- T
  } else {
    genes$DEVREG_T[x] <- F
  }
  pb <- txtProgressBar(min = 1, max = length(genes$gene_id), style = 3)        
  setTxtProgressBar(pb, x, title = NULL, label = NULL) 
}

##### GENES - DEVREG_K #####

for (x in 1:length(genes$gene_id)) {
  if (is.na((max(genes$P_FPKM[x], genes$YFB_FPKM[x], genes$FB_K_FPKM[x])/min(genes$P_FPKM[x], genes$YFB_FPKM[x], genes$FB_K_FPKM[x])))) {
    genes$DEVREG_K[x] <- NA 
  } else if ((max(genes$P_FPKM[x], genes$YFB_FPKM[x], genes$FB_K_FPKM[x])/min(genes$P_FPKM[x], genes$YFB_FPKM[x], genes$FB_K_FPKM[x])) >= 4) {
    genes$DEVREG_K[x] <- T
  } else {
    genes$DEVREG_K[x] <- F
  }
  pb <- txtProgressBar(min = 1, max = length(genes$gene_id), style = 3)        
  setTxtProgressBar(pb, x, title = NULL, label = NULL) 
}

##### GENES - DEVREG_FB #####

for (x in 1:length(genes$gene_id)) {
  if (is.na((max(genes$FB_T_FPKM[x], genes$FB_K_FPKM[x])/min(genes$FB_T_FPKM[x], genes$FB_K_FPKM[x])))) {
    genes$DEVREG_FB[x] <- NA 
  } else if ((max(genes$FB_T_FPKM[x], genes$FB_K_FPKM[x])/min(genes$FB_T_FPKM[x], genes$FB_K_FPKM[x])) >= 4) {
    genes$DEVREG_FB[x] <- T
  } else {
    genes$DEVREG_FB[x] <- F
  }
  pb <- txtProgressBar(min = 1, max = length(genes$gene_id), style = 3)        
  setTxtProgressBar(pb, x, title = NULL, label = NULL) 
}

##### GENES - FB-devreg #####

for (x in 1:length(genes$gene_id)) {
   if (genes$`>= 4 FPKM`[x] == T & (genes$DEVREG_T[x] == T | 
                                    genes$DEVREG_K[x] == T |
                                    genes$DEVREG_FB[x] == T) == T) {
    genes$"FB-devreg"[x] <- T
  } else {
    genes$"FB-devreg"[x] <- F
  }
  pb <- txtProgressBar(min = 1, max = length(genes$gene_id), style = 3)        
  setTxtProgressBar(pb, x, title = NULL, label = NULL) 
}

genes <- genes[,c(1:6,11)]

##### GENES - 4 FPKM VM-P #####

for (x in 1:length(genes$gene_id)) {
  if (max(genes[x,2:3]) >= 4) {
    genes$">= 4 FPKM"[x] <- T
  } else {
    genes$">= 4 FPKM"[x] <- F
  }
  pb <- txtProgressBar(min = 1, max = length(genes$gene_id), style = 3)        
  setTxtProgressBar(pb, x, title = NULL, label = NULL) 
}

##### GENES - FB-init #####

for (x in 1:length(genes$gene_id)) {
  if (genes$`>= 4 FPKM`[x] == T & ((genes$P_FPKM[x] / genes$VM_FPKM[x]) >= 4) == T) {
    genes$"FB-init"[x] <- T
  } else {
    genes$"FB-init"[x] <- F
  }
  pb <- txtProgressBar(min = 1, max = length(genes$gene_id), style = 3)        
  setTxtProgressBar(pb, x, title = NULL, label = NULL) 
}

genes <- genes[,c(1:7,9)]

##### ISOFORM - FPKM RANKING #####

isoforms.rank <- as_tibble()

for (x in 1:length(unique((isoforms$gene_id)))) {
  temp <- isoforms[isoforms$gene_id %in% unique(isoforms$gene_id)[x],]
  temp$VM_FPKM_Rank <- match(temp$VM_FPKM,sort(temp$VM_FPKM, decreasing = T))
  temp$P_FPKM_Rank <- match(temp$P_FPKM,sort(temp$P_FPKM, decreasing = T))
  temp$YFB_FPKM_Rank <- match(temp$YFB_FPKM,sort(temp$YFB_FPKM, decreasing = T))
  temp$FB_K_FPKM_Rank <- match(temp$FB_K_FPKM,sort(temp$FB_K_FPKM, decreasing = T))
  temp$FB_T_FPKM_Rank <- match(temp$FB_T_FPKM,sort(temp$FB_T_FPKM, decreasing = T))
  for(y in 1:length(temp$gene_id)) {
    temp$Rank_Sum[y] <- sum(temp$VM_FPKM_Rank[y],
                         temp$P_FPKM_Rank[y],
                         temp$YFB_FPKM_Rank[y],
                         temp$FB_K_FPKM_Rank[y],
                         temp$FB_T_FPKM_Rank[y])
    temp$Overall_FPKM[y] <- sum(temp$VM_FPKM[y],
                            temp$P_FPKM[y],
                            temp$YFB_FPKM[y],
                            temp$FB_K_FPKM[y],
                            temp$FB_T_FPKM[y])
  }
  temp <- arrange(temp, Rank_Sum, desc(Overall_FPKM))
  isoforms.rank <- rbind(isoforms.rank, temp)
  pb <- txtProgressBar(min = 1, max = length(unique((isoforms$gene_id))), style = 3)        
  setTxtProgressBar(pb, x, title = NULL, label = NULL) 
}

##### GENES - isoforms #####

genes.AS <- as_tibble()

for (x in 1:length(genes$gene_id)) {
  temp <- genes[genes$gene_id %in% genes$gene_id[x],]
  temp$isoforms <- nrow(isoforms.rank[isoforms.rank$gene_id %in% genes$gene_id[x],])
  genes.AS <- rbind(genes.AS, temp)
  pb <- txtProgressBar(min = 1, max = length(genes$gene_id), style = 3)        
  setTxtProgressBar(pb, x, title = NULL, label = NULL) 
}

##### ANNOTATION CORRECTION #####

annotation.corrected <- as_tibble()

for (x in 1:length(isoforms.rank$transcript_id)) {
  temp <- annotation[annotation$transcriptID %in% isoforms.rank$transcript_id[x],]
  annotation.corrected <- rbind(annotation.corrected, temp)
  pb <- txtProgressBar(min = 1, max = length(isoforms.rank$transcript_id), style = 3)        
  setTxtProgressBar(pb, x, title = NULL, label = NULL) 
}

sum(unique(annotation.corrected$transcriptID) != isoforms.rank$transcript_id)

##### WRITE FILES #####

write_tsv(isoforms.rank, "rmellea_isoforms.tsv")

write_tsv(genes.AS, "rmellea_genes.tsv")

annotation.corrected$transcriptID <- annotation.corrected$transcriptID %>% 
  str_replace("^", "\"") %>%
  str_replace("$", "\";")              # transcriptID-k átalakítása

merged <- unite(annotation.corrected, attributes, 9:12, sep = " ")

write.table(merged, file = "rmellea_corrected_annotation.gtf", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")



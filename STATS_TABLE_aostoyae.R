library(tidyr)
library(readr)
library(dplyr)
library(stringr)

##### LOAD FILES #####

setwd("~/Desktop/MTA_CMC_project/aostoyae/")

isoforms <- read_tsv("isoforms.fpkm_tracking")     # betöltjük az isoform FPKM táblát
isoforms <- isoforms %>%
  select(transcript_id = tracking_id, gene_id, VM_FPKM, 
         P1_FPKM, P2_C_FPKM, P2_S_FPKM, YFB_C_FPKM, 
         YFB_S_FPKM, FB_C_FPKM, FB_L_FPKM, FB_S_FPKM)            # kiválasztjuk a minket érdeklő oszlopokat
isoforms <- isoforms[,c(2,1,3:11)]
isoforms <- isoforms %>%
  arrange(gene_id, transcript_id)

annotation <- read_tsv("aostoyae_AS_annotation.gtf",  col_names = c("chr", "maker","type", "start", "end", "att1", "strand", "att2", "attributes")) %>%
  separate(attributes, c("transcriptID_label", "transcriptID", "geneID_label", "geneID"), sep = " ")
annotation$transcriptID <- annotation$transcriptID %>% 
  str_replace("\"", "") %>%
  str_replace("\";", "")              # transcriptID-k átalakítása

genes <- read_tsv("genes.fpkm_tracking")     # betöltjük az isoform FPKM táblát
genes <- genes %>%
  select(gene_id, VM_FPKM, 
         P1_FPKM, P2_C_FPKM, P2_S_FPKM, YFB_C_FPKM, 
         YFB_S_FPKM, FB_C_FPKM, FB_L_FPKM, FB_S_FPKM)            # kiválasztjuk a minket érdeklő oszlopokat

##### GENES - 4 FPKM P-FB #####

for (x in 1:length(genes$gene_id)) {
  if (max(genes[x,3:10]) >= 4) {
    genes$">= 4 FPKM"[x] <- T
  } else {
    genes$">= 4 FPKM"[x] <- F
  }
  pb <- txtProgressBar(min = 1, max = length(genes$gene_id), style = 3)        
  setTxtProgressBar(pb, x, title = NULL, label = NULL) 
}

##### GENES - DEVREG_S #####

for (x in 1:length(genes$gene_id)) {
  if (is.na((max(genes$P1_FPKM[x], genes$P2_S_FPKM[x], genes$YFB_S_FPKM[x], genes$FB_S_FPKM[x])/min(genes$P1_FPKM[x], genes$P2_S_FPKM[x], genes$YFB_S_FPKM[x], genes$FB_S_FPKM[x])))) {
    genes$DEVREG_S[x] <- NA 
  } else if ((max(genes$P1_FPKM[x], genes$P2_S_FPKM[x], genes$YFB_S_FPKM[x], genes$FB_S_FPKM[x])/min(genes$P1_FPKM[x], genes$P2_S_FPKM[x], genes$YFB_S_FPKM[x], genes$FB_S_FPKM[x])) >= 4) {
    genes$DEVREG_S[x] <- T
  } else {
    genes$DEVREG_S[x] <- F
  }
  pb <- txtProgressBar(min = 1, max = length(genes$gene_id), style = 3)        
  setTxtProgressBar(pb, x, title = NULL, label = NULL) 
}

##### GENES - DEVREG_C #####

for (x in 1:length(genes$gene_id)) {
  if (is.na((max(genes$P1_FPKM[x], genes$P2_C_FPKM[x], genes$YFB_C_FPKM[x], genes$FB_C_FPKM[x])/min(genes$P1_FPKM[x], genes$P2_C_FPKM[x], genes$YFB_C_FPKM[x], genes$FB_C_FPKM[x])))) {
    genes$DEVREG_C[x] <- NA 
  } else if ((max(genes$P1_FPKM[x], genes$P2_C_FPKM[x], genes$YFB_C_FPKM[x], genes$FB_C_FPKM[x])/min(genes$P1_FPKM[x], genes$P2_C_FPKM[x], genes$YFB_C_FPKM[x], genes$FB_C_FPKM[x])) >= 4) {
    genes$DEVREG_C[x] <- T
  } else {
    genes$DEVREG_C[x] <- F
  }
  pb <- txtProgressBar(min = 1, max = length(genes$gene_id), style = 3)        
  setTxtProgressBar(pb, x, title = NULL, label = NULL) 
}

##### GENES - DEVREG_L #####

for (x in 1:length(genes$gene_id)) {
  if (is.na((max(genes$P1_FPKM[x], genes$P2_C_FPKM[x], genes$YFB_C_FPKM[x], genes$FB_L_FPKM[x])/min(genes$P1_FPKM[x], genes$P2_C_FPKM[x], genes$YFB_C_FPKM[x], genes$FB_L_FPKM[x])))) {
    genes$DEVREG_L[x] <- NA 
  } else if ((max(genes$P1_FPKM[x], genes$P2_C_FPKM[x], genes$YFB_C_FPKM[x], genes$FB_L_FPKM[x])/min(genes$P1_FPKM[x], genes$P2_C_FPKM[x], genes$YFB_C_FPKM[x], genes$FB_L_FPKM[x])) >= 4) {
    genes$DEVREG_L[x] <- T
  } else {
    genes$DEVREG_L[x] <- F
  }
  pb <- txtProgressBar(min = 1, max = length(genes$gene_id), style = 3)        
  setTxtProgressBar(pb, x, title = NULL, label = NULL) 
}

##### GENES - DEVREG_P2 #####

for (x in 1:length(genes$gene_id)) {
  if (is.na((max(genes$P2_C_FPKM[x], genes$P2_S_FPKM[x])/min(genes$P2_C_FPKM[x], genes$P2_S_FPKM[x])))) {
    genes$DEVREG_P2[x] <- NA 
  } else if ((max(genes$P2_C_FPKM[x], genes$P2_S_FPKM[x])/min(genes$P2_C_FPKM[x], genes$P2_S_FPKM[x])) >= 4) {
    genes$DEVREG_P2[x] <- T
  } else {
    genes$DEVREG_P2[x] <- F
  }
  pb <- txtProgressBar(min = 1, max = length(genes$gene_id), style = 3)        
  setTxtProgressBar(pb, x, title = NULL, label = NULL) 
}

##### GENES - DEVREG_YFB #####

for (x in 1:length(genes$gene_id)) {
  if (is.na((max(genes$YFB_C_FPKM[x], genes$YFB_S_FPKM[x])/min(genes$YFB_C_FPKM[x], genes$YFB_S_FPKM[x])))) {
    genes$DEVREG_YFB[x] <- NA 
  } else if ((max(genes$YFB_C_FPKM[x], genes$YFB_S_FPKM[x])/min(genes$YFB_C_FPKM[x], genes$YFB_S_FPKM[x])) >= 4) {
    genes$DEVREG_YFB[x] <- T
  } else {
    genes$DEVREG_YFB[x] <- F
  }
  pb <- txtProgressBar(min = 1, max = length(genes$gene_id), style = 3)        
  setTxtProgressBar(pb, x, title = NULL, label = NULL) 
}

##### GENES - DEVREG_FB #####

for (x in 1:length(genes$gene_id)) {
  if (is.na((max(genes$FB_C_FPKM[x], genes$FB_S_FPKM[x], genes$FB_L_FPKM[x])/min(genes$FB_C_FPKM[x], genes$FB_S_FPKM[x], genes$FB_L_FPKM[x])))) {
    genes$DEVREG_FB[x] <- NA 
  } else if ((max(genes$FB_C_FPKM[x], genes$FB_S_FPKM[x], genes$FB_L_FPKM[x])/min(genes$FB_C_FPKM[x], genes$FB_S_FPKM[x], genes$FB_L_FPKM[x])) >= 4) {
    genes$DEVREG_FB[x] <- T
  } else {
    genes$DEVREG_FB[x] <- F
  }
  pb <- txtProgressBar(min = 1, max = length(genes$gene_id), style = 3)        
  setTxtProgressBar(pb, x, title = NULL, label = NULL) 
}

##### GENES - FB-devreg #####

for (x in 1:length(genes$gene_id)) {
   if (genes$`>= 4 FPKM`[x] == T & (genes$DEVREG_S[x] == T | 
                                    genes$DEVREG_C[x] == T |
                                    genes$DEVREG_L[x] == T | 
                                    genes$DEVREG_P2[x] == T |
                                    genes$DEVREG_YFB[x] == T |
                                    genes$DEVREG_FB[x] == T) == T) {
    genes$"FB-devreg"[x] <- T
  } else {
    genes$"FB-devreg"[x] <- F
  }
  pb <- txtProgressBar(min = 1, max = length(genes$gene_id), style = 3)        
  setTxtProgressBar(pb, x, title = NULL, label = NULL) 
}

genes <- genes[,c(1:10,18)]

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
  if (genes$`>= 4 FPKM`[x] == T & ((genes$P1_FPKM[x] / genes$VM_FPKM[x]) >= 4) == T) {
    genes$"FB-init"[x] <- T
  } else {
    genes$"FB-init"[x] <- F
  }
  pb <- txtProgressBar(min = 1, max = length(genes$gene_id), style = 3)        
  setTxtProgressBar(pb, x, title = NULL, label = NULL) 
}

genes <- genes[,c(1:11,13)]






##### ISOFORM - FPKM RANKING #####

isoforms.rank <- as_tibble()

for (x in 1:length(unique((isoforms$gene_id)))) {
  temp <- isoforms[isoforms$gene_id %in% unique(isoforms$gene_id)[x],]
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
    temp$Rank_Sum[y] <- sum(temp$VM_FPKM_Rank[y],
                         temp$P1_FPKM_Rank[y],
                         temp$P2_C_FPKM_Rank[y],
                         temp$P2_S_FPKM_Rank[y],
                         temp$YFB_C_FPKM_Rank[y],
                         temp$YFB_S_FPKM_Rank[y],
                         temp$FB_C_FPKM_Rank[y],
                         temp$FB_L_FPKM_Rank[y],
                         temp$FB_S_FPKM_Rank[y])
    temp$Overall_FPKM[y] <- sum(temp$VM_FPKM[y],
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

##### ANNOTATION CORRECTION #####

annotation.corrected <- as_tibble()

for (x in 1:length(isoforms.rank$transcript_id)) {
  temp <- annotation[annotation$transcriptID %in% isoforms.rank$transcript_id[x],]
  annotation.corrected <- rbind(annotation.corrected, temp)
  pb <- txtProgressBar(min = 1, max = length(isoforms.rank$transcript_id), style = 3)        
  setTxtProgressBar(pb, x, title = NULL, label = NULL) 
}



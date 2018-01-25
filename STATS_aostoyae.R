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

###########################
#     ISOFORM COUNTER     #
###########################

alternative_splicing_associated_genes <- function(isoforms.fpkm_FUN) {
  geneids <- unique(isoforms.fpkm_FUN$gene_id)                                           # kimentjük az egyedi geneID-kat
  isoforms.fpkm.AS <- data.frame()                                                       # megcsináljuk az üres data framet
  for (x in 1:length(geneids)) {                                                         # egy ciklussal végigmegyünk az gén ID-ken
    geneID <- geneids[x]                                                                 # kiszedjük az aktuális geneID-t
    isoforms_temp <- isoforms.fpkm_FUN[isoforms.fpkm_FUN$gene_id %in% geneID,]           # aktuális geneID-hoz tartozó isoformák kiszedése
    isoforms_temp$N_isoforms <- nrow(isoforms_temp)                                    # az adott génhez tartozó isoformák leszámolása
    isoforms.fpkm.AS <- rbind(isoforms.fpkm.AS, isoforms_temp)       
    pb <- txtProgressBar(min = 1, max = length(geneids), style = 3)                      # progress bar
    setTxtProgressBar(pb, x, title = NULL, label = NULL) 
  }
  return(isoforms.fpkm.AS)
}

isoforms.AS <- alternative_splicing_associated_genes(isoforms)       # lefuttatjuk az AS szűrést

##################
#     4 FPKM     #
##################

isoforms.AS.FPKM <- as_tibble()

for (x in 1:length(unique((isoforms.AS$gene_id)))) {
  temp <- isoforms.AS[isoforms.AS$gene_id %in% unique((isoforms.AS$gene_id))[x],]
  if (sum(temp$VM_FPKM) >= 4 | 
      sum(temp$P1_FPKM) >= 4 | 
      sum(temp$P2_C_FPKM) >= 4 | 
      sum(temp$P2_S_FPKM) >= 4 | 
      sum(temp$YFB_C_FPKM) >= 4 |
      sum(temp$YFB_S_FPKM) >= 4 | 
      sum(temp$FB_C_FPKM) >= 4 | 
      sum(temp$FB_L_FPKM) >= 4 | 
      sum(temp$FB_S_FPKM) >= 4) {
    temp$gene_4FPKM <- T
  } else {
    temp$gene_4FPKM <- F
  }
  isoforms.AS.FPKM <- rbind(isoforms.AS.FPKM, temp)
  pb <- txtProgressBar(min = 1, max = length(unique((isoforms.AS$gene_id))), style = 3)        
  setTxtProgressBar(pb, x, title = NULL, label = NULL) 
}

######################
#     FOLD CHANGE    #
######################

isoforms.AS.FPKM.FC <- as_tibble()

for (x in 1:length(unique((isoforms.AS.FPKM$gene_id)))) {
  temp <- isoforms.AS.FPKM[isoforms.AS.FPKM$gene_id %in% unique((isoforms.AS.FPKM$gene_id))[x],]
  temp$gene_VMtoP_FC <- sum(temp$P1_FPKM) / sum(temp$VM_FPKM)
  min_P_FB <- min(sum(temp$P1_FPKM), 
                  sum(temp$P2_C_FPKM),
                  sum(temp$P2_S_FPKM),
                  sum(temp$YFB_C_FPKM),
                  sum(temp$YFB_S_FPKM),
                  sum(temp$FB_C_FPKM),
                  sum(temp$FB_L_FPKM),
                  sum(temp$FB_S_FPKM))
  max_P_FB <- max(sum(temp$P1_FPKM), 
                  sum(temp$P2_C_FPKM),
                  sum(temp$P2_S_FPKM),
                  sum(temp$YFB_C_FPKM),
                  sum(temp$YFB_S_FPKM),
                  sum(temp$FB_C_FPKM),
                  sum(temp$FB_L_FPKM),
                  sum(temp$FB_S_FPKM))
  temp$gene_PtoFB_FC <- max_P_FB / min_P_FB    
  isoforms.AS.FPKM.FC <- rbind(isoforms.AS.FPKM.FC, temp)
  pb <- txtProgressBar(min = 1, max = length(unique((isoforms.AS.FPKM$gene_id))), style = 3)        
  setTxtProgressBar(pb, x, title = NULL, label = NULL) 
}

###############################
#     ISOFORM SIGNIFICANCY    #
###############################


isoforms.AS.FPKM.FC.Sign <- as_tibble()

for (y in 1:length(unique((isoforms.AS.FPKM.FC$gene_id)))) {
  # kiszedjük az aktuális gene_id-hoz tartozó részt
  temp <- isoforms.AS.FPKM.FC[isoforms.AS.FPKM.FC$gene_id %in% unique((isoforms.AS.FPKM.FC$gene_id))[y],]
  # megvizsgáljuk benne egymás után a sorokat
  for (x in 1:length(temp$gene_id)) {
    temp$isoform_Sign_VM[x] <- temp$VM_FPKM[x] >= (sum(temp$VM_FPKM) / temp$N_isoforms[1])
    temp$isoform_Sign_P1[x] <- temp$P1_FPKM[x] >= (sum(temp$P1_FPKM) / temp$N_isoforms[1])
    temp$isoform_Sign_P2_C[x] <- temp$P2_C_FPKM[x] >= (sum(temp$P2_C_FPKM) / temp$N_isoforms[1])
    temp$isoform_Sign_P2_S[x] <- temp$P2_S_FPKM[x] >= (sum(temp$P2_S_FPKM) / temp$N_isoforms[1])
    temp$isoform_Sign_YFB_C[x] <- temp$YFB_C_FPKM[x] >= (sum(temp$YFB_C_FPKM) / temp$N_isoforms[1])
    temp$isoform_Sign_YFB_S[x] <- temp$YFB_S_FPKM[x] >= (sum(temp$YFB_S_FPKM) / temp$N_isoforms[1])
    temp$isoform_Sign_FB_C[x] <- temp$FB_C_FPKM[x] >= (sum(temp$FB_C_FPKM) / temp$N_isoforms[1])
    temp$isoform_Sign_FB_L[x] <- temp$FB_L_FPKM[x] >= (sum(temp$FB_L_FPKM) / temp$N_isoforms[1])
    temp$isoform_Sign_FB_S[x] <- temp$FB_S_FPKM[x] >= (sum(temp$FB_S_FPKM) / temp$N_isoforms[1])
    temp$isoform_maxSign[x] <- sum(temp$isoform_Sign_VM[x], 
                                   temp$isoform_Sign_P1[x],
                                   temp$isoform_Sign_P2_C[x],
                                   temp$isoform_Sign_P2_S[x],
                                   temp$isoform_Sign_YFB_C[x],
                                   temp$isoform_Sign_YFB_S[x],
                                   temp$isoform_Sign_FB_C[x],
                                   temp$isoform_Sign_FB_L[x],
                                   temp$isoform_Sign_FB_S[x])
  }
  temp <- temp %>%
    arrange(desc(isoform_maxSign))
  isoforms.AS.FPKM.FC.Sign <- rbind(isoforms.AS.FPKM.FC.Sign, temp)
  pb <- txtProgressBar(min = 1, max = length(unique((isoforms.AS.FPKM.FC$gene_id))), style = 3)        
  setTxtProgressBar(pb, y, title = NULL, label = NULL) 
}










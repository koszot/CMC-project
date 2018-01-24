library(tidyr)
library(readr)
library(dplyr)
library(stringr)

setwd("~/Desktop/alternative_splicing/aostoyae/")

isoforms <- read_tsv("aostoyae_CUFFdiff/isoforms.fpkm_tracking")     # betöltjük az isoform FPKM táblát
genes <- read_tsv("aostoyae_CUFFdiff/genes.fpkm_tracking")           # betöltjük az genes FPKM táblát

isoforms <- isoforms %>%
  select(transcript_id = tracking_id, gene_id, RMA_FPKM, VM_FPKM, 
         P1_FPKM, P2_C_FPKM, P2_S_FPKM, YFB_C_FPKM, 
         YFB_S_FPKM, FB_C_FPKM, FB_L_FPKM, FB_S_FPKM)            # kiválasztjuk a minket érdeklő oszlopokat

genes <- genes %>%
  select(tracking_id, gene_id, RMA_FPKM, VM_FPKM, 
         P1_FPKM, P2_C_FPKM, P2_S_FPKM, YFB_C_FPKM, 
         YFB_S_FPKM, FB_C_FPKM, FB_L_FPKM, FB_S_FPKM)            # kiválasztjuk a minket érdeklő oszlopokat

###################################################
#     Alternaive Splicing gének és izoformáik     #
###################################################

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




#####################
#     R értékek     #
#####################

isoforms.AS.R <- data.frame()
ID <- unique(isoforms.AS$gene_id)
for (x in 1:length(ID)) {
    isoforms_temp <- isoforms.AS[isoforms.AS$gene_id %in% ID[x],]  
    isoforms_temp$RMA_Ri <- isoforms_temp$RMA_FPKM / sum(isoforms_temp$RMA_FPKM)
    isoforms_temp$VM_Ri <- isoforms_temp$VM_FPKM / sum(isoforms_temp$VM_FPKM)
    isoforms_temp$P1_Ri <- isoforms_temp$P1_FPKM / sum(isoforms_temp$P1_FPKM)
    isoforms_temp$P2_C_Ri <- isoforms_temp$P2_C_FPKM / sum(isoforms_temp$P2_C_FPKM)
    isoforms_temp$P2_S_Ri <- isoforms_temp$P2_S_FPKM / sum(isoforms_temp$P2_S_FPKM)
    isoforms_temp$YFB_C_Ri <- isoforms_temp$YFB_C_FPKM / sum(isoforms_temp$YFB_C_FPKM)
    isoforms_temp$YFB_S_Ri <- isoforms_temp$YFB_S_FPKM / sum(isoforms_temp$YFB_S_FPKM)
    isoforms_temp$FB_C_Ri <- isoforms_temp$FB_C_FPKM / sum(isoforms_temp$FB_C_FPKM)
    isoforms_temp$FB_L_Ri <- isoforms_temp$FB_L_FPKM / sum(isoforms_temp$FB_L_FPKM)
    isoforms_temp$FB_S_Ri <- isoforms_temp$FB_S_FPKM / sum(isoforms_temp$FB_S_FPKM)
    isoforms.AS.R <- rbind(isoforms.AS.R, isoforms_temp)
    pb <- txtProgressBar(min = 1, max = length(ID), style = 3)                      # progress bar
    setTxtProgressBar(pb, x, title = NULL, label = NULL) 
}

###############################
#     R érték fold change     #
###############################

for (x in 1:length(isoforms.AS.R$transcript_id)) {
  isoforms.AS.R$VM_Ri_to_P1_Ri_FoldChange[x] <- isoforms.AS.R[x,15] / isoforms.AS.R[x,16]          
  isoforms.AS.R$P1_Ri_to_FB_Ri_FoldChange[x] <- max(isoforms.AS.R[x,16:23]) / min(isoforms.AS.R[x,16:23])        
  isoforms.AS.R$VM_Ri_to_FB_Ri_FoldChange[x] <- max(isoforms.AS.R[x,15:23]) / min(isoforms.AS.R[x,15:23])        
  pb <- txtProgressBar(min = 1, max = length(isoforms.AS.R$transcript_id), style = 3)             
  setTxtProgressBar(pb, x, title = NULL, label = NULL)
}

isoforms.AS.R.arranged <- isoforms.AS.R %>%
  arrange(gene_id, transcript_id)



genes[genes$gene_id %in% "AROS_05942",] 










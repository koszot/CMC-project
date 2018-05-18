library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)

setwd("~/Desktop/MTA/CMC_project/ccinerea/")

# betöltjük a fájlokat
genes <- read_tsv("ccinerea_genes.tsv", col_names = T, cols(.default = col_guess(), gene_id = col_character()))
isoforms <- read_tsv("ccinerea_isoforms.tsv", col_names = T, cols(.default = col_guess(), gene_id = col_character()))

# kiválasztjuk azokat amik részt vesznek AS-ben
genes <- genes %>% filter(isoforms > 1)
isoforms <- isoforms %>% filter(`GENE-isoforms` > 1)

##### GENES TABLE #####

# megnézzük melyik gének lépek meg a 4 FPKM-et FB-devregben és kiszámoljuk a FC-eket
for (x in 1:length(genes$gene_id)) {
  if (max(genes[x,3:10]) >= 4) {
    genes$FB_devreg_4FPKM[x] <- T
  } else {
    genes$FB_devreg_4FPKM[x] <- F
  }
  genes$DEVREG_T[x] <- max(genes$H_FPKM[x], genes$P1_FPKM[x], genes$P2_FPKM[x], genes$YFB_T_FPKM[x], genes$FB_T_FPKM[x])/min(genes$H_FPKM[x], genes$P1_FPKM[x], genes$P2_FPKM[x], genes$YFB_T_FPKM[x], genes$FB_T_FPKM[x])
  genes$DEVREG_K[x] <- max(genes$H_FPKM[x], genes$P1_FPKM[x], genes$P2_FPKM[x], genes$YFB_K_FPKM[x], genes$FB_KL_FPKM[x])/min(genes$H_FPKM[x], genes$P1_FPKM[x], genes$P2_FPKM[x], genes$YFB_K_FPKM[x], genes$FB_KL_FPKM[x])
  genes$DEVREG_L[x] <- max(genes$H_FPKM[x], genes$P1_FPKM[x], genes$P2_FPKM[x], genes$YFB_L_FPKM[x], genes$FB_KL_FPKM[x])/min(genes$H_FPKM[x], genes$P1_FPKM[x], genes$P2_FPKM[x], genes$YFB_L_FPKM[x], genes$FB_KL_FPKM[x])
  genes$DEVREG_YFB[x] <- max(genes$YFB_K_FPKM[x], genes$YFB_L_FPKM[x], genes$YFB_T_FPKM[x])/min(genes$YFB_K_FPKM[x], genes$YFB_L_FPKM[x], genes$YFB_T_FPKM[x])
  genes$DEVREG_FB[x] <- max(genes$FB_T_FPKM[x], genes$FB_KL_FPKM[x])/min(genes$FB_T_FPKM[x], genes$FB_KL_FPKM[x])
}

# megnézzük melyik gének lépek meg a 4 FPKM-et FB-initben
for (x in 1:length(genes$gene_id)) {
  if (max(genes[x,2:4]) >= 4) {
    genes$FB_init_4FPKM[x] <- T
  } else {
    genes$FB_init_4FPKM[x] <- F
  }
  genes$INIT_FC[x] <- max(genes$H_FPKM[x], genes$P1_FPKM[x]) / genes$VM_FPKM[x]
}

# megkeressük a maximum FC-et
for (x in 1:length(genes$gene_id)) {
  genes$gene_FC[x] <- max(genes$DEVREG_T[x], genes$DEVREG_K[x], genes$DEVREG_L[x], genes$DEVREG_YFB[x], genes$DEVREG_FB[x], genes$INIT_FC[x])
  genes$DEVREG_FC[x] <- max(genes$DEVREG_T[x], genes$DEVREG_K[x], genes$DEVREG_L[x], genes$DEVREG_YFB[x], genes$DEVREG_FB[x])
}

# kiszűrjük azokat a géneket amik nem érik el a 4FPKM-et
genes <- genes %>% filter(FB_devreg_4FPKM == T | FB_init_4FPKM == T)

# takarítunk az oszlopok között, csak az marad ami kell
genes <- genes %>% select(gene_id, gene_FC)

##### ISOFORMS TABLE #####

# megnézzük melyik gének lépek meg a 4 FPKM-et FB-devregben és kiszámoljuk a FC-eket
for (x in 1:length(isoforms$gene_id)) {
  if (max(isoforms[x,4:11]) >= 4) {
    isoforms$FB_devreg_4FPKM[x] <- T
  } else {
    isoforms$FB_devreg_4FPKM[x] <- F
  }
  isoforms$DEVREG_T[x] <- max(isoforms$H_FPKM[x], isoforms$P1_FPKM[x], isoforms$P2_FPKM[x], isoforms$YFB_T_FPKM[x], isoforms$FB_T_FPKM[x])/min(isoforms$H_FPKM[x], isoforms$P1_FPKM[x], isoforms$P2_FPKM[x], isoforms$YFB_T_FPKM[x], isoforms$FB_T_FPKM[x])
  isoforms$DEVREG_K[x] <- max(isoforms$H_FPKM[x], isoforms$P1_FPKM[x], isoforms$P2_FPKM[x], isoforms$YFB_K_FPKM[x], isoforms$FB_KL_FPKM[x])/min(isoforms$H_FPKM[x], isoforms$P1_FPKM[x], isoforms$P2_FPKM[x], isoforms$YFB_K_FPKM[x], isoforms$FB_KL_FPKM[x])
  isoforms$DEVREG_L[x] <- max(isoforms$H_FPKM[x], isoforms$P1_FPKM[x], isoforms$P2_FPKM[x], isoforms$YFB_L_FPKM[x], isoforms$FB_KL_FPKM[x])/min(isoforms$H_FPKM[x], isoforms$P1_FPKM[x], isoforms$P2_FPKM[x], isoforms$YFB_L_FPKM[x], isoforms$FB_KL_FPKM[x])
  isoforms$DEVREG_YFB[x] <- max(isoforms$YFB_K_FPKM[x], isoforms$YFB_L_FPKM[x], isoforms$YFB_T_FPKM[x])/min(isoforms$YFB_K_FPKM[x], isoforms$YFB_L_FPKM[x], isoforms$YFB_T_FPKM[x])
  isoforms$DEVREG_FB[x] <- max(isoforms$FB_T_FPKM[x], isoforms$FB_KL_FPKM[x])/min(isoforms$FB_T_FPKM[x], isoforms$FB_KL_FPKM[x])
}

# megnézzük melyik gének lépek meg a 4 FPKM-et FB-initben
for (x in 1:length(isoforms$gene_id)) {
  if (max(isoforms[x,3:5]) >= 4) {
    isoforms$FB_init_4FPKM[x] <- T
  } else {
    isoforms$FB_init_4FPKM[x] <- F
  }
  isoforms$INIT_FC[x] <- max(isoforms$H_FPKM[x], isoforms$P1_FPKM[x]) / isoforms$VM_FPKM[x]
}

# megkeressük a maximum FC-et
for (x in 1:length(isoforms$gene_id)) {
  isoforms$isoform_FC[x] <- max(isoforms$DEVREG_T[x], isoforms$DEVREG_K[x], isoforms$DEVREG_L[x], isoforms$DEVREG_YFB[x], isoforms$DEVREG_FB[x], isoforms$INIT_FC[x])
}

# kiszűrjük azokat a géneket amik nem érik el a 4FPKM-et
isoforms <- isoforms %>% filter(FB_devreg_4FPKM == T | FB_init_4FPKM == T)

# devreg factorozás
for(x in 1:length(isoforms$gene_id)) {
  if (isoforms$`GENE-devreg-s-stricto`[x] == T & isoforms$`ISOFORM-devreg-s-stricto`[x] == T) {
    isoforms$devreg_factor[x] <- "Gene devreg-s-stricto / Isoform devreg-s-stricto"
  } else if (isoforms$`GENE-devreg-s-stricto`[x] == T & isoforms$`ISOFORM-devreg-s-stricto`[x] == F) {
    isoforms$devreg_factor[x] <- "Gene devreg-s-stricto / Isoform NOT devreg-s-stricto"
  } else if (isoforms$`GENE-devreg-s-stricto`[x] == F & isoforms$`ISOFORM-devreg-s-stricto`[x] == F) {
    isoforms$devreg_factor[x] <- "Gene NOT devreg-s-stricto / Isoform NOT devreg-s-stricto"
  } else if (isoforms$`GENE-devreg-s-stricto`[x] == F & isoforms$`ISOFORM-devreg-s-stricto`[x] == T) {
    isoforms$devreg_factor[x] <- "Gene NOT devreg-s-stricto / Isoform devreg-s-stricto"
  }
}

# takarítunk az oszlopok között, csak az marad ami kell
isoforms <- isoforms %>% select(gene_id, transcript_id, isoform_FC, devreg_factor)

##### PLOTOLÁS ######

full_FC <- left_join(isoforms, genes) %>% mutate(log2_isoform_FC = log2(isoform_FC)) %>% mutate(log2_gene_FC = log2(gene_FC))

# spec init miatt utófikszálás kell
FC_A <- full_FC %>% filter(devreg_factor == "Gene devreg-s-stricto / Isoform NOT devreg-s-stricto") %>% filter(isoform_FC > 4)
FC_B <- full_FC %>% filter(devreg_factor == "Gene NOT devreg-s-stricto / Isoform NOT devreg-s-stricto") %>% filter(isoform_FC > 4)

FC_FINAL <- setdiff(full_FC, FC_A) %>% setdiff(FC_B)

FC_FINAL_2 <- FC_FINAL %>% filter(gene_FC < 10000) %>% filter(isoform_FC < 10000)

e <- ggplot(FC_FINAL, aes(log2_isoform_FC, log2_gene_FC))
e + geom_point(aes(colour = factor(devreg_factor))) +
  theme(legend.position="right" ) +
  labs(x="Isoform's devreg-s-stricto Fold Change (log2)", y="Gene's devreg-s-stricto Fold Change (log2)")

e <- ggplot(FC_FINAL_2, aes(log2_isoform_FC, log2_gene_FC))
e + geom_point(aes(colour = factor(devreg_factor))) +
  theme(legend.position="right" ) +
  coord_fixed() + guides(col = guide_legend(nrow = 4)) +
  labs(x="Isoform's devreg-s-stricto Fold Change (log2)", y="Gene's devreg-s-stricto Fold Change (log2)")








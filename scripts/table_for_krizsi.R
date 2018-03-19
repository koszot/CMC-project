library(dplyr)
library(readr)
library(tidyr)
library(stringr)

##### aostoyae ####

setwd("~/Desktop/MTA/CMC_project/aostoyae/")

genes <- read_tsv("aostoyae_genes.tsv", col_names = T, cols(.default = col_guess(), gene_id = col_character()))
not_devreg <- genes %>% filter(`FB-devreg` == F) %>% dplyr::select(gene_id)
not_init <- genes %>% filter(`FB-init` == F) %>% dplyr::select(gene_id)
isoforms <- read_tsv("aostoyae_isoforms.tsv", col_names = T, cols(.default = col_guess(), gene_id = col_character()))

# not FB-devreg genes with FB-devreg isoforms

aostoyae_devreg <- left_join(not_devreg, isoforms) %>% filter(`GENE-isoforms` > 1) %>% filter(`ISOFORM-FB-devreg` == T) %>% select(gene_id) %>% unique()

# not FB-init genes with FB-init isoforms

aostoyae_init <- left_join(not_init, isoforms) %>% filter(`GENE-isoforms` > 1) %>% filter(`ISOFORM-FB-init` == T) %>% select(gene_id) %>% unique()

##### ccinerea ####

setwd("~/Desktop/MTA/CMC_project/ccinerea/")

genes <- read_tsv("ccinerea_genes.tsv", col_names = T, cols(.default = col_guess(), gene_id = col_character()))
not_devreg <- genes %>% filter(`FB-devreg` == F) %>% dplyr::select(gene_id)
not_init <- genes %>% filter(`FB-init` == F) %>% dplyr::select(gene_id)
isoforms <- read_tsv("ccinerea_isoforms.tsv", col_names = T, cols(.default = col_guess(), gene_id = col_character()))

# not FB-devreg genes with FB-devreg isoforms
ccinerea_devreg <- left_join(not_devreg, isoforms) %>% filter(`GENE-isoforms` > 1) %>% filter(`ISOFORM-FB-devreg` == T) %>% select(gene_id) %>% unique()
# not FB-init genes with FB-init isoforms
ccinerea_init <- left_join(not_init, isoforms) %>% filter(`GENE-isoforms` > 1) %>% filter(`ISOFORM-FB-init` == T) %>% select(gene_id) %>% unique()

# beolvassuk a headerfájlból készített szótárfájlt
fasta <- read_tsv("../../AS_project/FILES_ccinerea/ENRICHMENT_ccinerea/Copci_AmutBmut1_GeneCatalog_proteins_20130522.aa.fasta.headers", col_names = c("fasta_transcript_id", "gene_id"))
# beolvassuk az eredeti annotációt ami alapján az expressziós analyzis készült
original <- read_tsv("../../AS_project/FILES_ccinerea/GENOME_ccinerea/Copci_AmutBmut1_GeneCatalog_genes_20130522.gff", col_names = c("chr", "maker", "type", "start", "end", "att1", "strand",  "att2", "attributes"))
original <- original %>%
  filter(type == "exon") %>%
  separate(attributes, c("name", "gene_id", "tid", "annotation_transcript_id"), sep = " ") %>%
  select(gene_id, annotation_transcript_id)
original$gene_id <- original$gene_id %>% 
  str_replace("^\"", "") %>%
  str_replace("\";$", "")
original <- original[!duplicated(original), ]
final <- left_join(fasta, original, by = "gene_id")
names(final) <- c("CMC_gene_id", "gene_name", "gene_id")

ccinerea_devreg_fixed <- left_join(ccinerea_devreg, final, by = "gene_id") %>% select(CMC_gene_id)
names(ccinerea_devreg_fixed) <- "ccinerea_gene_id"
ccinerea_devreg_fixed$ccinerea_gene_id <- str_replace(ccinerea_devreg_fixed$ccinerea_gene_id, "^", "Copci_AmutBmut1_")

ccinerea_init_fixed <- left_join(ccinerea_init, final, by = "gene_id") %>% select(CMC_gene_id)
names(ccinerea_init_fixed) <- "ccinerea_gene_id"
ccinerea_init_fixed$ccinerea_gene_id <- str_replace(ccinerea_init_fixed$ccinerea_gene_id, "^", "Copci_AmutBmut1_")

##### NEMJÓ!!!! ltigrinus ####

setwd("~/Desktop/MTA/CMC_project/ltigrinus/")

genes <- read_tsv("ltigrinus_genes.tsv", col_names = T, cols(.default = col_guess(), gene_id = col_character()))
not_devreg <- genes %>% filter(`FB-devreg` == F) %>% dplyr::select(gene_id)
not_init <- genes %>% filter(`FB-init` == F) %>% dplyr::select(gene_id)
isoforms <- read_tsv("ltigrinus_isoforms.tsv", col_names = T, cols(.default = col_guess(), gene_id = col_character()))

# not FB-devreg genes with FB-devreg isoforms
ltigrinus_devreg <- left_join(not_devreg, isoforms) %>% filter(`GENE-isoforms` > 1) %>% filter(`ISOFORM-FB-devreg` == T) %>% select(gene_id) %>% unique()
# not FB-init genes with FB-init isoforms
ltigrinus_init <- left_join(not_init, isoforms) %>% filter(`GENE-isoforms` > 1) %>% filter(`ISOFORM-FB-init` == T) %>% select(gene_id) %>% unique()

# beolvassuk a headerfájlból készített szótárfájlt
fasta <- read_tsv("../../AS_project/FILES_ltigrinus/ENRICHMENT_ltigrinus/Sisbr1_GeneCatalog_proteins_20130805.aa.fasta.tsv.headers", col_names = c("fasta_transcript_id", "gene_id"))
# beolvassuk az eredeti annotációt ami alapján az expressziós analyzis készült
original <- read_tsv("../../AS_project/FILES_ltigrinus/GENOME_ltigrinus/Sisbr1_GeneCatalog_genes_20130805.gff", col_names = c("chr", "maker", "type", "start", "end", "att1", "strand",  "att2", "attributes"))
original <- original %>%
  filter(type == "exon") %>%
  separate(attributes, c("name", "gene_id", "tid", "annotation_transcript_id"), sep = " ") %>%
  select(gene_id, annotation_transcript_id)
original$gene_id <- original$gene_id %>% 
  str_replace("^\"", "") %>%
  str_replace("\";$", "")
original <- original[!duplicated(original), ]
final <- left_join(fasta, original, by = "gene_id")
names(final) <- c("CMC_gene_id", "gene_name", "gene_id")

ltigrinus_devreg_fixed <- left_join(ltigrinus_devreg, final, by = "gene_id") %>% select(CMC_gene_id)
names(ltigrinus_devreg_fixed) <- "ltigrinus_gene_id"
ltigrinus_devreg_fixed$ltigrinus_gene_id <- str_replace(ltigrinus_devreg_fixed$ltigrinus_gene_id, "^", "Lenti6_1_")

ltigrinus_init_fixed <- left_join(ltigrinus_init, final, by = "gene_id") %>% select(CMC_gene_id)
names(ltigrinus_init_fixed) <- "ltigrinus_gene_id"
ltigrinus_init_fixed$ltigrinus_gene_id <- str_replace(ltigrinus_init_fixed$ltigrinus_gene_id, "^", "Lenti6_1_")

test <- read_tsv("~/Desktop/testids")
names(test) <- "ltigrinus_gene_id"

intersect(ltigrinus_devreg_fixed, test)

##### rmellea ####

setwd("~/Desktop/MTA/CMC_project/rmellea/")

genes <- read_tsv("rmellea_genes.tsv", col_names = T, cols(.default = col_guess(), gene_id = col_character()))
not_devreg <- genes %>% filter(`FB-devreg` == F) %>% dplyr::select(gene_id)
not_init <- genes %>% filter(`FB-init` == F) %>% dplyr::select(gene_id)
isoforms <- read_tsv("rmellea_isoforms.tsv", col_names = T, cols(.default = col_guess(), gene_id = col_character()))

# not FB-devreg genes with FB-devreg isoforms
rmellea_devreg <- left_join(not_devreg, isoforms) %>% filter(`GENE-isoforms` > 1) %>% filter(`ISOFORM-FB-devreg` == T) %>% select(gene_id) %>% unique()
# not FB-init genes with FB-init isoforms
rmellea_init <- left_join(not_init, isoforms) %>% filter(`GENE-isoforms` > 1) %>% filter(`ISOFORM-FB-init` == T) %>% select(gene_id) %>% unique()

# beolvassuk a headerfájlból készített szótárfájlt
fasta <- read_tsv("../../AS_project/FILES_rmellea/ENRICHMENT_rmellea/Ricmel1_GeneCatalog_proteins_20151108.aa.fasta.headers", col_names = c("fasta_transcript_id", "gene_id"))
# beolvassuk az eredeti annotációt ami alapján az expressziós analyzis készült
original <- read_tsv("../../AS_project/FILES_rmellea/GENOME_rmellea/Ricmel1_GeneCatalog_genes_20151108.gff", col_names = c("chr", "maker", "type", "start", "end", "att1", "strand",  "att2", "attributes"))
original <- original %>%
  filter(type == "exon") %>%
  separate(attributes, c("name", "gene_id", "tid", "annotation_transcript_id"), sep = " ") %>%
  select(gene_id, annotation_transcript_id)
original$gene_id <- original$gene_id %>% 
  str_replace("^\"", "") %>%
  str_replace("\";$", "")
original <- original[!duplicated(original), ]
final <- left_join(fasta, original, by = "gene_id")
names(final) <- c("CMC_gene_id", "gene_name", "gene_id")

rmellea_devreg_fixed <- left_join(rmellea_devreg, final, by = "gene_id") %>% select(CMC_gene_id)
names(rmellea_devreg_fixed) <- "rmellea_gene_id"
rmellea_devreg_fixed$rmellea_gene_id <- str_replace(rmellea_devreg_fixed$rmellea_gene_id, "^", "Ricmel1_")

rmellea_init_fixed <- left_join(rmellea_init, final, by = "gene_id") %>% select(CMC_gene_id)
names(rmellea_init_fixed) <- "rmellea_gene_id"
rmellea_init_fixed$rmellea_gene_id <- str_replace(rmellea_init_fixed$rmellea_gene_id, "^", "Ricmel1_")

test <- read_tsv("~/Desktop/testids")
names(test) <- "rmellea_gene_id"

intersect(rmellea_init_fixed, test)

##### scommune ####

setwd("~/Desktop/MTA/CMC_project/scommune/")

genes <- read_tsv("scommune_genes.tsv", col_names = T, cols(.default = col_guess(), gene_id = col_character()))
not_devreg <- genes %>% filter(`FB-devreg` == F) %>% dplyr::select(gene_id)
not_init <- genes %>% filter(`FB-init` == F) %>% dplyr::select(gene_id)
isoforms <- read_tsv("scommune_isoforms.tsv", col_names = T, cols(.default = col_guess(), gene_id = col_character()))

# not FB-devreg genes with FB-devreg isoforms
scommune_devreg <- left_join(not_devreg, isoforms) %>% filter(`GENE-isoforms` > 1) %>% filter(`ISOFORM-FB-devreg` == T) %>% select(gene_id) %>% unique()
# not FB-init genes with FB-init isoforms
scommune_init <- left_join(not_init, isoforms) %>% filter(`GENE-isoforms` > 1) %>% filter(`ISOFORM-FB-init` == T) %>% select(gene_id) %>% unique()

# beolvassuk a headerfájlból készített szótárfájlt
fasta <- read_tsv("../../AS_project/FILES_scommune/ENRICHMENT_scommune/Schco3_GeneCatalog_proteins_20130812.aa.fasta.headers", col_names = c("fasta_transcript_id", "gene_id"))
# beolvassuk az eredeti annotációt ami alapján az expressziós analyzis készült
original <- read_tsv("../../AS_project/FILES_scommune/GENOME_scommune/Schco3_GeneCatalog_genes_20130812.gff", col_names = c("chr", "maker", "type", "start", "end", "att1", "strand",  "att2", "attributes"))
original <- original %>%
  filter(type == "exon") %>%
  separate(attributes, c("name", "gene_id", "tid", "annotation_transcript_id"), sep = " ") %>%
  select(gene_id, annotation_transcript_id)
original$gene_id <- original$gene_id %>% 
  str_replace("^\"", "") %>%
  str_replace("\";$", "")
original <- original[!duplicated(original), ]
final <- left_join(fasta, original, by = "gene_id")
names(final) <- c("CMC_gene_id", "gene_name", "gene_id")

scommune_devreg_fixed <- left_join(scommune_devreg, final, by = "gene_id") %>% select(CMC_gene_id)
names(scommune_devreg_fixed) <- "scommune_gene_id"
scommune_devreg_fixed$scommune_gene_id <- str_replace(scommune_devreg_fixed$scommune_gene_id, "^", "Schco3_")

scommune_init_fixed <- left_join(scommune_init, final, by = "gene_id") %>% select(CMC_gene_id)
names(scommune_init_fixed) <- "scommune_gene_id"
scommune_init_fixed$scommune_gene_id <- str_replace(scommune_init_fixed$scommune_gene_id, "^", "Schco3_")

test <- read_tsv("~/Desktop/testids")
names(test) <- "scommune_gene_id"

intersect(scommune_devreg_fixed, test)

##### pchrysosporium ####

setwd("~/Desktop/MTA/CMC_project/pchrysosporium/")

genes <- read_tsv("pchrysosporium_genes.tsv", col_names = T, cols(.default = col_guess(), gene_id = col_character()))
not_devreg <- genes %>% filter(`FB-devreg` == F) %>% dplyr::select(gene_id)
not_init <- genes %>% filter(`FB-init` == F) %>% dplyr::select(gene_id)
isoforms <- read_tsv("pchrysosporium_isoforms.tsv", col_names = T, cols(.default = col_guess(), gene_id = col_character()))

# not FB-devreg genes with FB-devreg isoforms
pchrysosporium_devreg <- left_join(not_devreg, isoforms) %>% filter(`GENE-isoforms` > 1) %>% filter(`ISOFORM-FB-devreg` == T) %>% select(gene_id) %>% unique()
# not FB-init genes with FB-init isoforms
pchrysosporium_init <- left_join(not_init, isoforms) %>% filter(`GENE-isoforms` > 1) %>% filter(`ISOFORM-FB-init` == T) %>% select(gene_id) %>% unique()

# beolvassuk a headerfájlból készített szótárfájlt
fasta <- read_tsv("../../AS_project/FILES_pchrysosporium/ENRICHMENT_pchrysosporium/Phchr2_GeneCatalog_proteins_20131210.aa.fasta.headers", col_names = c("fasta_transcript_id", "gene_id"))
# beolvassuk az eredeti annotációt ami alapján az expressziós analyzis készült
original <- read_tsv("../../AS_project/FILES_pchrysosporium/GENOME_pchrysosporium/Phchr2_GeneCatalog_genes_20131210.gff", col_names = c("chr", "maker", "type", "start", "end", "att1", "strand",  "att2", "attributes"))
original <- original %>%
  filter(type == "exon") %>%
  separate(attributes, c("name", "gene_id", "tid", "annotation_transcript_id"), sep = " ") %>%
  select(gene_id, annotation_transcript_id)
original$gene_id <- original$gene_id %>% 
  str_replace("^\"", "") %>%
  str_replace("\";$", "")
original <- original[!duplicated(original), ]
final <- left_join(fasta, original, by = "gene_id")
names(final) <- c("CMC_gene_id", "gene_name", "gene_id")

pchrysosporium_devreg_fixed <- left_join(pchrysosporium_devreg, final, by = "gene_id") %>% select(CMC_gene_id)
names(pchrysosporium_devreg_fixed) <- "pchrysosporium_gene_id"
pchrysosporium_devreg_fixed$pchrysosporium_gene_id <- str_replace(pchrysosporium_devreg_fixed$pchrysosporium_gene_id, "^", "Phchr2_")

pchrysosporium_init_fixed <- left_join(pchrysosporium_init, final, by = "gene_id") %>% select(CMC_gene_id)
names(pchrysosporium_init_fixed) <- "pchrysosporium_gene_id"
pchrysosporium_init_fixed$pchrysosporium_gene_id <- str_replace(pchrysosporium_init_fixed$pchrysosporium_gene_id, "^", "Phchr2_")

test <- read_tsv("~/Desktop/testids")
names(test) <- "pchrysosporium_gene_id"

intersect(pchrysosporium_init_fixed, test)

##### TABLE kiírása #####

setwd("~/Desktop/MTA/CMC_project/not_devreg_gene_devreg_isoform/")
write_tsv(aostoyae_devreg, "aostoyae_devreg.tsv")
write_tsv(ccinerea_devreg_fixed, "ccinerea_devreg.tsv")
write_tsv(rmellea_devreg_fixed, "rmellea_devreg.tsv")
write_tsv(scommune_devreg_fixed, "scommune_devreg.tsv")
write_tsv(pchrysosporium_devreg_fixed, "pchrysosprium_devreg.tsv")

setwd("~/Desktop/MTA/CMC_project/not_init_gene_init_isoform/")
write_tsv(aostoyae_init, "aostoyae_init.tsv")
write_tsv(ccinerea_init_fixed, "ccinerea_init.tsv")
write_tsv(rmellea_init_fixed, "rmellea_init.tsv")
write_tsv(scommune_init_fixed, "scommune_init.tsv")
write_tsv(pchrysosporium_init_fixed, "pchrysosprium_init.tsv")


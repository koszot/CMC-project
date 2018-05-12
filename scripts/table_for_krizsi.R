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



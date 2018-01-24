library(dplyr)
library(readr)
library(seqinr)
library(stringr)
library(tidyr)

FASTA_file <- read.fasta(file = "~/Desktop/75genomes_plus_DMel.input.fasta",  seqtype = "AA", set.attributes = FALSE)


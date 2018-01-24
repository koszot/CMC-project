library(readr)
library(dplyr)
library(ape)
library(stringr)

dict <- read_tsv("~/Desktop/eniko/compare_reference", col_names = F)
species <- read.tree("~/Desktop/eniko/trees_for_Notung_01_09_rem_MS_result/species.tree.nwk") # behívjuk a labelezett species tree fájlt

for (i in 1:length(dict$X1)) {
  x <- as.character(dict[i,3])
  x <- unlist(str_split(x, " "))
  if (length(x) == 1) {
    dict[i,1] <- dict[i,3]
  } else {
    MRCAnode <- getMRCA(species,x)
    realnode <- MRCAnode - length(species$tip.label)
    node <- species[["node.label"]][realnode]
    dict[i,1] <- node
  }
}
  
write_tsv(dict, "~/Desktop/eniko/notung_compare_dict_FINAL.tsv")

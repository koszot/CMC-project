library(readr)
library(dplyr)
library(stringr)

dict <- read_tsv("Desktop/eniko/notung_compare_dict_FINAL.tsv")
species <- read_lines("Desktop/eniko/species.tree.nwk") # behívjuk a labelezett species tree fájlt

species_temp <- species

for (i in 1:141) {
  to_change <- as.character(dict[i,1])
  change_for <- as.character(dict[i,2])
  species_temp <- str_replace(species_temp, to_change, change_for)
}

write_lines(species_temp, "Desktop/eniko/species_FIXED.tree")
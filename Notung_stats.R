library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(tibble)

##### LOAD RAW DATA #####

setwd("~/Desktop/SMC_project/Notung/Notung_results/")          # megadjuk a working directorynak ahol a NOTUNG statok vannak
temp <- list.files(pattern = "*.txt")            # csinálunk egy listát ami tartalmazza a fájlneveket
tree_stats <- lapply(temp, read_lines)           # behívjuk az összes fájlt a mappából egy listába

tree_stats_names <- temp %>%
  str_replace("_boot.tree_mod.tree.rearrange.0.parsable.txt", "") %>%      # átalakítjuk a temp lista elemeit, hogy csak a Cluster[0-9] maradjon meg
  str_replace("RAxML_bipartitions_", "")
names(tree_stats) <- tree_stats_names                       # átnevezzük a lista elemeit annak megfelelően hogy melyik clusterből származik a NOTUNG output

dict <- read_tsv("~/Desktop/SMC_project/Notung/notung_compare_dict_FINAL.tsv")      # betöltjük a COMPARE/NOTUNG szótárat
names(dict) <- c("notung", "compare", "MRCA")

##### PROCESS RAW DATA #####

for (i in 1:length(tree_stats)) {                                                 # for ciklussal végiglépkedünk a tree_stats-on
  parsable <- tree_stats[[i]]                                                     # kiemeljük a cluster statisztikáját
  parsable_formated <- parsable %>%
    str_extract("^#S.*") %>%                                                      # kiemeljük csak a minket érdeklő sorokat
    str_replace("#S\t", "") %>%                                                   # kitöröljük a nem táblázatosítható előtagokat
    str_replace("Species\t.*", "header") %>%                                      # a header sort átnevezzük header-re a könnyebb törlés végett
    as_tibble() %>%                                                               # átalakítjuk tibble-re
    filter(!is.na(value)) %>%                                                     # kitöröljük az üres sorokat
    filter(value != "header") %>%                                                 # kitöröljük a headert
    separate(value, c("NOTUNG_like", "Duplications",  "cds", "Trans_from", "Trans_to", "Losses"), sep = "\t") %>%   # szétszedjük oszlpokra az adatsort
    select(NOTUNG_like, Duplications, Losses)                                     # kiemeljük a minket érdeklő oszlopokat
  for (f in 1:length(parsable_formated$NOTUNG_like)) {                            # új for ciklus végig a szótár hosszán
    dict_temp <- dict[dict$notung %in% as.character(parsable_formated[f,1]),]     # kimentjük az aktuális szótárelemet
    parsable_formated[f,4] <- dict_temp[1,2]                                      # belementjük az új oszlopba a NOTUNG-nak megfelelő COMPARE node ID-t
  }
  parsable_final <- select(parsable_formated, NOTUNG_like, COMPARE_like = compare, Duplications, Losses)  # átalakítjuk a táblázatot a végső formátumra
  tree_stats[[i]] <- parsable_final                                               # belementjük a táblázatot a tree_stats lista azon elemébe ahonnan a raw datát kiemeltük
  pb <- txtProgressBar(min = 1, max = length(tree_stats), style = 3)              # progress bar
  setTxtProgressBar(pb, i, title = NULL, label = NULL)
}

rm(dict_temp,parsable,temp,f,i,parsable_final,parsable_formated,tree_stats_names, pb)     # kitöröljük a felesleges változókat

###### LOSSES ######

losses_details <- vector("list", length(dict$compare))    # létrehozunk egy listát amibe a cluster/loss-okat mentjük
node_names <- dict$compare
names(losses_details) <- node_names                       # elnevezzük a lista elemeit a COMPARE node ID-nak megfelelően
losses_clustercount <- losses_details                     # létrehozunk egy listát amibe az a klaszterek száma lesz amikben loss van nodeonként
losses_count <- losses_details                            # létrehozunk egy listát amibe az a lossok száma lesz node-onként
 
rm(node_names,dict)                                       # kitöröljük a felesleges változókat

# az adott node-ra vonatkozó lossok összeszámolása

for (i in 1:length(losses_count)) {                                           # for ciklussal végiglépkedünk a losses_count-on
  temp_node_name <- names(losses_count[i])                                    # kiszedjük az adott listaelem nevét ami a node ID-t tartalmazza
  losses_count_for_node <- 0                                                  # létrehozunk egy integert amibe az adott cluster-hez tartozó összlosst mentjük
  for (f in 1:length(tree_stats)) {                                           # for ciklussal végiglépkedünk a tree_stats-on
    temp_losses_count <- tree_stats[[f]] %>%                                  # kiemeljük az aktuális clusterhez tartozó adatokat
      filter(COMPARE_like == temp_node_name) %>%                              # tovább filterezzük az aktuális node alapján
      select(Losses)                                                          # végül kiemeljük az adott node-ra vonatkozó loss-számot
    temp_losses_count <- as.integer(temp_losses_count)                        # átalakítjuk integerre
    losses_count_for_node <- losses_count_for_node + temp_losses_count        # összeadjuk a különböző clusterek adatait az adott node-ra
  }
  losses_count[[i]] <- losses_count_for_node                                  # bementjük a megfelelő node-hoz az összeget
  rm(losses_count_for_node)                                                   # takarítás a következő loop előtt (nem szükséges)
  pb <- txtProgressBar(min = 1, max = length(losses_count), style = 3)        # progress bar
  setTxtProgressBar(pb, i, title = NULL, label = NULL)
}

tbl_losses_count <- as_tibble(unlist(losses_count))                       # tibble formátumra alakítjuk
tbl_losses_count <- rownames_to_column(tbl_losses_count, var = "node")    # a sorok neveit (node_ID) átemeljük oszlopnak
colnames(tbl_losses_count) <- c("node_ID", "all_losses")                  # újranevezzük az oszlopokat

rm(losses_count,f,i,pb,temp_losses_count,temp_node_name)                  # takarítás

# az adott node-ra vonatkozó lossok kigyűjtése clusterenként --> clusterID/loss_number

for (i in 1:length(losses_details)) {                                          # for ciklussal végiglépkedünk a losses_details-en
  temp_node_name <- names(losses_details[i])                                   # kiszedjük az adott listaelem nevét ami a node ID-t tartalmazza
  losses_details_for_node <- as.character(vector())                            # létrehozunk egy karaktervektort amibe az adott cluster-hez tartozó lossokat mentjük
  for (f in 1:length(tree_stats)) {                                            # for ciklussal végiglépkedünk a tree_stats-on
    temp_loss <- tree_stats[[f]] %>%                                           # kiemeljük az aktuális clusterhez tartozó adatokat
      filter(COMPARE_like == temp_node_name) %>%                               # tovább filterezzük az aktuális node alapján
      select(Losses)                                                           # végül kiemeljük az adott node-ra vonatkozó loss-számot
    temp_clusterloss <- c(names(tree_stats[f]), "/", temp_loss)                # létrehozzuk a clusterID/loss_number szisztémát
    temp_clusterloss <- as.character(temp_clusterloss)
    temp_clusterloss <- str_c(temp_clusterloss, collapse = "")
    losses_details_for_node <- c(losses_details_for_node, temp_clusterloss)    # összegyüjtjük egymás mellé a node-hoz tartozó lossokat clusterenként
  }
  losses_details[[i]] <- losses_details_for_node                               # bementjük a megfelelő node-hoz az összeget
  rm(losses_details_for_node)                                                  # takarítás a következő loop előtt 
  pb <- txtProgressBar(min = 1, max = length(losses_details), style = 3)       # progress bar
  setTxtProgressBar(pb, i, title = NULL, label = NULL)
}

for (i in 1:length(losses_details)) {                                           # for ciklussal végiglépkedünk a losses_details-en
  losses_details[[i]] <- str_replace(losses_details[[i]], ".*/0", "empty")      # kitörölgetjük azokat a clustereket amikhez nincs érték
  losses_details[[i]] <- losses_details[[i]][losses_details[[i]] != "empty"]   
  losses_details[[i]] <- str_c(losses_details[[i]], collapse = ", ")            # összeolvasztjuk az eredményt egy string-be
  pb <- txtProgressBar(min = 1, max = length(losses_details), style = 3)        # progress bar
  setTxtProgressBar(pb, i, title = NULL, label = NULL)
}

tbl_losses_details <- as_tibble(unlist(losses_details))                          # tibble formátumra alakítjuk
tbl_losses_details <- rownames_to_column(tbl_losses_details, var = "node")       # a sorok neveit (node_ID) átemeljük oszlopnak
colnames(tbl_losses_details) <- c("node_ID", "losses_details")                   # újranevezzük az oszlopokat

rm(temp_loss,f,i,pb,temp_clusterloss,temp_node_name)                             # takarítás

# az adott node-ra vonatkozó clusterek összeszámolása amik losst tartalmaznak

for (i in 1:length(losses_clustercount)) {
  temp_cc <- losses_details[[i]]
  losses_clustercount[[i]] <- length(unlist(str_split(temp_cc, ", ")))
  pb <- txtProgressBar(min = 1, max = length(losses_details), style = 3)        # progress bar
  setTxtProgressBar(pb, i, title = NULL, label = NULL)
}

tbl_losses_clustercount <- as_tibble(unlist(losses_clustercount))                        # tibble formátumra alakítjuk
tbl_losses_clustercount <- rownames_to_column(tbl_losses_clustercount, var = "node")     # a sorok neveit (node_ID) átemeljük oszlopnak
colnames(tbl_losses_clustercount) <- c("node_ID", "all_clusters")                        # újranevezzük az oszlopokat

rm(i,losses_clustercount,losses_details,pb,temp_cc)            # takarítás

tbl_losses <- tbl_losses_count %>%                             # végső táblázat összeállítása
  left_join(tbl_losses_clustercount, by = "node_ID") %>%
  left_join(tbl_losses_details, by = "node_ID")
tbl_losses$event <- "DELETION"
tbl_losses <- select(tbl_losses, event, node_ID, event_count = all_losses,  cluster_count = all_clusters, details = losses_details)

rm(tbl_losses_count,tbl_losses_clustercount,tbl_losses_details)            # takarítás

###### DUPLICATION #####

dict <- read_tsv("~/Desktop/SMC_project/Notung/notung_compare_dict_FINAL.tsv")      # betöltjük a COMPARE/NOTUNG szótárat
names(dict) <- c("notung", "compare", "MRCA")

duplication_details <- vector("list", length(dict$compare))   # létrehozunk egy listát amibe a cluster/duplication-okat mentjük
node_names <- dict$compare
names(duplication_details) <- node_names                      # elnevezzük a lista elemeit a COMPARE node ID-nak megfelelően
duplication_clustercount <- duplication_details                    # létrehozunk egy listát amibe az a klaszterek száma lesz amikben duplication van nodeonként
duplication_count <- duplication_details                           # létrehozunk egy listát amibe az a duplication-ok száma lesz node-onként

rm(node_names,dict)                                           # kitöröljük a felesleges változókat

# az adott node-ra vonatkozó duplication-ok összeszámolása

for (i in 1:length(duplication_count)) {                                                # for ciklussal végiglépkedünk a duplication_count-on
  temp_node_name <- names(duplication_count[i])                                         # kiszedjük az adott listaelem nevét ami a node ID-t tartalmazza
  duplication_count_for_node <- 0                                                       # létrehozunk egy integert amibe az adott cluster-hez tartozó összduplicationt mentjük
  for (f in 1:length(tree_stats)) {                                                     # for ciklussal végiglépkedünk a tree_stats-on
    temp_duplication_count <- tree_stats[[f]] %>%                                       # kiemeljük az aktuális clusterhez tartozó adatokat
      filter(COMPARE_like == temp_node_name) %>%                                        # tovább filterezzük az aktuális node alapján
      select(Duplications)                                                              # végül kiemeljük az adott node-ra vonatkozó loss-számot
    temp_duplication_count <- as.integer(temp_duplication_count)                        # átalakítjuk integerre
    duplication_count_for_node <- duplication_count_for_node + temp_duplication_count   # összeadjuk a különböző clusterek adatait az adott node-ra
  }
  duplication_count[[i]] <- duplication_count_for_node                                  # bementjük a megfelelő node-hoz az összeget
  rm(duplication_count_for_node)                                                        # takarítás a következő loop előtt (nem szükséges)
  pb <- txtProgressBar(min = 1, max = length(duplication_count), style = 3)             # progress bar
  setTxtProgressBar(pb, i, title = NULL, label = NULL)
}

tbl_duplication_count <- as_tibble(unlist(duplication_count))                       # tibble formátumra alakítjuk
tbl_duplication_count <- rownames_to_column(tbl_duplication_count, var = "node")    # a sorok neveit (node_ID) átemeljük oszlopnak
colnames(tbl_duplication_count) <- c("node_ID", "all_duplication")                  # újranevezzük az oszlopokat

rm(duplication_count,f,i,pb,temp_duplication_count,temp_node_name)                  # takarítás

# az adott node-ra vonatkozó duplicationok kigyűjtése clusterenként --> clusterID/duplication_number

for (i in 1:length(duplication_details)) {                                          # for ciklussal végiglépkedünk a duplication_details-en
  temp_node_name <- names(duplication_details[i])                                   # kiszedjük az adott listaelem nevét ami a node ID-t tartalmazza
  duplication_details_for_node <- as.character(vector())                            # létrehozunk egy karaktervektort amibe az adott cluster-hez tartozó duplicationokat mentjük
  for (f in 1:length(tree_stats)) {                                                 # for ciklussal végiglépkedünk a tree_stats-on
    temp_dups <- tree_stats[[f]] %>%                                                # kiemeljük az aktuális clusterhez tartozó adatokat
      filter(COMPARE_like == temp_node_name) %>%                                    # tovább filterezzük az aktuális node alapján
      select(Duplications)                                                          # végül kiemeljük az adott node-ra vonatkozó duplication-számot
    temp_clusterdups <- c(names(tree_stats[f]), "/", temp_dups)                     # létrehozzuk a clusterID/duplication_number szisztémát
    temp_clusterdups <- as.character(temp_clusterdups)
    temp_clusterdups <- str_c(temp_clusterdups, collapse = "")
    duplication_details_for_node <- c(duplication_details_for_node, temp_clusterdups)    # összegyüjtjük egymás mellé a node-hoz tartozó duplicationokat clusterenként
  }
  duplication_details[[i]] <- duplication_details_for_node                          # bementjük a megfelelő node-hoz az összeget
  rm(duplication_details_for_node)                                                  # takarítás a következő loop előtt 
  pb <- txtProgressBar(min = 1, max = length(duplication_details), style = 3)       # progress bar
  setTxtProgressBar(pb, i, title = NULL, label = NULL)
}

for (i in 1:length(duplication_details)) {                                                   # for ciklussal végiglépkedünk a duplication_details-en
  duplication_details[[i]] <- str_replace(duplication_details[[i]], ".*/0", "empty")         # kitörölgetjük azokat a clustereket amikhez nincs érték
  duplication_details[[i]] <- duplication_details[[i]][duplication_details[[i]] != "empty"]   
  duplication_details[[i]] <- str_c(duplication_details[[i]], collapse = ", ")               # összeolvasztjuk az eredményt egy string-be
  pb <- txtProgressBar(min = 1, max = length(duplication_details), style = 3)                # progress bar
  setTxtProgressBar(pb, i, title = NULL, label = NULL)
}

tbl_duplication_details <- as_tibble(unlist(duplication_details))                          # tibble formátumra alakítjuk
tbl_duplication_details <- rownames_to_column(tbl_duplication_details, var = "node")       # a sorok neveit (node_ID) átemeljük oszlopnak
colnames(tbl_duplication_details) <- c("node_ID", "duplication_details")                   # újranevezzük az oszlopokat

rm(temp_clusterdups,f,i,pb,temp_node_name,temp_dups)                             # takarítás

# az adott node-ra vonatkozó clusterek összeszámolása amik duplicationt tartalmaznak

for (i in 1:length(duplication_clustercount)) {
  temp_cc <- duplication_details[[i]]
  duplication_clustercount[[i]] <- length(unlist(str_split(temp_cc, ", ")))
  pb <- txtProgressBar(min = 1, max = length(duplication_details), style = 3)        # progress bar
  setTxtProgressBar(pb, i, title = NULL, label = NULL)
}

tbl_duplication_clustercount <- as_tibble(unlist(duplication_clustercount))                        # tibble formátumra alakítjuk
tbl_duplication_clustercount <- rownames_to_column(tbl_duplication_clustercount, var = "node")     # a sorok neveit (node_ID) átemeljük oszlopnak
colnames(tbl_duplication_clustercount) <- c("node_ID", "all_clusters")                             # újranevezzük az oszlopokat

rm(i,duplication_clustercount,duplication_details,pb,temp_cc)            # takarítás

tbl_duplication <- tbl_duplication_count %>%                             # végső táblázat összeállítása
  left_join(tbl_duplication_clustercount, by = "node_ID") %>%
  left_join(tbl_duplication_details, by = "node_ID")
tbl_duplication$event <- "DUPLICATION"
tbl_duplication <- select(tbl_duplication, event, node_ID, event_count = all_duplication,  cluster_count = all_clusters, details = duplication_details)

rm(tbl_duplication_count,tbl_duplication_clustercount,tbl_duplication_details)            # takarítás

final <- bind_rows(tbl_duplication,tbl_losses)

##### WRITE FILE #####

write_tsv(final, "~/Desktop/SMC_project/Notung/NOTUNG_stats.tsv")



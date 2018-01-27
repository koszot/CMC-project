# Complex Multicellularity Project

# Armillaria ostoyae

## Táblázat elkészítése

### Script:
- TABLE_aostoyae.R

### Input:
- __isoforms.fpkm_tracking__
  - CuffDiff output FPKM érték táblázat az izoformákra
- __genes.fpkm_tracking__
  - CuffDiff output FPKM érték táblázat az génekre
- __aostoyae_AS_annotation.gtf__
  - a végső annotációs fájl amivel a CuffDiff analízis is futott

### Output:
- __aostoyae_isoforms.tsv__
  - __VM_FPKM_Rank__ -> __FB_S_FPKM_Rank__ __oszlopok__: megadja, hogy az adott izoforma expressziós értéke hanyadik helyre van sorolva az összes izoforma expressziós értéke közül az adott fejlődési fázisban/szövettípusban, megmutatja hogy az adott izoforma mennyire domináns az adott fejlődési fázisban/szövettípusban
  - __Rank_Sum__ __oszlop__: az izoformák dominancia értékei/rank-jei összeadva az összes fejlődési fázisban/szövettípusban
  - __Overall_FPKM__ __oszlok__: az izoformák teljes expressziós értéke az összes fejlődési fázisban/szövettípusban
  - __isoforms__ : a gén izoformáinak mennyisége
  - __isoform__ __switch__ : a rank-ok különböző fejlődési fázisok/szövettípusok közötti váltakozása miatt, megfigyelhető hogy egyes izoformák dominánsak-e egy fázisban majd a dominanciájukat elveszítik-e egy másikban
  - __primary__ __transcript__ : azt a transzkriptet nevezzük primary transzkriptnek amelyik összeadott dominancia értéke (Rank_Sum) a legkisebb, vagyis a legdominánsabb az izoformák közül, amennyiben ezt az értéket több izoforma is felveszi, az izoforma összeexpressziójának értéke (Overall_FPKM) alapján döntünk, vagyis az lesz a primary transcript amelyiknek a legnagyobb az expressziós értéke az azonos rank-ú izoformák közül
- __aostoyae_genes.tsv__
  - __FB-devreg__
    - legalább 4 FPKM értékkel rendelkezik P - FB között valamelyik fejlődési fázisban/szövettípusban
    - legalább 4x Fold Change-el rendelkezik valamelyik fejlődési útvonalon
      - DEVREG_S : P1 -> P2_S -> YFB_S -> FB_S
      - DEVREG_C : P1 -> P2_C -> YFB_C -> FB_C
      - DEVREG_L : P1 -> P2_C -> YFB_C -> FB_L
      - DEVREG_P2 : P2_S -> P2_C
      - DEVREG_YFB : YFB_S -> YFB_C
      - DEVREG_FB : FB_S -> FB_C -> FB_L
  - __FB-init__
    - legalább 4 FPKM értékkel rendelkezik VM - P között valamelyik fejlődési fázisban/szövettípusban
    - legalább 4x pozitív Fold Change-el rendelkezik VM - P átmenetben
- __aostoyae_corrected_annotation.gtf__
  - újrarendezzük az annotációs fájl sorait az izoforma dominancia alapján, hogy a primary transcript legyen elől, az ASpli ezután a primary transcripthez viszonyítva adja meg az alternative splicing eseményeket

## Statisztikák kinyerése

### Script:
- STATS_aostoyae.R

### Input:
- __aostoyae_genes.tsv__
- __aostoyae_isoforms.tsv__
- __aostoyae_corrected_annotation.gtf__

### Output:
- __aostoyae_stats.tsv__
  - a statisztikai adatokat tartalmazó fájl, részletezése az utolsó lépésnél
- __aostoyae_AS_genes.gtf__
  - alternative splicingban résztvevő géneket tartalmazó annotációs fájl az ASpli számára
- __aostoyae_FB_DEVREG_genes.gtf__
  - alternative splicingban résztvevő FB-devreg géneket tartalmazó annotációs fájl az ASpli számára
- __aostoyae_FB_INIT_genes.gtf__
  - alternative splicingban résztvevő FB-init géneket tartalmazó annotációs fájl az ASpli számára

## Alternative Splicing statisztikák kinyerése

### Script:
- AS_aostoyae.R

### Input:
- __aostoyae_stats.tsv__
- __aostoyae_AS_genes.gtf__
- __aostoyae_FB_DEVREG_genes.gtf__
- __aostoyae_FB_INIT_genes.gtf__

### Output:
- __aostoyae_stats.tsv__
  - __genes__ : összes annotált gén
  - __FB-devreg__ __genes__ : FB-devreg gének száma | FB-devreg gének százaléka az összes annotált génhez viszonyítva
  - __FB-init__ __genes__ : FB-init gének száma | FB-init gének százaléka az összes annotált génhez viszonyítva
  - __genes__ __with__ __AS__ : alternative splicingban résztvevő gének száma | alternative splicingban résztvevő gének százaléka az összes annotált génhez viszonyítva
  - __FB-devreg__ __genes__ __with__ __AS__ : alternative splicingban résztvevő FB-devreg gének száma | alternative splicingban résztvevő FB-devreg gének százaléka az összes annotált génhez viszonyítva
  - __FB-init__ __genes__ __with__ __AS__ : alternative splicingban résztvevő FB-init gének száma | alternative splicingban résztvevő FB-init gének százaléka az összes annotált génhez viszonyítva
  - __genes__ __with__ __AS__ __ES/IR/ALT5SS/ALT3SS__ : alternative splicingban résztvevő gének AS eseményeinek száma | alternative splicingban résztvevő gének AS eseményeinek százalékos aránya
  - __FB-devreg__ __genes__ __with__ __AS__ __ES/IR/ALT5SS/ALT3SS__ : alternative splicingban résztvevő FB-devreg gének AS eseményeinek száma | alternative splicingban résztvevő FB-devreg gének AS eseményeinek százalékos aránya
  - __FB-init__ __genes__ __with__ __AS__ __ES/IR/ALT5SS/ALT3SS__ : alternative splicingban résztvevő FB-init gének AS eseményeinek száma | alternative splicingban résztvevő FB-init gének AS eseményeinek százalékos aránya


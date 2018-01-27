# Complex Multicellularity Project

# Armillaria ostoyae

## Táblázat elkészítése

### Input:
- isoforms.fpkm_tracking
  - CuffDiff output FPKM érték táblázat az izoformákra
- genes.fpkm_tracking
  - CuffDiff output FPKM érték táblázat az génekre
- aostoyae_AS_annotation.gtf
  - a végső annotációs fájl amivel a CuffDiff analízis is futott

### Output:
- aostoyae_isoforms.tsv
  - __VM_FPKM_Rank__ -> __FB_S_FPKM_Rank__ __oszlopok__: megadja, hogy az adott izoforma expressziós értéke hanyadik helyre van sorolva az összes izoforma expressziós értéke közül az adott fejlődési fázisban/szövettípusban, megmutatja hogy az adott izoforma mennyire domináns az adott fejlődési fázisban/szövettípusban
  - __Rank_Sum__ __oszlop__: az izoformák dominancia értékei/rank-jei összeadva az összes fejlődési fázisban/szövettípusban
  - __Overall_FPKM__ __oszlok__: az izoformák teljes expressziós értéke az összes fejlődési fázisban/szövettípusban
  - __isoform__ __switch__ : a rank-ok különböző fejlődési fázisok/szövettípusok közötti váltakozása miatt, megfigyelhető hogy egyes izoformák dominánsak-e egy fázisban majd a dominanciájukat elveszítik-e egy másikban
  - __primary__ __transcript__ : azt a transzkriptet nevezzük primary transzkriptnek amelyik összeadott dominancia értéke (Rank_Sum) a legkisebb, vagyis a legdominánsabb az izoformák közül, amennyiben ezt az értéket több izoforma is felveszi, az izoforma összeexpressziójának értéke (Overall_FPKM) alapján döntünk, vagyis az lesz a primary transcript amelyiknek a legnagyobb az expressziós értéke az azonos rank-ú izoformák közül
- aostoyae_genes.tsv
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

  
- aostoyae_corrected_annotation.gtf

### Script:
- STATS_TABLE_aostoyae.R


# Complex Multicellularity Project

# Armillaria ostoyae

## Táblázat elkészítése

### Script:
- TABLE_aostoyae.R

### Input:
- __isoforms.fpkm_tracking__ : CuffDiff output FPKM érték táblázat az izoformákra
- __genes.fpkm_tracking__ : CuffDiff output FPKM érték táblázat az génekre
- __aostoyae_AS_annotation.gtf__ : a végső annotációs fájl amivel a CuffDiff analízis is futott

### Output: aostoyae_genes.tsv
- __gene_id__ (character) : az adott gén azonosítója
- __VM/P1/P2_S/P2_C/YFB_S/YFB_C/FB_S/FB_L/FB_C_FPKM__ (integer) : az adott génhez tartozó expressziós érték minden fejlődési fázisban/szövettípusban
- __FB-devreg__ (boolean)
  - legalább 4 FPKM értékkel rendelkezik P1 - FB között valamelyik fejlődési fázisban/szövettípusban
  - legalább 4x Fold Change-el rendelkezik valamelyik fejlődési útvonalon
    - DEVREG_S : P1 -> P2_S -> YFB_S -> FB_S
    - DEVREG_C : P1 -> P2_C -> YFB_C -> FB_C
    - DEVREG_L : P1 -> P2_C -> YFB_C -> FB_L
    - DEVREG_P2 : P2_S -> P2_C
    - DEVREG_YFB : YFB_S -> YFB_C
    - DEVREG_FB : FB_S -> FB_C -> FB_L
- __FB-init__ (boolean)
  - legalább 4 FPKM értékkel rendelkezik VM - P1 között valamelyik fejlődési fázisban/szövettípusban
  - legalább 4x pozitív Fold Change-el rendelkezik VM - P1 átmenetben
- __devreg-s-stricto__ (boolean) : vagy az FB-devreg, vagy az FB-init TRUE-t vesz fel
- __isoforms__ (integer) : a gén izoformáinak mennyisége

### Output: aostoyae_isoforms.tsv
- __gene_id__ (character) : az adott gén azonosítója
- __transcript_id__ (character) : az adott génhez tartozó izoformák azonosítója
- __VM/P1/P2_S/P2_C/YFB_S/YFB_C/FB_S/FB_L/FB_C_FPKM__ (integer) : az adott izoformához tartozó expressziós érték minden fejlődési fázisban/szövettípusban
- __VM/P1/P2_S/P2_C/YFB_S/YFB_C/FB_S/FB_L/FB_C_FPKM_Rank__ (integer) : megadja, hogy az adott izoforma expressziós értéke hanyadik helyre van sorolva az összes izoforma expressziós értéke közül az adott fejlődési fázisban/szövettípusban, megmutatja hogy az adott izoforma mennyire domináns az adott fejlődési fázisban/szövettípusban
- __Rank_Sum__ (integer) : az izoformák dominancia értékei/rank-jei összeadva az összes fejlődési fázisban/szövettípusban
- __Overall_FPKM__ (integer) : az izoformák összesített expressziós értéke az összes fejlődési fázisban/szövettípusban
- __GENE-FB-devreg__ (boolean)
  - a gén legalább 4 FPKM értékkel rendelkezik P1 - FB között valamelyik fejlődési fázisban/szövettípusban
  - a gén legalább 4x Fold Change-el rendelkezik valamelyik fejlődési útvonalon
    - DEVREG_S : P1 -> P2_S -> YFB_S -> FB_S
    - DEVREG_C : P1 -> P2_C -> YFB_C -> FB_C
    - DEVREG_L : P1 -> P2_C -> YFB_C -> FB_L
    - DEVREG_P2 : P2_S -> P2_C
    - DEVREG_YFB : YFB_S -> YFB_C
    - DEVREG_FB : FB_S -> FB_C -> FB_L
- __GENE-FB-init__ (boolean)
  - a gén legalább 4 FPKM értékkel rendelkezik VM - P1 között valamelyik fejlődési fázisban/szövettípusban
  - a gén legalább 4x pozitív Fold Change-el rendelkezik VM - P1 átmenetben
- __GENE-devreg-s-stricto__ (boolean) : vagy a GENE-FB-devreg, vagy a GENE-FB-init TRUE-t vesz fel
- __GENE-isoforms__ (integer) : a gén izoformáinak mennyisége

### Output: aostoyae_corrected_annotation.gtf
- újrarendezzük az annotációs fájl sorait az izoforma dominancia alapján, hogy a primary transcript legyen elől, az ASpli ezután a primary transcripthez viszonyítva adja meg az alternative splicing eseményeket

## Statisztikák kinyerése

### Script:
- ANNOTATIONS_aostoyae.R

### Input:
- __aostoyae_genes.tsv__
- __aostoyae_isoforms.tsv__
- __aostoyae_corrected_annotation.gtf__

### Output:
- __aostoyae_AS_genes.gtf__
  - alternative splicingban résztvevő géneket tartalmazó annotációs fájl az ASpli számára
- __aostoyae_FB_DEVREG_genes.gtf__
  - alternative splicingban résztvevő FB-devreg géneket tartalmazó annotációs fájl az ASpli számára
- __aostoyae_FB_INIT_genes.gtf__
  - alternative splicingban résztvevő FB-init géneket tartalmazó annotációs fájl az ASpli számára

## Alternative Splicing statisztikák kinyerése

### Script:
- STATISTICS_aostoyae.R

### Input:
- __aostoyae_AS_genes.gtf__
- __aostoyae_FB_DEVREG_genes.gtf__
- __aostoyae_FB_INIT_genes.gtf__

### Output:
- __aostoyae_stats.tsv__
  - __All__ __Genes__ : összes annotált gén
  - __FB-devreg__ __Genes__ : FB-devreg gének száma | FB-devreg gének százaléka az összes annotált génhez viszonyítva
  - __FB-init__ __Genes__ : FB-init gének száma | FB-init gének százaléka az összes annotált génhez viszonyítva
  - __AS__ __Genes__ : alternative splicingban résztvevő gének száma | alternative splicingban résztvevő gének százaléka az összes annotált génhez viszonyítva
  - __AS__ __FB-devreg__ __Genes__ : alternative splicingban résztvevő FB-devreg gének száma | alternative splicingban résztvevő FB-devreg gének százaléka az összes annotált génhez viszonyítva
  - __AS__ __FB-init__ __Genes__ : alternative splicingban résztvevő FB-init gének száma | alternative splicingban résztvevő FB-init gének százaléka az összes annotált génhez viszonyítva
  - __FB-devreg__ __and__ __FB-init__ __Genes__ : FB-devreg és FB-init gének száma | FB-devreg és FB-init gének százaléka az összes annotált génhez viszonyítva
  - __AS__ __FB-devreg__ __and__ __AS__ __FB-init__ __Genes__ : alternative splicingban résztvevő FB-devreg és FB-init gének száma | alternative splicingban résztvevő FB-devreg és FB-init gének százaléka az összes annotált génhez viszonyítva
  - __AS__ __Genes__ __--__ __ES__/__IR__/__Alt__ __5'__ __SS__/__Alt__ __3'__ __SS__ __bins__ : alternative splicingban résztvevő gének AS eseményeinek száma | alternative splicingban résztvevő gének AS eseményeinek százalékos aránya
  - __AS__ __FB-devreg__ __Genes__ __--__ __ES__/__IR__/__Alt__ __5'__ __SS__/__Alt__ __3'__ __SS__ __bins__ : alternative splicingban résztvevő FB-devreg gének AS eseményeinek száma | alternative splicingban résztvevő FB-devreg gének AS eseményeinek százalékos aránya
  - __AS__ __FB_init__ __Genes__ __--__ __ES__/__IR__/__Alt__ __5'__ __SS__/__Alt__ __3'__ __SS__ __bins__ : alternative splicingban résztvevő FB-init gének AS eseményeinek száma | alternative splicingban résztvevő FB-init gének AS eseményeinek százalékos aránya

# Coprinopsis cinerea

## Táblázat elkészítése

### Script:
- TABLE_ccinerea.R

### Input:
- __isoforms.fpkm_tracking__
  - CuffDiff output FPKM érték táblázat az izoformákra
- __genes.fpkm_tracking__
  - CuffDiff output FPKM érték táblázat az génekre
- __ccinerea_AmutBmut_AS_annotation.gtf__
  - a végső annotációs fájl amivel a CuffDiff analízis is futott

### Output:
- __ccinerea_isoforms.tsv__
  - __VM_FPKM_Rank__ -> __FB_T_FPKM_Rank__ __oszlopok__: megadja, hogy az adott izoforma expressziós értéke hanyadik helyre van sorolva az összes izoforma expressziós értéke közül az adott fejlődési fázisban/szövettípusban, megmutatja hogy az adott izoforma mennyire domináns az adott fejlődési fázisban/szövettípusban
  - __Rank_Sum__ __oszlop__: az izoformák dominancia értékei/rank-jei összeadva az összes fejlődési fázisban/szövettípusban
  - __Overall_FPKM__ __oszlok__: az izoformák teljes expressziós értéke az összes fejlődési fázisban/szövettípusban
  - __isoforms__ : a gén izoformáinak mennyisége
  - __isoform__ __switch__ : a rank-ok különböző fejlődési fázisok/szövettípusok közötti váltakozása miatt, megfigyelhető hogy egyes izoformák dominánsak-e egy fázisban majd a dominanciájukat elveszítik-e egy másikban
  - __primary__ __transcript__ : azt a transzkriptet nevezzük primary transzkriptnek amelyik összeadott dominancia értéke (Rank_Sum) a legkisebb, vagyis a legdominánsabb az izoformák közül, amennyiben ezt az értéket több izoforma is felveszi, az izoforma összeexpressziójának értéke (Overall_FPKM) alapján döntünk, vagyis az lesz a primary transcript amelyiknek a legnagyobb az expressziós értéke az azonos rank-ú izoformák közül
- __ccinerea_genes.tsv__
  - __FB-devreg__
    - legalább 4 FPKM értékkel rendelkezik H - FB között valamelyik fejlődési fázisban/szövettípusban
    - legalább 4x Fold Change-el rendelkezik valamelyik fejlődési útvonalon
      - DEVREG_T : H -> P1 -> P2 -> YFB_T -> FB_T
      - DEVREG_K : H -> P1 -> P2 -> YFB_K -> FB_KL
      - DEVREG_L : H -> P1 -> P2 -> YFB_L -> FB_KL
      - DEVREG_YFB : YFB_T -> YFB_K -> YFB_L
      - DEVREG_FB : FB_T -> FB_KL
  - __FB-init__
    - legalább 4 FPKM értékkel rendelkezik VM - P1 között valamelyik fejlődési fázisban/szövettípusban
    - legalább 4x pozitív Fold Change-el rendelkezik VM - H+P1 átmenetben
- __ccinerea_corrected_annotation.gtf__
  - újrarendezzük az annotációs fájl sorait az izoforma dominancia alapján, hogy a primary transcript legyen elől, az ASpli ezután a primary transcripthez viszonyítva adja meg az alternative splicing eseményeket

## Statisztikák kinyerése

### Script:
- ANNOTATIONS_ccinerea.R

### Input:
- __ccinerea_genes.tsv__
- __ccinerea_isoforms.tsv__
- __ccinerea_corrected_annotation.gtf__

### Output:
- __ccinerea_AS_genes.gtf__
  - alternative splicingban résztvevő géneket tartalmazó annotációs fájl az ASpli számára
- __ccinerea_FB_DEVREG_genes.gtf__
  - alternative splicingban résztvevő FB-devreg géneket tartalmazó annotációs fájl az ASpli számára
- __ccinerea_FB_INIT_genes.gtf__
  - alternative splicingban résztvevő FB-init géneket tartalmazó annotációs fájl az ASpli számára

## Alternative Splicing statisztikák kinyerése

### Script:
- STATISTICS_ccinerea.R

### Input:
- __ccinerea_AS_genes.gtf__
- __ccinerea_FB_DEVREG_genes.gtf__
- __ccinerea_FB_INIT_genes.gtf__

### Output:
- __ccinerea_stats.tsv__
  - __All__ __Genes__ : összes annotált gén
  - __FB-devreg__ __Genes__ : FB-devreg gének száma | FB-devreg gének százaléka az összes annotált génhez viszonyítva
  - __FB-init__ __Genes__ : FB-init gének száma | FB-init gének százaléka az összes annotált génhez viszonyítva
  - __AS__ __Genes__ : alternative splicingban résztvevő gének száma | alternative splicingban résztvevő gének százaléka az összes annotált génhez viszonyítva
  - __AS__ __FB-devreg__ __Genes__ : alternative splicingban résztvevő FB-devreg gének száma | alternative splicingban résztvevő FB-devreg gének százaléka az összes annotált génhez viszonyítva
  - __AS__ __FB-init__ __Genes__ : alternative splicingban résztvevő FB-init gének száma | alternative splicingban résztvevő FB-init gének százaléka az összes annotált génhez viszonyítva
  - __FB-devreg__ __and__ __FB-init__ __Genes__ : FB-devreg és FB-init gének száma | FB-devreg és FB-init gének százaléka az összes annotált génhez viszonyítva
  - __AS__ __FB-devreg__ __and__ __AS__ __FB-init__ __Genes__ : alternative splicingban résztvevő FB-devreg és FB-init gének száma | alternative splicingban résztvevő FB-devreg és FB-init gének százaléka az összes annotált génhez viszonyítva
  - __AS__ __Genes__ __--__ __ES__/__IR__/__Alt__ __5'__ __SS__/__Alt__ __3'__ __SS__ __bins__ : alternative splicingban résztvevő gének AS eseményeinek száma | alternative splicingban résztvevő gének AS eseményeinek százalékos aránya
  - __AS__ __FB-devreg__ __Genes__ __--__ __ES__/__IR__/__Alt__ __5'__ __SS__/__Alt__ __3'__ __SS__ __bins__ : alternative splicingban résztvevő FB-devreg gének AS eseményeinek száma | alternative splicingban résztvevő FB-devreg gének AS eseményeinek százalékos aránya
  - __AS__ __FB_init__ __Genes__ __--__ __ES__/__IR__/__Alt__ __5'__ __SS__/__Alt__ __3'__ __SS__ __bins__ : alternative splicingban résztvevő FB-init gének AS eseményeinek száma | alternative splicingban résztvevő FB-init gének AS eseményeinek százalékos aránya

# Lentinus tigrinus

## Táblázat elkészítése

### Script:
- TABLE_ltigrinus.R

### Input:
- __isoforms.fpkm_tracking__
  - CuffDiff output FPKM érték táblázat az izoformákra
- __genes.fpkm_tracking__
  - CuffDiff output FPKM érték táblázat az génekre
- __ltigrinus_AS_annotation.gtf__
  - a végső annotációs fájl amivel a CuffDiff analízis is futott

### Output:
- __ltigrinus_isoforms.tsv__
  - __VM_FPKM_Rank__ -> __FB_T_FPKM_Rank__ __oszlopok__: megadja, hogy az adott izoforma expressziós értéke hanyadik helyre van sorolva az összes izoforma expressziós értéke közül az adott fejlődési fázisban/szövettípusban, megmutatja hogy az adott izoforma mennyire domináns az adott fejlődési fázisban/szövettípusban
  - __Rank_Sum__ __oszlop__: az izoformák dominancia értékei/rank-jei összeadva az összes fejlődési fázisban/szövettípusban
  - __Overall_FPKM__ __oszlok__: az izoformák teljes expressziós értéke az összes fejlődési fázisban/szövettípusban
  - __isoforms__ : a gén izoformáinak mennyisége
  - __isoform__ __switch__ : a rank-ok különböző fejlődési fázisok/szövettípusok közötti váltakozása miatt, megfigyelhető hogy egyes izoformák dominánsak-e egy fázisban majd a dominanciájukat elveszítik-e egy másikban
  - __primary__ __transcript__ : azt a transzkriptet nevezzük primary transzkriptnek amelyik összeadott dominancia értéke (Rank_Sum) a legkisebb, vagyis a legdominánsabb az izoformák közül, amennyiben ezt az értéket több izoforma is felveszi, az izoforma összeexpressziójának értéke (Overall_FPKM) alapján döntünk, vagyis az lesz a primary transcript amelyiknek a legnagyobb az expressziós értéke az azonos rank-ú izoformák közül
- __ltigrinus_genes.tsv__
  - __FB-devreg__
    - legalább 4 FPKM értékkel rendelkezik P1 - FB között valamelyik fejlődési fázisban/szövettípusban
    - legalább 4x Fold Change-el rendelkezik valamelyik fejlődési útvonalon
      - DEVREG_T : VM -> P1 -> P2_T -> YFB_T -> FB_T
      - DEVREG_K : VM -> P1 -> P2_K -> YFB_K -> FB_K
      - DEVREG_L : VM -> P1 -> P2_K -> YFB_K -> FB_L
      - DEVREG_P2 : P2_T -> P2_K
      - DEVREG_YFB : YFB_T -> YFB_K
      - DEVREG_FB : FB_T -> FB_K -> FB_L
  - __FB-init__
    - legalább 4 FPKM értékkel rendelkezik VM - P1 között valamelyik fejlődési fázisban/szövettípusban
    - legalább 4x pozitív Fold Change-el rendelkezik VM - P1 átmenetben
- __ltigrinus_corrected_annotation.gtf__
  - újrarendezzük az annotációs fájl sorait az izoforma dominancia alapján, hogy a primary transcript legyen elől, az ASpli ezután a primary transcripthez viszonyítva adja meg az alternative splicing eseményeket

## Statisztikák kinyerése

### Script:
- ANNOTATIONS_ltigrinus.R

### Input:
- __ltigrinus_genes.tsv__
- __ltigrinus_isoforms.tsv__
- __ltigrinus_corrected_annotation.gtf__

### Output:
- __ltigrinus_AS_genes.gtf__
  - alternative splicingban résztvevő géneket tartalmazó annotációs fájl az ASpli számára
- __ltigrinus_FB_DEVREG_genes.gtf__
  - alternative splicingban résztvevő FB-devreg géneket tartalmazó annotációs fájl az ASpli számára
- __ltigrinus_FB_INIT_genes.gtf__
  - alternative splicingban résztvevő FB-init géneket tartalmazó annotációs fájl az ASpli számára

## Alternative Splicing statisztikák kinyerése

### Script:
- STATISTICS_ltigrinus.R

### Input:
- __ltigrinus_AS_genes.gtf__
- __ltigrinus_FB_DEVREG_genes.gtf__
- __ltigrinus_FB_INIT_genes.gtf__

### Output:
- __ltigrinus_stats.tsv__
  - __All__ __Genes__ : összes annotált gén
  - __FB-devreg__ __Genes__ : FB-devreg gének száma | FB-devreg gének százaléka az összes annotált génhez viszonyítva
  - __FB-init__ __Genes__ : FB-init gének száma | FB-init gének százaléka az összes annotált génhez viszonyítva
  - __AS__ __Genes__ : alternative splicingban résztvevő gének száma | alternative splicingban résztvevő gének százaléka az összes annotált génhez viszonyítva
  - __AS__ __FB-devreg__ __Genes__ : alternative splicingban résztvevő FB-devreg gének száma | alternative splicingban résztvevő FB-devreg gének százaléka az összes annotált génhez viszonyítva
  - __AS__ __FB-init__ __Genes__ : alternative splicingban résztvevő FB-init gének száma | alternative splicingban résztvevő FB-init gének százaléka az összes annotált génhez viszonyítva
  - __FB-devreg__ __and__ __FB-init__ __Genes__ : FB-devreg és FB-init gének száma | FB-devreg és FB-init gének százaléka az összes annotált génhez viszonyítva
  - __AS__ __FB-devreg__ __and__ __AS__ __FB-init__ __Genes__ : alternative splicingban résztvevő FB-devreg és FB-init gének száma | alternative splicingban résztvevő FB-devreg és FB-init gének százaléka az összes annotált génhez viszonyítva
  - __AS__ __Genes__ __--__ __ES__/__IR__/__Alt__ __5'__ __SS__/__Alt__ __3'__ __SS__ __bins__ : alternative splicingban résztvevő gének AS eseményeinek száma | alternative splicingban résztvevő gének AS eseményeinek százalékos aránya
  - __AS__ __FB-devreg__ __Genes__ __--__ __ES__/__IR__/__Alt__ __5'__ __SS__/__Alt__ __3'__ __SS__ __bins__ : alternative splicingban résztvevő FB-devreg gének AS eseményeinek száma | alternative splicingban résztvevő FB-devreg gének AS eseményeinek százalékos aránya
  - __AS__ __FB_init__ __Genes__ __--__ __ES__/__IR__/__Alt__ __5'__ __SS__/__Alt__ __3'__ __SS__ __bins__ : alternative splicingban résztvevő FB-init gének AS eseményeinek száma | alternative splicingban résztvevő FB-init gének AS eseményeinek százalékos aránya

# Rickenella mellea

## Táblázat elkészítése

### Script:
- TABLE_rmellea.R

### Input:
- __isoforms.fpkm_tracking__
  - CuffDiff output FPKM érték táblázat az izoformákra
- __genes.fpkm_tracking__
  - CuffDiff output FPKM érték táblázat az génekre
- __rmellea_AS_annotation.gtf__
  - a végső annotációs fájl amivel a CuffDiff analízis is futott

### Output:
- __rmellea_isoforms.tsv__
  - __VM_FPKM_Rank__ -> __FB_T_FPKM_Rank__ __oszlopok__: megadja, hogy az adott izoforma expressziós értéke hanyadik helyre van sorolva az összes izoforma expressziós értéke közül az adott fejlődési fázisban/szövettípusban, megmutatja hogy az adott izoforma mennyire domináns az adott fejlődési fázisban/szövettípusban
  - __Rank_Sum__ __oszlop__: az izoformák dominancia értékei/rank-jei összeadva az összes fejlődési fázisban/szövettípusban
  - __Overall_FPKM__ __oszlok__: az izoformák teljes expressziós értéke az összes fejlődési fázisban/szövettípusban
  - __isoforms__ : a gén izoformáinak mennyisége
  - __isoform__ __switch__ : a rank-ok különböző fejlődési fázisok/szövettípusok közötti váltakozása miatt, megfigyelhető hogy egyes izoformák dominánsak-e egy fázisban majd a dominanciájukat elveszítik-e egy másikban
  - __primary__ __transcript__ : azt a transzkriptet nevezzük primary transzkriptnek amelyik összeadott dominancia értéke (Rank_Sum) a legkisebb, vagyis a legdominánsabb az izoformák közül, amennyiben ezt az értéket több izoforma is felveszi, az izoforma összeexpressziójának értéke (Overall_FPKM) alapján döntünk, vagyis az lesz a primary transcript amelyiknek a legnagyobb az expressziós értéke az azonos rank-ú izoformák közül
- __rmellea_genes.tsv__
  - __FB-devreg__
    - legalább 4 FPKM értékkel rendelkezik P - FB között valamelyik fejlődési fázisban/szövettípusban
    - legalább 4x Fold Change-el rendelkezik valamelyik fejlődési útvonalon
      - DEVREG_T : VM -> P -> YFB -> FB_T
      - DEVREG_K : VM -> P -> YFB -> FB_K
      - DEVREG_FB : FB_T -> FB_K
  - __FB-init__
    - legalább 4 FPKM értékkel rendelkezik VM - P között valamelyik fejlődési fázisban/szövettípusban
    - legalább 4x pozitív Fold Change-el rendelkezik VM - P átmenetben
- __rmellea_corrected_annotation.gtf__
  - újrarendezzük az annotációs fájl sorait az izoforma dominancia alapján, hogy a primary transcript legyen elől, az ASpli ezután a primary transcripthez viszonyítva adja meg az alternative splicing eseményeket

## Statisztikák kinyerése

### Script:
- ANNOTATIONS_rmellea.R

### Input:
- __rmellea_genes.tsv__
- __rmellea_isoforms.tsv__
- __rmellea_corrected_annotation.gtf__

### Output:
- __rmellea_AS_genes.gtf__
  - alternative splicingban résztvevő géneket tartalmazó annotációs fájl az ASpli számára
- __rmellea_FB_DEVREG_genes.gtf__
  - alternative splicingban résztvevő FB-devreg géneket tartalmazó annotációs fájl az ASpli számára
- __rmellea_FB_INIT_genes.gtf__
  - alternative splicingban résztvevő FB-init géneket tartalmazó annotációs fájl az ASpli számára

## Alternative Splicing statisztikák kinyerése

### Script:
- STATISTICS_rmellea.R

### Input:
- __rmellea_AS_genes.gtf__
- __rmellea_FB_DEVREG_genes.gtf__
- __rmellea_FB_INIT_genes.gtf__

### Output:
- __rmellea_stats.tsv__
  - __All__ __Genes__ : összes annotált gén
  - __FB-devreg__ __Genes__ : FB-devreg gének száma | FB-devreg gének százaléka az összes annotált génhez viszonyítva
  - __FB-init__ __Genes__ : FB-init gének száma | FB-init gének százaléka az összes annotált génhez viszonyítva
  - __AS__ __Genes__ : alternative splicingban résztvevő gének száma | alternative splicingban résztvevő gének százaléka az összes annotált génhez viszonyítva
  - __AS__ __FB-devreg__ __Genes__ : alternative splicingban résztvevő FB-devreg gének száma | alternative splicingban résztvevő FB-devreg gének százaléka az összes annotált génhez viszonyítva
  - __AS__ __FB-init__ __Genes__ : alternative splicingban résztvevő FB-init gének száma | alternative splicingban résztvevő FB-init gének százaléka az összes annotált génhez viszonyítva
  - __FB-devreg__ __and__ __FB-init__ __Genes__ : FB-devreg és FB-init gének száma | FB-devreg és FB-init gének százaléka az összes annotált génhez viszonyítva
  - __AS__ __FB-devreg__ __and__ __AS__ __FB-init__ __Genes__ : alternative splicingban résztvevő FB-devreg és FB-init gének száma | alternative splicingban résztvevő FB-devreg és FB-init gének százaléka az összes annotált génhez viszonyítva
  - __AS__ __Genes__ __--__ __ES__/__IR__/__Alt__ __5'__ __SS__/__Alt__ __3'__ __SS__ __bins__ : alternative splicingban résztvevő gének AS eseményeinek száma | alternative splicingban résztvevő gének AS eseményeinek százalékos aránya
  - __AS__ __FB-devreg__ __Genes__ __--__ __ES__/__IR__/__Alt__ __5'__ __SS__/__Alt__ __3'__ __SS__ __bins__ : alternative splicingban résztvevő FB-devreg gének AS eseményeinek száma | alternative splicingban résztvevő FB-devreg gének AS eseményeinek százalékos aránya
  - __AS__ __FB_init__ __Genes__ __--__ __ES__/__IR__/__Alt__ __5'__ __SS__/__Alt__ __3'__ __SS__ __bins__ : alternative splicingban résztvevő FB-init gének AS eseményeinek száma | alternative splicingban résztvevő FB-init gének AS eseményeinek százalékos aránya

# Schizophyllum commune

## Táblázat elkészítése

### Script:
- TABLE_scommune.R

### Input:
- __isoforms.fpkm_tracking__
  - CuffDiff output FPKM érték táblázat az izoformákra
- __genes.fpkm_tracking__
  - CuffDiff output FPKM érték táblázat az génekre
- __scommune_AS_annotation.gtf__
  - a végső annotációs fájl amivel a CuffDiff analízis is futott

### Output:
- __scommune_isoforms.tsv__
  - __VM_FPKM_Rank__ -> __FB_FPKM_Rank__ __oszlopok__: megadja, hogy az adott izoforma expressziós értéke hanyadik helyre van sorolva az összes izoforma expressziós értéke közül az adott fejlődési fázisban/szövettípusban, megmutatja hogy az adott izoforma mennyire domináns az adott fejlődési fázisban/szövettípusban
  - __Rank_Sum__ __oszlop__: az izoformák dominancia értékei/rank-jei összeadva az összes fejlődési fázisban/szövettípusban
  - __Overall_FPKM__ __oszlok__: az izoformák teljes expressziós értéke az összes fejlődési fázisban/szövettípusban
  - __isoforms__ : a gén izoformáinak mennyisége
  - __isoform__ __switch__ : a rank-ok különböző fejlődési fázisok/szövettípusok közötti váltakozása miatt, megfigyelhető hogy egyes izoformák dominánsak-e egy fázisban majd a dominanciájukat elveszítik-e egy másikban
  - __primary__ __transcript__ : azt a transzkriptet nevezzük primary transzkriptnek amelyik összeadott dominancia értéke (Rank_Sum) a legkisebb, vagyis a legdominánsabb az izoformák közül, amennyiben ezt az értéket több izoforma is felveszi, az izoforma összeexpressziójának értéke (Overall_FPKM) alapján döntünk, vagyis az lesz a primary transcript amelyiknek a legnagyobb az expressziós értéke az azonos rank-ú izoformák közül
- __scommune_genes.tsv__
  - __FB-devreg__
    - legalább 4 FPKM értékkel rendelkezik P1 - FB között valamelyik fejlődési fázisban/szövettípusban
    - legalább 4x Fold Change-el rendelkezik valamelyik fejlődési útvonalon
      - DEVREG : VM -> P1 -> YFB -> FB
  - __FB-init__
    - legalább 4 FPKM értékkel rendelkezik VM - P1 között valamelyik fejlődési fázisban/szövettípusban
    - legalább 4x pozitív Fold Change-el rendelkezik VM - P1 átmenetben
- __scommune_corrected_annotation.gtf__
  - újrarendezzük az annotációs fájl sorait az izoforma dominancia alapján, hogy a primary transcript legyen elől, az ASpli ezután a primary transcripthez viszonyítva adja meg az alternative splicing eseményeket

## Statisztikák kinyerése

### Script:
- ANNOTATIONS_scommune.R

### Input:
- __scommune_genes.tsv__
- __scommune_isoforms.tsv__
- __scommune_corrected_annotation.gtf__

### Output:
- __scommune_AS_genes.gtf__
  - alternative splicingban résztvevő géneket tartalmazó annotációs fájl az ASpli számára
- __scommune_FB_DEVREG_genes.gtf__
  - alternative splicingban résztvevő FB-devreg géneket tartalmazó annotációs fájl az ASpli számára
- __scommune_FB_INIT_genes.gtf__
  - alternative splicingban résztvevő FB-init géneket tartalmazó annotációs fájl az ASpli számára

## Alternative Splicing statisztikák kinyerése

### Script:
- STATISTICS_scommune.R

### Input:
- __scommune_AS_genes.gtf__
- __scommune_FB_DEVREG_genes.gtf__
- __scommune_FB_INIT_genes.gtf__

### Output:
- __scommune_stats.tsv__
  - __All__ __Genes__ : összes annotált gén
  - __FB-devreg__ __Genes__ : FB-devreg gének száma | FB-devreg gének százaléka az összes annotált génhez viszonyítva
  - __FB-init__ __Genes__ : FB-init gének száma | FB-init gének százaléka az összes annotált génhez viszonyítva
  - __AS__ __Genes__ : alternative splicingban résztvevő gének száma | alternative splicingban résztvevő gének százaléka az összes annotált génhez viszonyítva
  - __AS__ __FB-devreg__ __Genes__ : alternative splicingban résztvevő FB-devreg gének száma | alternative splicingban résztvevő FB-devreg gének százaléka az összes annotált génhez viszonyítva
  - __AS__ __FB-init__ __Genes__ : alternative splicingban résztvevő FB-init gének száma | alternative splicingban résztvevő FB-init gének százaléka az összes annotált génhez viszonyítva
  - __FB-devreg__ __and__ __FB-init__ __Genes__ : FB-devreg és FB-init gének száma | FB-devreg és FB-init gének százaléka az összes annotált génhez viszonyítva
  - __AS__ __FB-devreg__ __and__ __AS__ __FB-init__ __Genes__ : alternative splicingban résztvevő FB-devreg és FB-init gének száma | alternative splicingban résztvevő FB-devreg és FB-init gének százaléka az összes annotált génhez viszonyítva
  - __AS__ __Genes__ __--__ __ES__/__IR__/__Alt__ __5'__ __SS__/__Alt__ __3'__ __SS__ __bins__ : alternative splicingban résztvevő gének AS eseményeinek száma | alternative splicingban résztvevő gének AS eseményeinek százalékos aránya
  - __AS__ __FB-devreg__ __Genes__ __--__ __ES__/__IR__/__Alt__ __5'__ __SS__/__Alt__ __3'__ __SS__ __bins__ : alternative splicingban résztvevő FB-devreg gének AS eseményeinek száma | alternative splicingban résztvevő FB-devreg gének AS eseményeinek százalékos aránya
  - __AS__ __FB_init__ __Genes__ __--__ __ES__/__IR__/__Alt__ __5'__ __SS__/__Alt__ __3'__ __SS__ __bins__ : alternative splicingban résztvevő FB-init gének AS eseményeinek száma | alternative splicingban résztvevő FB-init gének AS eseményeinek százalékos aránya
  
# Phanerochaete chrysosporium

## Táblázat elkészítése

### Script:
- TABLE_pchrysosporium.R

### Input:
- __isoforms.fpkm_tracking__
  - CuffDiff output FPKM érték táblázat az izoformákra
- __genes.fpkm_tracking__
  - CuffDiff output FPKM érték táblázat az génekre
- __pchrysosporium_AS_annotation.gtf__
  - a végső annotációs fájl amivel a CuffDiff analízis is futott

### Output:
- __pchrysosporium_isoforms.tsv__
  - __VM_FPKM_Rank__ -> __FB_FPKM_Rank__ __oszlopok__: megadja, hogy az adott izoforma expressziós értéke hanyadik helyre van sorolva az összes izoforma expressziós értéke közül az adott fejlődési fázisban/szövettípusban, megmutatja hogy az adott izoforma mennyire domináns az adott fejlődési fázisban/szövettípusban
  - __Rank_Sum__ __oszlop__: az izoformák dominancia értékei/rank-jei összeadva az összes fejlődési fázisban/szövettípusban
  - __Overall_FPKM__ __oszlok__: az izoformák teljes expressziós értéke az összes fejlődési fázisban/szövettípusban
  - __isoforms__ : a gén izoformáinak mennyisége
  - __isoform__ __switch__ : a rank-ok különböző fejlődési fázisok/szövettípusok közötti váltakozása miatt, megfigyelhető hogy egyes izoformák dominánsak-e egy fázisban majd a dominanciájukat elveszítik-e egy másikban
  - __primary__ __transcript__ : azt a transzkriptet nevezzük primary transzkriptnek amelyik összeadott dominancia értéke (Rank_Sum) a legkisebb, vagyis a legdominánsabb az izoformák közül, amennyiben ezt az értéket több izoforma is felveszi, az izoforma összeexpressziójának értéke (Overall_FPKM) alapján döntünk, vagyis az lesz a primary transcript amelyiknek a legnagyobb az expressziós értéke az azonos rank-ú izoformák közül
- __pchrysosporium_genes.tsv__
  - __FB-devreg__
    - legalább 4 FPKM értékkel rendelkezik YFB - FB között valamelyik fejlődési fázisban/szövettípusban
    - legalább 4x Fold Change-el rendelkezik valamelyik fejlődési útvonalon
      - DEVREG : VM -> YFB -> FB
  - __FB-init__
    - legalább 4 FPKM értékkel rendelkezik VM - YFB között valamelyik fejlődési fázisban/szövettípusban
    - legalább 4x pozitív Fold Change-el rendelkezik VM - YFB átmenetben
- __pchrysosporium_corrected_annotation.gtf__
  - újrarendezzük az annotációs fájl sorait az izoforma dominancia alapján, hogy a primary transcript legyen elől, az ASpli ezután a primary transcripthez viszonyítva adja meg az alternative splicing eseményeket

## Statisztikák kinyerése

### Script:
- ANNOTATIONS_pchrysosporium.R

### Input:
- __pchrysosporium_genes.tsv__
- __pchrysosporium_isoforms.tsv__
- __pchrysosporium_corrected_annotation.gtf__

### Output:
- __pchrysosporium_AS_genes.gtf__
  - alternative splicingban résztvevő géneket tartalmazó annotációs fájl az ASpli számára
- __pchrysosporium_FB_DEVREG_genes.gtf__
  - alternative splicingban résztvevő FB-devreg géneket tartalmazó annotációs fájl az ASpli számára
- __pchrysosporium_FB_INIT_genes.gtf__
  - alternative splicingban résztvevő FB-init géneket tartalmazó annotációs fájl az ASpli számára

## Alternative Splicing statisztikák kinyerése

### Script:
- STATISTICS_pchrysosporium.R

### Input:
- __pchrysosporium_AS_genes.gtf__
- __pchrysosporium_FB_DEVREG_genes.gtf__
- __pchrysosporium_FB_INIT_genes.gtf__

### Output:
- __pchrysosporium_stats.tsv__
  - __All__ __Genes__ : összes annotált gén
  - __FB-devreg__ __Genes__ : FB-devreg gének száma | FB-devreg gének százaléka az összes annotált génhez viszonyítva
  - __FB-init__ __Genes__ : FB-init gének száma | FB-init gének százaléka az összes annotált génhez viszonyítva
  - __AS__ __Genes__ : alternative splicingban résztvevő gének száma | alternative splicingban résztvevő gének százaléka az összes annotált génhez viszonyítva
  - __AS__ __FB-devreg__ __Genes__ : alternative splicingban résztvevő FB-devreg gének száma | alternative splicingban résztvevő FB-devreg gének százaléka az összes annotált génhez viszonyítva
  - __AS__ __FB-init__ __Genes__ : alternative splicingban résztvevő FB-init gének száma | alternative splicingban résztvevő FB-init gének százaléka az összes annotált génhez viszonyítva
  - __FB-devreg__ __and__ __FB-init__ __Genes__ : FB-devreg és FB-init gének száma | FB-devreg és FB-init gének százaléka az összes annotált génhez viszonyítva
  - __AS__ __FB-devreg__ __and__ __AS__ __FB-init__ __Genes__ : alternative splicingban résztvevő FB-devreg és FB-init gének száma | alternative splicingban résztvevő FB-devreg és FB-init gének százaléka az összes annotált génhez viszonyítva
  - __AS__ __Genes__ __--__ __ES__/__IR__/__Alt__ __5'__ __SS__/__Alt__ __3'__ __SS__ __bins__ : alternative splicingban résztvevő gének AS eseményeinek száma | alternative splicingban résztvevő gének AS eseményeinek százalékos aránya
  - __AS__ __FB-devreg__ __Genes__ __--__ __ES__/__IR__/__Alt__ __5'__ __SS__/__Alt__ __3'__ __SS__ __bins__ : alternative splicingban résztvevő FB-devreg gének AS eseményeinek száma | alternative splicingban résztvevő FB-devreg gének AS eseményeinek százalékos aránya
  - __AS__ __FB_init__ __Genes__ __--__ __ES__/__IR__/__Alt__ __5'__ __SS__/__Alt__ __3'__ __SS__ __bins__ : alternative splicingban résztvevő FB-init gének AS eseményeinek száma | alternative splicingban résztvevő FB-init gének AS eseményeinek százalékos aránya
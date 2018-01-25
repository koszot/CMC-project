# Complex Multicellularity Project

# Feladatok:
- NEM KELL AZ AMPLA, SE A NEM SAJÁT ANYAGOK: cneoformans, umaydis
- Splice landscape (hány gén, hány esemény, milyen megoszlás) - alapadatokból kinyerhető, de finomítani kell
- 1. DEVREG - 4 FPKM VM - P 4x fold change
- 2. DEVREG - 4 FPKM P - FB 4x fold change 
- szignifikancia threshold felállítása izoformáknál: expresszió/izoformák száma -> minimális expressziós szint amire azt mondhatjuk hogy akkor most ez az izoforma szignifikáns
- splice események összehasonlítása a devreg - nem devreg között -- ehhez egy script ami kinyeri minden izoforma alternative splicing eseményeit
- megnézni hogy milyen változások vannak a DEVREG géneknél az izoformák megoszlásába, isoform-switch -> R érték számítása minden izoformára és erre egy R érték fold change kiszámolása
- UTR régiók összehasonlítása, 1-2 devreg, nem devreg

# Táblázat elkészítése

Az adatokat amik szükségesek egy táblázatba gyűjtjük bele, a kiindulási táblázat a CUFFDIFF outputjaiban lévő isoforms.fpkm táblázat amit egy R script segítségével számolunk tovább.
Táblázat leírása:
- __gene_id__ : mivel a táblázata alapvetően izoformák adatait tartalmazza ezért az első oszlop annak a génnek az ID-ja amiből leszármazik az aktuális izoforma
- __transcript_id__ : az aktuális izoforma ID-ja amihez az adatok tartoznak
- __VM_FPKM__ __(RMA_FPKM)__ __->__ __FB_FPKM__ : az izoforma FPKM értékei a különböző fejlődési fázisokban
- __N_isoforms__ : megadja, hogy az anyagénnek összesen hány izoformája van az aktuális izoformával együtt
- __gene_4FPKM__ : megadja, hogy az anyagén valamelyik fejlődési fázisban eléri-e a legalább 4FPKM értéket
- __gene_VMtoP_FC__ : megadja, hogy az anyagén mekkora fold change-el rendelkezik a VM-P átmentetben (a 4x fold changehez itt az kell hogy az érték 4 felett legyen)
- __gene_PtoFB_FC__ : megadja, hogy az anyagén mekkora fold change-el rendelkezik P-től FB-ig (ehhez vesszük a P-FB közötti fejlődési fázisok maximumát és minimumát és azok alapján határozzuk meg hogy mekkora a fold change a két legtávolabbi érték között)




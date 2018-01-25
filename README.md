# Complex Multicellularity Project

# Feladatok:
- NEM KELL AZ AMPLA, SE A NEM SAJÁT ANYAGOK: cneoformans, umaydis
- __DONE__ 1. DEVREG - 4 FPKM VM - P 4x fold change
- __DONE__ 2. DEVREG - 4 FPKM P - FB 4x fold change 
- __DONE__ szignifikancia threshold felállítása izoformáknál: expresszió/izoformák száma -> minimális expressziós szint amire azt mondhatjuk hogy akkor most ez az izoforma szignifikáns
- megnézni hogy milyen változások vannak a DEVREG géneknél az izoformák megoszlásába, isoform-switch -> R érték számítása minden izoformára és erre egy R érték fold change kiszámolása
- splice események összehasonlítása a devreg - nem devreg között -- ehhez egy script ami kinyeri minden izoforma alternative splicing eseményeit
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
- __gene_DEVREG_type1__ : a gén FPKM értéke valamelyik fázisban meghaladja a 4 FPKM-et és a VM - P átmenetben legalább 4x-es fold change-el rendelkezik
- __gene_DEVREG_type2__ : a gén FPKM értéke valamelyik fázisban meghaladja a 4 FPKM-et és a P - FB átmenetben legalább 4x-es fold change-el rendelkezik
- __isoform_Sign_VM__ -> __isoform_Sign_FB__ : megadja, hogy az adott isoforma szignifikáns-e az adott fázisban, ezt úgy számoljuk ki, hogy felállítunk egy thresholdot amit úgy kapunk meg hogy elosztjuk a gén teljes epresszióját, az izoformák számával így megkapjuk, hogyha egyenlő FPKM megoszlása lenne minden izoformának akkor mekkora FPKM értékkel kellene rendelkeznie, majd megnézzük hogy az egyes izoformák expressziós értéke meghaladja-e ezt a thresholdot, amennyiben igen TRUE-t kap, amennyiben nem FALSE-ot
- __isoform_maxSign__ : összeadja hogy az egyes izoformák összesen hány fejlődési állapotban tekinthetőek szignifikánsnak a threshold alapján, ezután sorbaállítjuk ez alapján az adott gén izoformáit és az első lesz így a primary isoform





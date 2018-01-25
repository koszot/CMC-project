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
- gene_id : mivel a táblázata alapvetően izoformák adatait tartalmazza ezért az első oszlop annak a génnek az ID-ja amiből leszármazik az aktuális izoforma
- transcript_id : az aktuális izoforma ID-ja amihez az adatok tartoznak
- VM_FPKM (RMA_FPKM) -> FB_FPKM : az izoforma FPKM értékei a különböző fejlődési fázisokban
 



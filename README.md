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
  - __VM_FPKM_Rank__ -> __FB_S_FPKM_Rank__
- aostoyae_genes.tsv
- aostoyae_corrected_annotation.gtf

### Script:
- STATS_TABLE_aostoyae.R


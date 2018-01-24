#######################
#    SILIX, HIFIX     #
#######################

# elvégezzük az all vs all blastot

nohup mpiformatdb -N 64 -i SMC_project_allgenomes.fasta -o T &
mv SMC_project_allgenomes.fasta.* ../share/
nohup mpirun -n 32 mpiblast -d SMC_project_allgenomes.fasta -i SMC_project_allgenomes.fasta -p blastp -m 8 -o SMC_project_allgenomes_blasted --copy-via=none &

# selenocysteine-es baszakodás

U_to_X.R

# Silix a Hifix-hez

silix SMC_project_allgenomes.fasta SMC_project_allgenomes_blasted --net > SMC_project_SLX.fnodes

# Hifix

hifix -t 64 SMC_project_allgenomes.fasta SMC_project_allgenomes_blasted.net SMC_project_SLX.fnodes> SMC_project_HFX.fnodes

# nem egyezik a szekvenciaszám a FASTA fileba és a HFX-be ---> több fehérje többször is előfordul

#################
#    NOTUNG     #
#################

# először megnyitjuk Notungba a fajfát és  Display Options -> Display Internal Node Names és elmentjük NEWICK formátumba (species.tree.nwk)
# A compare outputból kimásoljuk a reference részt ami a node neveket és az alatta található fajokat tartalmazza, majd egyeztetjük a NOTUNG-os fajfával a MRCA keresés alapján, ezzel megkapjuk, hogy melyik COMPARE node névnek mi a NOTUNGOS megfelelője

Dict.R 
!!!!!! header <- notung, compare, species

# bemásoljuk a species.tree.nwk-t a fákat tartalmazó mappába
# mivel a Cryne_H99_1 bezavar ezért ezt át kell nevezni az összes fába CryneH991-re

perl -pi -e 's/Cryne_H99_1/CryneH991/g' ./*

# elkészítjük a batch fájlt és behelyezzük a fákat tartalmazó mappába

java -jar /Users/genomelab3/programs/Notung-2.9/Notung-2.9.jar -b /Users/genomelab3/Desktop/SMC_project/Notung/trees_for_Notung/batch.run --reconcile --speciestag postfix --parsable --progressbar --rearrange --edgeweights length --threshold 70% --log --treeoutput newick --nolosses

# --maxtrees --> ha az összes fa kell

# parsable statisztikák átalakítása és egyesítése COMPARE formátumra

Notung_stats.R

Notung_stats_functional_groups.R

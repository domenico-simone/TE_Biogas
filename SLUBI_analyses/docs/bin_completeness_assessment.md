```bash
cd /home/domeni/theCellarDoor/SLUBI_Project_TE_Biogas_snic2019-30-23/SLUBI_analyses

mkdir bins_fasta
cd bins_fasta/
ln -sf ../../Binning/Bins/*.fa .
cd ..

# checkm taxonomy_wf domain Bacteria ./bins ./output
checkm taxonomy_wf -t 20 \
-x fa \
--force_overwrite \
-f bins_taxonomy_wf/bacteria/summary.tab \
domain Bacteria ./bins_fasta bins_taxonomy_wf/bacteria

checkm taxonomy_wf -t 20 \
-x fa \
--force_overwrite \
-f bins_taxonomy_wf/archaea/summary.tab \
domain Archaea ./bins_fasta bins_taxonomy_wf/archaea
```

# Gene enrichment analysis: EGGNOG mapping

<!-- TOC START min:2 max:5 link:true asterisk:false update:true -->
- [Collect all protein sequences](#collect-all-protein-sequences)
- [Functional annotation](#functional-annotation)
    - [Some stats](#some-stats)
- [Copy outputs locally...](#copy-outputs-locally)
- [Parse annotations](#parse-annotations)
    - [Get annotation file for GO](#get-annotation-file-for-go)
<!-- TOC END -->

Down to the `emapper.py` run, the following analyses are run on UPPMAX.

## Collect all protein sequences

```bash
export PATH=/home/domeni/theCellarDoor/tools/eggnog-mapper:$PATH

cd /home/domeni/theCellarDoor/SLUBI_Project_TE_Biogas_snic2019-30-23/SLUBI_analyses
mkdir -p gene_enrichment
cd gene_enrichment

cat /home/domeni/theCellarDoor/SLUBI_Project_TE_Biogas_snic2019-30-23/Binning/Bins/*/*.faa > all.proteins.fa
```

Get their lengths. We'll need them for enrichment analysis.

```python
from Bio import SeqIO
proteins = SeqIO.parse("all.proteins.fa", 'fasta')
out = open("all.proteins.lengths", 'w')
_ = out.write("gene\tlength\n")
for i in proteins:
    _ = out.write("{gene}\t{length}\n".format(gene=i.id, length=len(i)*3))

out.close()
```

## Functional annotation

Fetch eggnog mapper

```bash
git clone \
https://github.com/eggnogdb/eggnog-mapper.git \
~/home/domeni/theCellarDoor/tools
```

Run it

```bash
proj=snic2019-3-557

sbatch -A $proj -M snowy -t 30:00:00 -p node -J emapper -o /home/domeni/theCellarDoor/SLUBI_Project_TE_Biogas_snic2019-30-23/SLUBI_analyses/logs/emapper.log<<'EOF'
#!/bin/bash

cp -R ~/theCellarDoor/tools/eggnog-mapper ${SNIC_TMP}
cp all.proteins.fa ${SNIC_TMP}
mkdir -p ${SNIC_TMP}/tmp

${SNIC_TMP}/eggnog-mapper/emapper.py \
-i ${SNIC_TMP}/all.proteins.fa \
--output ${SNIC_TMP}/all.proteins.emapper.out \
--temp_dir ${SNIC_TMP}/tmp \
-m diamond --cpu 16

cp ${SNIC_TMP}/all.proteins.emapper.out* .

EOF
```

### Some stats

Total proteins:

```bash
grep -c ">" all.proteins.fa
```

```
114574
```

How many proteins with KO annotations?

```bash
grep -v "^#" all.proteins.emapper.out.emapper.annotations | awk 'BEGIN{FS="\t"}{print $1, $7, $9}' | grep -c ko
```

```
72866
```

How many proteins with GO annotations?

```bash
grep -v "^#" all.proteins.emapper.out.emapper.annotations | awk 'BEGIN{FS="\t"}{print $1, $7, $9}' | grep -c GO
```

```
25411
```

Then we should probably run the GSEA with both annotation sets.

## Copy outputs locally...

## Parse annotations

What we have from eggnog mapper:

```
AS_TE.10_00002   GO:0000271,GO:0003674,GO:0003824,GO:0005975
```

What we need:

```
AS_TE.10_00002   GO:0000271
AS_TE.10_00002   GO:0003674
AS_TE.10_00002   GO:0003824
AS_TE.10_00002   GO:0005975
```

Same for KO ids.

### Get annotation file for GO

```python
eggnog_out = open('gene_enrichment/all.proteins.emapper.out.emapper.annotations', 'r')

ko_annotations = open('gene_enrichment/all.proteins.emapper.ko', 'w')
go_annotations = open('gene_enrichment/all.proteins.emapper.go', 'w')

for l in eggnog_out:
    if l.startswith("#"):
        continue
    l = l.split("\t")
    gene_id = l[0]
    go_list = l[6].split(",")
    ko_list = l[8].split(",")
    if go_list[0] != "":
        for go in go_list:
            _ = go_annotations.write("{gene_id}\t{go}\n".format(gene_id=gene_id,
                                                            go=go))
    if ko_list[0] != "":
        for ko in ko_list:
            _ = ko_annotations.write("{gene_id}\t{ko}\n".format(gene_id=gene_id,
                                                            ko=ko))

ko_annotations.close()
go_annotations.close()
```

# Bin quantification in metagenomic samples

## Salmon

Installed v1.3 as conda env.

## Use only binned contigs

### Map reads with bowtie2

```bash
cd /home/domeni/theCellarDoor/SLUBI_Project_TE_Biogas_snic2019-30-23/SLUBI_analyses
mkdir -p bin_quantification
```

Collect all bin sequences, with contig names as _eg_ `AS_TE.75__k141_123`

```python
import glob, os
from Bio import SeqIO

# cat all bin sequences, but prepend bin name to contig names
bin_files = glob.glob("../Binning/Bins/AS*.fa")
outhandle = open("bin_quantification/all_bins.fa", 'w')
for bin in bin_files:
    print(bin)
    bin_name = os.path.split(bin)[1].replace(".fa", "")
    bin_handle = SeqIO.parse(bin, 'fasta')
    for contig in bin_handle:
        id_out = "{}__{}".format(bin_name, contig.id)
        _ = outhandle.write(">{}\n{}\n".format(id_out, str(contig.seq)))

outhandle.close()
```

Create bowtie2 index

```bash
bowtie2-build \
    --threads 20 \
    bin_quantification/all_bins.fa \
    bin_quantification/all_bins
```

Map reads

```bash
for i in GQ2 GR2; do
    export i=$i
sbatch -A snic2019-3-557 -t 3:00:00 -p core -n 10 \
-o ${i}_alignment.out -e ${i}_alignment.err<<'EOF'
#!/usr/bin/env bash

module load bioinfo-tools
module load bowtie2

time bowtie2 \
    -x bin_quantification/all_bins \
    -1 ../Metagenome/${i}_R1.fastq.gz \
    -2 ../Metagenome/${i}_R2.fastq.gz \
    -S bin_quantification/${i}.sam \
    --end-to-end --no-unal

EOF
done
```

Convert to BAM

```bash
module load samtools

# bam should be unsorted
samtools view -b bin_quantification/GQ2.sam > bin_quantification/GQ2.bam
samtools view -b bin_quantification/GR2.sam > bin_quantification/GR2.bam
```

Quantify with salmon

```bash
salmon quant -t bin_quantification/all_bins.fa \
    -l A \
    -a bin_quantification/GQ2.bam \
    -o bin_quantification/GQ2_salmon \
    --meta \
    -p 20
```

## Use all contigs

### Map reads with bowtie2

```bash
cd /home/domeni/theCellarDoor/SLUBI_Project_TE_Biogas_snic2019-30-23/SLUBI_analyses
mkdir -p bin_quantification_with_unbinned
```

Collect all sequences (binned/unbinned).
Contig names as _eg_ `AS_TE.75__k141_123` if they are binned; otherwise `Unbinned__k141_234`.

```python
import glob, os
from Bio import SeqIO

# create dictionary contigID: bin
bin_files = glob.glob("../Binning/Bins/AS*.fa")
outhandle = open("bin_quantification_with_unbinned/all_bins_with_unbinned.fa", 'w')
bin_dict = {}
for bin in bin_files:
    bin_name = os.path.split(bin)[1].replace(".fa", "")
    bin_handle = SeqIO.parse(bin, 'fasta')
    for contig in bin_handle:
        bin_dict[contig.id] = bin_name

# iterate over contigs
contig_seqs = SeqIO.parse("../Assembly/AS_TE.contigs.fa", 'fasta')
for contig in contig_seqs:
    if contig.id in bin_dict:
        id_out = "{}__{}".format(bin_name, contig.id)
        _ = outhandle.write(">{}\n{}\n".format(id_out, str(contig.seq)))
    else:
        id_out = "{}__{}".format("Unbinned", contig.id)
        _ = outhandle.write(">{}\n{}\n".format(id_out, str(contig.seq)))


# # cat all bin sequences, but prepend bin name to contig names
# bin_files = glob.glob("../Binning/Bins/AS*.fa")
# outhandle = open("bin_quantification/all_bins.fa", 'w')
# for bin in bin_files:
#     print(bin)
#     bin_name = os.path.split(bin)[1].replace(".fa", "")
#     bin_handle = SeqIO.parse(bin, 'fasta')
#     for contig in bin_handle:
#         id_out = "{}__{}".format(bin_name, contig.id)
#         _ = outhandle.write(">{}\n{}\n".format(id_out, str(contig.seq)))

outhandle.close()
```

Create bowtie2 index

```bash
sbatch -A snic2019-3-557 -t 3:00:00 -p core -n 10 \
-o bin_quantification_with_unbinned/all_bins_with_unbinned.bowtie2_build.out<<'EOF'
#!/usr/bin/env bash

module load bioinfo-tools
module load bowtie2

bowtie2-build \
    --threads 20 \
    bin_quantification_with_unbinned/all_bins_with_unbinned.fa \
    bin_quantification_with_unbinned/all_bins_with_unbinned

EOF
```

Map reads

```bash
for i in GQ2 GR2; do
    export i=$i
sbatch -A snic2019-3-557 -t 3:00:00 -p core -n 10 \
-o bin_quantification_with_unbinned/${i}_alignment.out -e bin_quantification_with_unbinned/${i}_alignment.err<<'EOF'
#!/usr/bin/env bash

module load bioinfo-tools
module load bowtie2

time bowtie2 \
    -x bin_quantification_with_unbinned/all_bins_with_unbinned \
    -1 ../Metagenome/${i}_R1.fastq.gz \
    -2 ../Metagenome/${i}_R2.fastq.gz \
    -S bin_quantification_with_unbinned/${i}.sam \
    --end-to-end --no-unal

EOF
done
```

Convert to BAM

```bash
module load samtools

# bam should be unsorted
samtools view -b bin_quantification_with_unbinned/GQ2.sam > bin_quantification_with_unbinned/GQ2.bam
samtools view -b bin_quantification_with_unbinned/GR2.sam > bin_quantification_with_unbinned/GR2.bam
```

Quantify with salmon

```bash
salmon quant -t bin_quantification_with_unbinned/all_bins_with_unbinned.fa \
    -l A \
    -a bin_quantification_with_unbinned/GQ2.bam \
    -o bin_quantification_with_unbinned/GQ2_salmon \
    --meta \
    -p 20

salmon quant -t bin_quantification_with_unbinned/all_bins_with_unbinned.fa \
    -l A \
    -a bin_quantification_with_unbinned/GR2.bam \
    -o bin_quantification_with_unbinned/GR2_salmon \
    --meta \
    -p 20

```

Map metatranscriptome reads on these new bowtie2 dbs.

```bash
mkdir -p bin_quantification_with_unbinned_metat
# 
# for i in GQ2 GR2; do
#     export i=$i
sbatch -A snic2019-3-557 -t 3:00:00 -p core -n 20 \
--array=2-$(wc -l < ../Metatranscriptome/metadata.tsv) \
-o bin_quantification_with_unbinned_metat/alignment_%a.out -e bin_quantification_with_unbinned_metat/alignment_%a.err<<'EOF'
#!/usr/bin/env bash

module load bioinfo-tools
module load bowtie2
module load samtools

#let LINE=${SLURM_ARRAY_TASK_ID}+1
LINE=${SLURM_ARRAY_TASK_ID}
SAMPLE=$(cat ../Metatranscriptome/metadata.tsv | sed -n ${LINE}p | awk '{print $1}')

time bowtie2 \
    -x bin_quantification_with_unbinned/all_bins_with_unbinned \
    -1 ../Metatranscriptome/${SAMPLE}_R1_trimmed.fastq.gz \
    -2 ../Metatranscriptome/${SAMPLE}_R2_trimmed.fastq.gz \
    -S bin_quantification_with_unbinned_metat/${SAMPLE}.sam \
    --end-to-end --no-unal -p 20

samtools view -b bin_quantification_with_unbinned_metat/${SAMPLE}.sam > bin_quantification_with_unbinned_metat/${SAMPLE}.bam

EOF
```

Quantify with salmon

```bash
salmon quant -t bin_quantification_with_unbinned/all_bins_with_unbinned.fa \
    -l A \
    -a bin_quantification_with_unbinned/GQ2.bam \
    -o bin_quantification_with_unbinned/GQ2_salmon \
    --meta \
    -p 20

salmon quant -t bin_quantification_with_unbinned/all_bins_with_unbinned.fa \
    -l A \
    -a bin_quantification_with_unbinned_metat/GQ2-22-10-1.bam \
    -o bin_quantification_with_unbinned_metat/GQ2-22-10-1 \
    --meta \
    -p 20

```

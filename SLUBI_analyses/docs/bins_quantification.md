# Bin quantification in metagenomic samples

## Salmon

Installed v1.3 as conda env.

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

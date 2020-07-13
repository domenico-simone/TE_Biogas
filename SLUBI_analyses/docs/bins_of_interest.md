```bash
cd /home/domeni/theCellarDoor/SLUBI_Project_TE_Biogas_snic2019-30-23/SLUBI_analyses
mkdir -p bins_of_interest/fastANI

for i in $(ls bins_of_interest/*.fa)
do
    echo $i
    # bins_of_interest/bin.23.fa
    outname_i=$(basename $i)
    for j in $(ls bins_fasta/AS_TE.*.fa)
    do
        echo $j
        outname_j=$(basename $j)
        outname=${outname_i}_${outname_j}
        # fastANI -q bins_of_interest/bin.23.fa -r bins_fasta/AS_TE.15.fa -o bins_of_interest/bin.23_AS_TE.15.out
        fastANI -q $i -r $j -o bins_of_interest/fastANI/${outname}.out -t 20
    done
done
```

Relevant results?

```bash
cat bins_of_interest/fastANI/*.out
# bins_of_interest/bin.23.fa	bins_fasta/AS_TE.34.fa	98.7925	447	547
# bins_of_interest/bin.23.fa	bins_fasta/AS_TE.63.fa	78.8538	278	547
# bins_of_interest/bin.28.fa	bins_fasta/AS_TE.22.fa	98.1476	354	561
# bins_of_interest/bin.28.fa	bins_fasta/AS_TE.42.fa	84.6216	71	561
# bins_of_interest/bin.54.fa	bins_fasta/AS_TE.22.fa	86.2978	298	603
# bins_of_interest/bin.54.fa	bins_fasta/AS_TE.42.fa	95.6625	94	603
# bins_of_interest/bin.64.fa	bins_fasta/AS_TE.34.fa	83.2743	326	527
# bins_of_interest/bin.64.fa	bins_fasta/AS_TE.63.fa	79.4681	283	527
```

Re-run for those with relevant results to get the visual output

```bash
i=bins_of_interest/bin.23.fa
j=bins_fasta/AS_TE.34.fa
outname_i=$(basename $i)
outname_j=$(basename $j)
outname=${outname_i}_${outname_j}
fastANI \
-q $i \
-r $j \
-o bins_of_interest/fastANI/${outname}.out \
-t 20 \
--visualize

i=bins_of_interest/bin.28.fa
j=bins_fasta/AS_TE.22.fa
outname_i=$(basename $i)
outname_j=$(basename $j)
outname=${outname_i}_${outname_j}
fastANI \
-q $i \
-r $j \
-o bins_of_interest/fastANI/${outname}.out \
-t 20 \
--visualize

i=bins_of_interest/bin.54.fa
j=bins_fasta/AS_TE.42.fa
outname_i=$(basename $i)
outname_j=$(basename $j)
outname=${outname_i}_${outname_j}
fastANI \
-q $i \
-r $j \
-o bins_of_interest/fastANI/${outname}.out \
-t 20 \
--visualize

```




```python
A = open("fastANI/test1.out.visual", 'r')
c = 0
for i in A:
    i = i.split()
    c += int(i[7]) - int(i[6]) + 1

    >>> c
    1341000
    >>> from Bio import SeqIO
    >>> Q = SeqIO.parse("bin.23.fa", 'fasta')
    >>> l
    1995130
    >>> Q = SeqIO.parse("../bins_fasta/AS_TE.34.fa", 'fasta')
    2310158

```
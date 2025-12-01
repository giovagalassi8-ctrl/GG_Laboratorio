 # Polishing

### Mapping short and long reads

#### Mapping short reads

```bash
#Short reads
minimap2 -ax sr --MD -t 6 Anoste_raw.fasta SRR11672503_1_paired.fastq SRR11672503_2_paired.fastq > Anoste_raw_sr.sam
samtools view -Sb Anoste_raw_sr.sam > Anoste_raw.bam
rm Anoste_raw_sr.sam
samtools sort -@6 -o Anoste_raw_sr_sorted.bam Anoste_raw_sr.bam
samtools index Anoste_raw_sr.bam
rm Anoste_raw_sr.bam
```

#### Mapping long reads

```bash
#Long reads
minimap2 -ax map-pb --MD -t 6 Anoste_raw.fasta SRR11672506.fastq.gz > Anoste_raw_lr.sam
samtools view -Sb Anoste_raw_lr.sam > Anoste_raw_lr.bam
rm Anoste_raw_lr.sam
samtools sort -@6 -o Anoste_raw_lr_sorted.bam Anoste_raw_lr.bam
samtools index Anoste_raw_lr_sorted.bam
rm Anoste_raw_lr.bam
```

### Hypo
#### Mean coverage of short reads

```bash
mosdepth -n --fast-mode --by 500 Anoste_raw_sr Anoste_raw_sr_sorted.bam
```

#### Hypo
```bash
echo -e "$R1\n$R2" > Sr.Path
hypo -d Anoste_raw.fasta -r @Sr.path -s <APPROXIMATE_GENOMESIZE> -c <SHORT_READSCOVERAGE> -b <SORTED_BAM_SR> -B <SORTED_BAM_PB> -t <NUMBER_THREADS>
```

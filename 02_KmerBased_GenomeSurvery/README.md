# Kmer - Based Genome Survery

## Reads quality check

```bash
#[assembly]
fastqc SRR11672503_1.fastq.gz SRR11672503_2.fastq.gz
```
-----

## Trimming

```bash
#[assembly]
trimmomatic PE -threads 8 -phred33 SRR11672503_1.fastq.gz SRR11672503_2.fastq.gz SRR11672503_1_paired.fastq SRR11672503_1_unpaired.fastq SRR11672503_2_paired.fastq SRR11672503_2_unpaired.fastq ILLUMINACLIP:/opt/miniforge3/envs/assembly/share/trimmomatic-0.40-0/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 2> stats_trimmomatic.log
```
-----

## Compute *k-mer* frequency

```bash
#[kat]
kat hist -t 8 -m 27 -o Anoste SRR11672503_1_paired_fastqc.html SRR11672503_2_paired_fastqc.html
```

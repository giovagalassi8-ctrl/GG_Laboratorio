# Raw assembly

### Contig-level assembly 
#### Assembling reads

```bash
wtdbg2 -x rs -g 227054799 -t 8 -i SRR11672506.fastq.gz -fo Anoste_raw
```
-----

#### Final consensus in fasta format

```bash
wtpoa-cns -t 7 -i Anoste_raw.ctg.lay.gz -fo Anoste_raw
```
-----

### Quality check

#### N50
```bash
assembly-stats Anoste_raw.fasta > Anoste_raw.stats
```

#### Busco

```bash
export NUMEXPR_MAX_THREADS= 80
busco -m geno -l $BUSCO/culicidae_odb12 -c 6 -o Anoste_raw.busco -i Anoste_raw.fasta
```

#### Spectra-cn (KAT)

```bash
kat comp -t 8 -o Anoste_pol 'SRR11672503_1_paired.fastq SRR11672503_2_fastq' Anoste_pol.fasta
```
-----

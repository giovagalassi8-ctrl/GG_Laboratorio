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


-----

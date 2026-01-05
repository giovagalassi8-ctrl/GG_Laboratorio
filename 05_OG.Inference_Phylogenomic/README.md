# Ortholgy inference and Phylogenomic

This represents the comparative analysis proper, where the functional portion of the genome and selective pressures are examined to link specific evolutionary patterns to corresponding phenotypic patterns.
The analysis aims to identify genes that exhibit similar selective pressures across species sharing the same phenotype.

## Dataset creation

To perform the comparison, it is essential to have a comprehensive genomic dataset, in this case provided in the [mosquitoes.tsv](../00_practice/00_dataset/mosquitoes.tsv) file.

To retrieve the genomes selected for the comparative analysis, a dedicated script was developed, as shown below:

```bash
#!/bin/bash

#This script downloads genomes and relative gff from NCBI datasets creating ready to use folders

AN2name=$1

mkdir 00_genome
mkdir 01_gff

while IFS=$'\t' read -r AN sname ID; do
	echo $AN
	#download specifying the name of the folder
	datasets download genome accession "$AN" --filename "$ID".zip --include genome,gff3
	#unzip specifying the name of the folder
	unzip "$ID".zip -d "$ID"
	#rename the two file of interest
	mv "$ID"/ncbi_dataset/data/"$AN"/*.fna 00_genome/"$ID".fna
	mv "$ID"/ncbi_dataset/data/"$AN"/*.gff 01_gff/"$ID".gff
	#delete the folder
	rm -r "$ID"/
done < "$AN2name"
```

## Longest isoform and translation

Once the genomes have been downloaded, the raw data consists of nucleotide sequences. However, in comparative genomics, analyses are typically conducted at the amino acid level rather than using raw nucleotides. To perform this conversion, the AGAT toolkit is utilized. Specifically, once the script results are obtained, a for-cicle can be executed to process the files through this command.

```bash
#[sequence]
for gff in *.gff; do agat_sp_keep_longest_isoform.pl --gff "$gff" -o ${gff/.gff/_longest.gff}; done
```

This initial script is designed to select only the longest isoform for each gene. Including all alternative isoforms would artificially inflate the proteome size; by retaining only the longest sequence, we ensure that the maximum number of exons is represented for each locus.

```bash
#[sequence]
for gff in *_longest.gff; do agat_sp_extract_sequences.pl -g "$gff" -f ../00_genome/${gff/_longest.gff/.fna} -t cds -p --cfs --output ../02_raw_proteoms/${gff/_longest.gff/.faa}; done
```

This second script is designed to extract CDS (Coding Sequence) features and translate them into proteins using the standard eukaryotic genetic code.

## Removal of pseudogenes

Despite being reference-grade sequences, these genomes may still contain pseudogenes (defined here as any sequence containing internal premature stop codons). Consequently, these must be removed. The following script was employed for this purpose:

```bash
#!/bin/bash/

# list each plausible pseudogene present. 

mkdir raw_proteomes
mv *.faa raw_proteomes/

#make proteome one line
cd raw_proteomes
for proteome in *.faa; do
	awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' "$proteome" > ../${proteome/.faa}".faa"
done
cd ..

#Extract all pseudogene names
mkdir 00_pseudogene_name

for proteome in *.faa; do
	species=$(basename -s .faa "$proteome")
	grep -B1 '*' "$proteome" | grep ">" >> 00_pseudogene_name/"$species"_pseudogenes_name.txt
done

#removes sequences identified as pseudogenes

for pseudo_file in 00_pseudogene_name/*_pseudogenes_name.txt; do
	species=$(basename -s _pseudogenes_name.txt "$pseudo_file")
	while IFS=$'\t' read -r header; do
		sed -E -i "/${header}/{N;d;}" "$species".faa # N option loads the next line found after the pattern and put it into pattern space too; d delete the pattern space
	done < "$pseudo_file" 
done

mv 00_pseudogene_name ../00_genome
```
### Header modification

At this stage, the FASTA headers must be modified to ensure they are clear and functional for downstream analyses. Specifically, the headers were truncated to retain only the species identifier and the unique sequence ID. This step is indispensable for the subsequent stages of the analysis.

```bash
for prot in *.faa; do ID=$(basename -s .faa "$prot"); sed -i.old -E "s/>(rna-XM_[0-9]+\.[0-9]) (gene=gene-.[^ ]+) name=(.[^ ]+) .+$/ >${ID}\|\3/" "$prot"; done
```

## Orthofinder

OrthoFinder was employed to identify orthologous groups. The software identifies orthologs through a pipeline that searches for maximum sequence similarity to infer orthology. It utilizes a tree-based method, defining orthology through speciation events; consequently, it distinguishes between duplication events, which create paralogs, and speciation events, which result in orthologs, based on the reconstructed gene trees.

```bash
orthofinder -t 8 -a 8 -f .
```

The software generates a comprehensive output organized into several directories, each containing specific types of evolutionary and genomic information.

## Paralog filtering

To refine the dataset, it was necessary to filter the orthogroups to remove paralogs. For this purpose, we utilized DISCO, which decomposes complex orthogroups into distinct, orthology-consistent sub-clusters. DISCO script requires a specific header syntax to properly work. Luckily, the supported syntax is the one we already implemented and used so far (SPECIES|SEQUENCE_ID).

```bash
while IFS=' ' read -r OG tree; do python3 ~/GG_Laboratorio/99_scripts/disco.py -i <(echo "$tree") -o ../../../01_disco/${OG/:/}.nwk -d "|" -m 4 --remove_in_paralogs --keep-labels --verbose >> ../../../01_disco/disco.log; done < <(sed -E 's/[A-Z][a-z]{5}_//g; s/\)n[0-9]*+/\)/g' Resolved_Gene_Trees.txt)
```

> The disco.py script is available at the following [link](https://github.com/jsdoublel/DISCO/blob/master/disco.py)

To optimize computational resources and storage, empty DISCO output files (zero-size files) were identified and removed from the dataset.

```bash
find . -size 0 -print > empty_disco.txt
find . -size 0 -delete
```

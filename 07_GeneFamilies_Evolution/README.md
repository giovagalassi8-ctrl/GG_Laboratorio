# Gene Families Evolution

## CAFE 

CAFE is one of the most widely used software tools based on gene and species trees. Its primary purpose is to analyze changes in gene family sizes, accounting for phylogenetic history to provide a statistical foundation for evolutionary inferences.

The program employs a birth-and-death process to model gene gain and loss across a user-specified phylogenetic tree. The distribution of family sizes generated under this model serves as a baseline to assess the statistical significance of observed size differences among taxa. In brief, CAFE estimates the ancestral states of gene families by comparing the tip states (extant family sizes) across the phylogeny.

----

### Data Preparation

CAFE requires two main inputs: the ultrametric tree (generated in the previous section) and a gene count table. The table must be formatted as follows:
NONE Orthogroup ID ID ID ID
NONE OG0000000  93 0  5  9 

OrthoFinder generates a similar table, located in the results folder as Orthogroups.GeneCount.tsv. However, this file requires preprocessing to match CAFE's input format:

-The first column required by CAFE is missing.

-There is an extra column that is unnecessary and must be removed.

This file was processed to meet CAFE's input requirements, ensuring the correct column structure and formatting.

```bash
sed $'s/^/NONE\t/g' Orthogroups.GeneCount.tsv | rev Orthogroups.GeneCount.tsv  | cut -f 2- | rev > ../../../../07_GeneFamilies_Evolution/GeneCount_CAFE.tsv 
```
Finally, please verify that your phylogenetic tree is formatted as Newick rather than Nexus before proceeding with the analysis.

----

### Generating the Error Model

While the program can immediately infer one or more birth-death parameters (lambda), it is standard practice to first estimate a general error model. 
This statistical model describes the probability that the observed gene count in a family differs from the true count due to technical artifacts (such as genome assembly or annotation errors). The statistics derived from this model will be incorporated into subsequent analysis steps to ensure more accurate evolutionary inferences.

```bash
cafe5 -i GenCount_CAFE.tsv -t timetree.nwk -o Error_model -e
```

----

### Single-Lambda Analysis
Then we can run the program specifying the error model just inferred. An initial CAFE5 analysis was performed using a single $\lambda$ (lambda) parameter for the entire tree. This model assumes a uniform turnover rate, implying that the tendency for gene families to expand or contract remains constant across all species.

```bash
for k in {1..5}; do for n in {1..10}; do mkdir -p 00_1L/${k}K/${n}N; cafe5 -i GenCount_ -t timetree.nwk -o 00_1L/${k}K/${n}N -e./Error_model/Base_error_model.txt -k ${k}; done; done
```
----

## Two-Lambda Analysis

To allow for rate heterogeneity across the phylogeny, we performed a two-lambda analysis. This model assigns distinct turnover rates to specific lineages. To compute multiple $\lambda$ rates, you must define the number of distinct rate categories and assign them to specific species or clades.These assignments are specified in a separate tree file (distinct from the main ultrametric tree) passed via the -y flag. In this file, integers are used to map specific branches to their corresponding $\lambda$ category.

```bash
(Anofun:1,(Anoste:1,((Culqui:2,Culpip:2):1,(Aedalb:1,Aedaeg:1):1):1):1);
```

The execution script for the two-lambda analysis is provided below:

```bash
for k in {1..5}; do for n in {1..10}; do mkdir -p 00_2L/${k}K/${n}N; cafe5 -i GeneCount_ -t timetree.nwk -o 00_2L/${k}K/${n}N -y timetree_2Lambda.nwk -e./Error_model/Base_error_model.txt -k ${k}; done; done
```

----

## Evolutionary Rate Analysis

Once the gene family dynamics have been established, the next step is to characterize the specific genes and their evolutionary history. Various tools, such as CodeML and HyPhy, are available for this purpose.Both programs utilize the dN/dS ratio ($\omega$), a metric of natural selection acting on protein-coding sequences. The primary goal of this analysis is to identify signatures of positive selection that may be restricted to specific lineages (branches) or specific amino acid sites.

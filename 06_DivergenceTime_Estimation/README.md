# Divergence Time Estimation

We proceeded with the dating of the species tree obtained in the previous step to generate an ultrametric tree. The tree was initially calibrated using [iTol](https://itol.embl.de/), specifically by selecting the branch containing the outgroup.

At this stage, it is necessary to date the nodes to obtain a functional time-tree. This step was performed using the LSD2 (Least Squares Dating) algorithm in IQ-TREE, although it has since been superseded by more advanced methods. To complete the calibration, we utilized [TimeTree.org](https://timetree.org/), a comprehensive resource containing divergence times and paleo-ecological data for all known organisms.

On TimeTree.org, the [mosquitoes.tsv](../00_practice/00_dataset/mosquitoes.tsv) dataset can be uploaded; however, it is crucial to ensure that the first column is modified so that underscores are replaced by spaces. If underscores are retained, the database will fail to recognize the taxonomic names.

# Profiling_metagenomes_and_mags

This repository encapsulates a workflow to profile metagenomes and their associated MAGs.

More specifically, the workflow involves:
1. [Build a gene catalogue from metagenomic assemblies](https://github.com/tpriest0/Profiling_metagenomes_and_mags/wiki/01-Constructing-gene-catalogue)
2. Profile the gene catalogue: estimate the coverage per haploid genome of each gene in each metagenome (and metatranscriptomes if available - in which case expression is also calculated)
3. Profile protein functions: estimate the coverage per haploid genome of each protein function in each metagenome (and metatranscriptomes if available - in which case expression is also calculated)
4. Profile metagenome-assembled genomes and their protein functional content: estimate the proportion of haploid genomes in each metagenome (and metatranscriptomes if available - in which case expression is also calculated) attributed to each MAG. In addition, estimate the coverage per haploid genome of each protein function in each MAG and the proportion of that function in relation to the whole community (also computed on metatranscriptomes if available).


# DUGMO
Detection of unknown genetically modified organisms

## Introduction

DUGMO is a bioinformatics pipeline for the detection of unknown GMO based on paired end data from Illumina sequencing. This tool performs a cleaning of high throughput sequencing data to identify coding sequences (CDS) belonging to the wild genome from potential GMO coding sequences (pgmCDS). These CDS are aligned in two successive blasts against the host pangenome with relevant tuned parameters. Then, Bray-Curtis distances are calculated between the CDSs and the pgmCDSs based on the difference of genomic vocabulary between these two sets. Finally, a machine learning method, random forest, is carried out to target true GMO CDS(s) with six different variables among which Bray-Curtis distances and GC content.

## Installation with conda 

Create a conda environnement 
```
conda create -n DUGMO snakemake
conda activate DUGMO
```
## Quick start

Launch the pipeline snakemake
```
snakemake -rp --use-conda
```
## Paramet

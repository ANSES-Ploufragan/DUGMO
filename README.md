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
## Params in config.json

Params name | Description
------------|------------
`R1 :` | Input file Read 1 in fq.gz format
`R2 :` | Input file Read 2 in fq.gz format
`RefGenome` | Input file with 
`PangenomeCDS` | Input file with 
`PangenomeWhole` | Input file with 
`outputFolder` | Output folder for all results
`threads` | Number of threads to use (default: 10)

## Output files

Foldername | Description
-----------|------------
`Cleaning_pipeline` | Cleaning pipeline with Shovill, Trimmomatic, BWA, Bowtie2, Blast and Prokka
`HostGenome` | Host genome files
`GMO` | GMO potential files
`TrainingGMO` | Databank files
`MachineLearning` | Results at the final step


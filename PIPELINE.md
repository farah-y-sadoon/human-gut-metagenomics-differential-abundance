# Diet-driven differences in the human gut microbiome: A bacterial community analysis

## 0. Environment Setup and Data Aquisition
Create metadata CSV (manual)
```bash
nano data/metadata/metadata.csv
sample_id,diet,srr,location
1,omnivore,SRR8146935,Turin
2,omnivore,SRR8146936,Turin
3,omnivore,SRR8146938,Turin
4,vegan,SRR8146944,Turin
5,vegan,SRR8146951,Bari
6,vegan,SRR8146952,Bari
```
Download database
```bash
# Submit job to download database
sbatch scripts/01_download_database.sh
```
Download reads
```bash
# Submit job to download reads
sbatch scripts/02_download_reads.sh
```
Convert SRA files to FASTQ format
```bash
# Submit job to convert SRA files to FASTQ format
sbatch scripts/03_convert_sra_fastq.sh
```
## 1. Quality Control
Assess quality of raw reads with FastQC and MultiQC
```bash
# Submit job to run FastQC
sbatch scripts/04_fastqc_raw_reads.sh

# Create mamba environment for MultiQC
mamba install -n multiqc -c bioconda multiqc

# Submit job to run MultiQC
sbatch scripts/05_multiqc_raw_reads_report.sh
```
Trim adapter sequences using Fastp
```bash
# Submit job to trim adapters with Fastp
sbatch scripts/06_trim_adapters.sh
```

## 2. Taxonomic Classification with Kraken2 / KrakenUniq

## 3. Abundance Restimation with Bracken

## 4. Data Analysis and Exploration in R
Rarefaction to observe number of taxa
Visualize relative abundance
Alpha diversity - 2 measures
Beta diversity - 2 measures

## Differential Abundance Analysis in R with ANCOM-BC2 / Aldex3
Heatmap
Volcano plots or Bar plots
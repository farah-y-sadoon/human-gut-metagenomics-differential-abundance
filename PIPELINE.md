# Diet-driven differences in the human gut microbiome: A bacterial community analysis

## 0. Environment Setup and Data Aquisition
A folder was created in the scratch directory of the Digital Research Alliance Canada's Fir cluster for this analysis.

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
mamba create -n multiqc -c bioconda multiqc -y

# Submit job to run MultiQC
sbatch scripts/05_multiqc_raw_reads_report.sh
```
Trim adapter sequences using Fastp
```bash
# Submit job to trim adapters with Fastp
# Note: Fastp generates its own QC report for each sample in results/fastp/
sbatch scripts/06_trim_adapters.sh
```
## 2. Taxonomic Classification with Kraken2
Classify reads with Kraken2
```bash
# Submit job to perform taxonomic classification with Kraken2
sbatch scripts/07_kraken_classify.sh
```
## 3. Abundance Restimation with Bracken
Re-estimate abundance with Bracken
```bash
# Submit job to perform abundance re-estimation with Bracken
sbatch scripts/08_bracken_estimate.sh

# Create mamba kraken-bion environment
mamba create -n kraken-biom -c bioconda kraken-biom -y

# Run script to use kraken-biom to convert the Bracken report output to BIOM format
./scripts/09_kraken_biom_convert.sh
``` 
## 4. Diversity and Differential Abundance Analysis in R with Aldex3
```bash
# Run R script for diversity and differential abundance analysis
Rscript scripts/10_differential_abundance_analysis.R
```
#!/bin/bash
# Script to generate multiqc report for raw reads

#SBATCH --job-name=multiqc_raw
#SBATCH --time=01:00:00
#SBATCH --account=def-cottenie
#SBATCH --cpus-per-task=1
#SBATCH --error=logs/%j_error.txt
#SBATCH --output=logs/%j_output.txt
#SBATCH --mem=4G

# Activate mamba environment
source ~/.bashrc
mamba activate multiqc

# Define variables
INPUT_DIR="./results/qc_results/raw_reads"
OUTPUT_DIR="./results/qc_results/raw_reads/multiqc/"

# Make output directory
mkdir -p "$OUTPUT_DIR"

# Run fastqc for all files in data/raw directory
multiqc "$INPUT_DIR" -o "$OUTPUT_DIR"

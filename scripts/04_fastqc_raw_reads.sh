#!/bin/bash
# Script to check read quality for samples with FastQC 

#SBATCH --job-name=fastqc_raw
#SBATCH --time=02:00:00
#SBATCH --account=def-cottenie
#SBATCH --cpus-per-task=12
#SBATCH --error=logs/%j_error.txt
#SBATCH --output=logs/%j_output.txt
#SBATCH --mem=24G

# Load module(s)
module load StdEnv/2023 fastqc/0.12.1

# Define variables
INPUT_DIR="./data/raw"
OUTPUT_DIR="./results/qc_results/raw_reads"

# Make output directory
mkdir -p "$OUTPUT_DIR"

# Run fastqc for all files in data/raw directory
fastqc "$INPUT_DIR"/*.fastq -o "$OUTPUT_DIR" -t 12

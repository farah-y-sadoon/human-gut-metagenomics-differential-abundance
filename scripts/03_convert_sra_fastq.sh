#!/bin/bash
# Script to convert SRA files to FASTQ files

#SBATCH --job-name=sra_to_fastq
#SBATCH --time=04:00:00
#SBATCH --account=def-cottenie
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=logs/fastq_conv_%j.txt

module load StdEnv/2023 sra-toolkit/3.0.9

# Set directories
WORKING_DIR="./data/raw"
TEMP_DIR="/scratch/fsadoon/sra_cache/tmp" # to deal with temporary files created by fasterq-dump

# Create temporary directory
mkdir -p "$TEMP_DIR"

for srr_dir in "$WORKING_DIR"/*/; do
    # Get the accession
    srr=$(basename "$srr_dir")
    
    echo "Processing $srr in $WORKING_DIR..."
    
    fasterq-dump "$srr_dir/${srr}.sra" \
        -O "$WORKING_DIR" \
        --split-3 \
        --temp "$TEMP_DIR" \
        --threads 8
done

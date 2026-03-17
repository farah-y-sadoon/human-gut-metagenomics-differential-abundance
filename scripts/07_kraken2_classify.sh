#!/bin/bash
# Script to classify metagenome with Kraken2

#SBATCH --job-name=kraken_metagenomics
#SBATCH --time=8:00:00
#SBATCH --account=def-cottenie
#SBATCH --cpus-per-task=32
#SBATCH --mem=400G
#SBATCH --error=logs/%j_kraken.err
#SBATCH --output=logs/%j_kraken.out

# Load module(s) 
module load kraken2/2.1.6

# Define path variables
DB_PATH="/scratch/fsadoon/kraken2_db/"
INPUT_DIR="./data/filtered"
OUTPUT_DIR="./results/kraken2"

# Make output directory
mkdir -p "$OUTPUT_DIR"

# Iterate through samples

for R1 in "$INPUT_DIR"/*_1.clean.fastq.gz; do
    SAMPLE=$(basename "$R1" _1.clean.fastq.gz)
    R2="${INPUT_DIR}/${SAMPLE}_2.clean.fastq.gz"

    echo "Starting classification for sample $SAMPLE..."

    # Run Kraken2 on each sample
    kraken2 --db "$DB_PATH" \
        --confidence 0.15 \
        --threads 32 \
        --paired \
        --gzip-compressed \
        --output "${OUTPUT_DIR}/${SAMPLE}.kraken2" \
        --report "${OUTPUT_DIR}/${SAMPLE}_report.tsv" \
        "$R1" "$R2"

    echo "Finished classification for sample $SAMPLE."
done
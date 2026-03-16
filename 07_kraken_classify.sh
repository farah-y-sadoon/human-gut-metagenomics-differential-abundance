#!/bin/bash
# Script to classify metagenome with KrakenUniq

#SBATCH --job-name=kraken_metagenomics
#SBATCH --time=08:00:00
#SBATCH --account=def-cottenie
#SBATCH --cpus-per-task=32
#SBATCH --mem=400G
#SBATCH --error=logs/%j_kraken.err
#SBATCH --output=logs/%j_kraken.out

# Load the module
module load StdEnv/2023 krakenuniq/1.0.4

# Define path variables
DB_PATH="../kraken_uniq_db"
INPUT_DIR="./data/filtered"
OUTPUT_DIR="./results/krakenuniq"

# Make output directory
mkdir -p "$OUTPUT_DIR"

# Iterate through samples and run KrakenUniq for taxonomic classification
for R1 in "$INPUT_DIR"/*_1.clean.fastq.gz; do
    SAMPLE=$(basename "$R1" _1.clean.fastq.gz)
    R2="${INPUT_DIR}/${SAMPLE}_2.clean.fastq.gz"

    echo "Starting classification for sample $SAMPLE..."

    krakenuniq --db "$DB_PATH" \
        --threads 32 \
        --preload \
        --paired "$R1" "$R2" \
        --output "${OUTPUT_DIR}/${SAMPLE}.kraken" \
        --report-file "${OUTPUT_DIR}/${SAMPLE}_report.tsv" \
        --gzip-compressed

    echo "Finished classification for sample $SAMPLE."
done
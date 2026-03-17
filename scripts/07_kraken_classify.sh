#!/bin/bash
# Script to classify metagenome with KrakenUniq

#SBATCH --job-name=kraken_metagenomics
#SBATCH --time=10:00:00
#SBATCH --account=def-cottenie
#SBATCH --cpus-per-task=32
#SBATCH --mem=450G
#SBATCH --error=logs/%j_kraken.err
#SBATCH --output=logs/%j_kraken.out

# Define path variables
ROOT_DIR="/scratch/fsadoon/"
SIF_PATH="${ROOT_DIR}/human-gut-metagenomics-differential-abundance/tools/krakenuniq.sif"
DB_PATH="${ROOT_DIR}/krakenuniq_db"
INPUT_DIR="${ROOT_DIR}/human-gut-metagenomics-differential-abundance/data/filtered"
OUTPUT_DIR="${ROOT_DIR}/human-gut-metagenomics-differential-abundance/results/krakenuniq"

# Make output directory
mkdir -p "$OUTPUT_DIR"

# Iterate through samples
for R1 in "$INPUT_DIR"/*_1.clean.fastq.gz; do
    SAMPLE=$(basename "$R1" _1.clean.fastq.gz)
    R2="${INPUT_DIR}/${SAMPLE}_2.clean.fastq.gz"

    echo "Starting classification for sample $SAMPLE..."

    # Run KrakenUniq on each sample from apptainer
    apptainer exec --bind /scratch:/scratch "$SIF_PATH" \
    krakenuniq --db "$DB_PATH" \
        --threads 32 \
        --preload \
        --paired "$R1" "$R2" \
        --output "${OUTPUT_DIR}/${SAMPLE}.kraken" \
        --report-file "${OUTPUT_DIR}/${SAMPLE}_report.tsv"

    echo "Finished classification for sample $SAMPLE."
done

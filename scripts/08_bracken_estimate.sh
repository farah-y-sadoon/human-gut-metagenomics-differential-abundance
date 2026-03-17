#!/bin/bash
# Script to perform abundance estimation

#SBATCH --job-name=bracken_estimate
#SBATCH --time=02:00:00
#SBATCH --account=def-cottenie
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --error=logs/%j_bracken.err
#SBATCH --output=logs/%j_bracken.out

# Load module(s)
module load StdEnv/2023 bracken/3.0

# Define path variables
DB_PATH="/scratch/fsadoon/kraken2_db/"
INPUT_DIR="./results/kraken2"
OUTPUT_DIR="./results/bracken"

# Make output directory
mkdir -p "$OUTPUT_DIR"

# Iterate through Kraken2 reports and re-estimate abundance
for REPORT in "$INPUT_DIR"/*_report.tsv; do
    SAMPLE=$(basename "$REPORT" _report.tsv)

    echo "Starting $SAMPLE..."
    bracken -d "$DB_PATH" \
            -i "$REPORT" \
            -o "${OUTPUT_DIR}/${SAMPLE}_species.bracken" \
            -r 150 \
            -l S \
            -t 10
    echo "Finished $SAMPLE."
done

echo "Abundance estimation complete."
ls -lh $OUTPUT_DIR

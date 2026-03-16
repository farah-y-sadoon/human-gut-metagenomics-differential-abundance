#!/bin/bash
# Script to trim adapters with Fastp

# CHANGE THIS TO MATCH JOB SPECS
#SBATCH --job-name=trim_adapters
#SBATCH --time=01:00:00
#SBATCH --account=def-cottenie
#SBATCH --cpus-per-task=8
#SBATCH --error=logs/%j_error.txt
#SBATCH --output=logs/%j_output.txt
#SBATCH --mem=16G

# Load module(s)
module load fastp/1.0.1 python/3.14.2

# Define variables
INPUT_DIR="./data/raw"
OUTPUT_DIR="./data/filtered"
REPORT_DIR="./results/filtered_reads"
TOOLS_DIR="./tools"
URL="https://raw.githubusercontent.com/OpenGene/fastp/master/parallel.py"

# Make output directory
mkdir -p "$OUTPUT_DIR" "$REPORT_DIR" "$TOOLS_DIR"

# Download parallel.py script to run fastp all at once
echo "Trimming and filtering reads..."
wget "$URL" -O "$TOOLS_DIR"/parallel.py

python3 "$TOOLS_DIR/parallel.py" \
    -i "$INPUT_DIR" -o "$OUTPUT_DIR" \
    -r "$REPORT_DIR" \
    --args="--detect_adapter_for_pe \
    --length_required 50 \
    --qualified_quality_phred 20 \
    --unqualified_percent_limit 20"

echo "Trimming and filtering complete."
ls -lh "$OUTPUT_DIR"


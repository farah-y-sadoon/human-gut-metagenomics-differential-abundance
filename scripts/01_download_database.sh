#!/bin/bash
# Script to download Kraken2 standard database

#SBATCH --job-name=download_db
#SBATCH --time=08:00:00
#SBATCH --account=def-cottenie
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=logs/download_db_%j.txt

# Define directory path
KRAKEN2_DB="../kraken2_db"

# Make kraken2_db directory and change directories
mkdir -p "$KRAKEN2_DB"
cd "$KRAKEN2_DB"

# Download database
echo "Downloading the database..."
wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20251015.tar.gz

# Extract the archive
echo "Extracting ..."
tar -xzf k2_standard_20251015.tar.gz
rm k2_standard_20251015.tar.gz

echo "Download and extraction complete."

# Print database directory contents
echo "Verifying database contents:"
ls -lh "$KRAKEN2_DB"

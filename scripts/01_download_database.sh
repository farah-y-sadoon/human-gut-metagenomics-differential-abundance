#!/bin/bash
# Script to download KrakenUniq standard database

#SBATCH --job-name=download_db
#SBATCH --time=08:00:00
#SBATCH --account=def-cottenie
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=logs/download_db_%j.txt

# Define directory paths
KRAKEN_UNIQ_DB="../krakenuniq_db"

# Make krakenuniq_db directory and change directories
mkdir -p "$KRAKEN_UNIQ_DB"
cd "$KRAKEN_UNIQ_DB"

# Download database
echo "Downloading the database.kdb..."
wget https://genome-idx.s3.amazonaws.com/kraken/uniq/krakendb-2022-06-16-STANDARD/database.kdb

echo "Downloading the archive..."
wget https://genome-idx.s3.amazonaws.com/kraken/uniq/krakendb-2022-06-16-STANDARD/kuniq_standard_minus_kdb.20220616.tgz

# Extract the archive
echo "Extracting archive..."
tar -xzf kuniq_standard_minus_kdb.20220616.tgz

echo "Download and extraction complete."

# Print database directory contents
echo "Verifying database contents:"
ls -lh "$KRAKEN_UNIQ_DB"

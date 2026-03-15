#!/bin/bash
# Script to download reads

#SBATCH --job-name=download_srr
#SBATCH --time=02:00:00
#SBATCH --account=def-cottenie
#SBATCH --cpus-per-task=1
#SBATCH --error=logs/%j_error.txt
#SBATCH --output=logs/%j_output.txt
#SBATCH --mem=4G

# Load module(s)
module load StdEnv/2023 sra-toolkit/3.0.9

# Create output directory
OUT_DIR="./data/raw/"
mkdir -p "$OUT_DIR"

for i in {72..77};
do
echo "Dowloading SRR81469"${i}
prefetch SRR81469${i} -O "$OUT_DIR"
done

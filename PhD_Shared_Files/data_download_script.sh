#!/bin/bash
#SBATCH --job-name=Zhe_Hi-C
#SBATCH --time=40:00:00
#SBATCH --account=def-zoejl
#SBATCH --mem=32G
#SBATCH --cpus-per-task=32
#SBATCH --output=Hi-C.out
#SBATCH --error=Hi-C.err

# Load SRA Toolkit module (if needed on your system)
module load sra-toolkit/3.0.9

# Define array of accession numbers
accessions=(
SRR10991565
SRR10991566
SRR10991567
SRR10991568
SRR10991569
SRR10991570
SRR10991571
#SRR10991572
)

# Create temp directory once
mkdir -p tmp

# Loop through accessions
for acc in "${accessions[@]}"; do
    echo "[INFO] Downloading $acc..."
    
    mkdir -p sra_fastq_files/$acc

    fasterq-dump $acc \
        --split-files \
        --outdir sra_fastq_files/$acc \
        --temp tmp \
        --progress

    echo "[INFO] Finished downloading $acc"
done

echo "[DONE] All downloads complete."




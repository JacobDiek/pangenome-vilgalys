#!/bin/bash
#SBATCH --job-name=dl_pilot
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=04:00:00
#SBATCH --output=/home/jad33864/pangenome/apes/logs/%x_%j.out
#SBATCH --error=/home/jad33864/pangenome/apes/logs/%x_%j.err

# Pilot download: human T2T-CHM13 + both chimpanzee (mPanTro3) haplotypes.
# Versions pinned to the ones used in Yoo et al. 2025 for reproducibility.

set -euo pipefail   # stop on any error, unset variable, or failed pipe step

OUTDIR=/home/jad33864/pangenome/apes/assemblies
mkdir -p "$OUTDIR"
cd "$OUTDIR"

# name -> NCBI GenBank accession
declare -A ACC=(
  [chm13]=GCA_009914755.4     # human T2T-CHM13v2.0
  [ptr_pri]=GCA_028858775.2   # chimpanzee mPanTro3 primary haplotype
  [ptr_alt]=GCA_028858805.2   # chimpanzee mPanTro3 alternate haplotype
)

# ---- Step 1: download ----
module purge
module load NCBI-Datasets-CLI/18.26.0

for name in "${!ACC[@]}"; do
  acc=${ACC[$name]}
  if [[ -s "${name}.fa" ]]; then
    echo "[$(date)] ${name}.fa already exists, skipping"
    continue
  fi
  echo "[$(date)] downloading $name ($acc)"
  datasets download genome accession "$acc" --include genome --filename "${name}.zip"
  unzip -o -q "${name}.zip" -d "${name}_tmp"
  # the zip contains exactly one genomic .fna file; move it out and rename it
  find "${name}_tmp" -name "*.fna" -exec mv {} "${name}.fa" \;
  rm -rf "${name}_tmp" "${name}.zip"
done

# ---- Step 2: index each FASTA ----
module purge
module load SAMtools/1.23.1-GCC-13.3.0

for name in "${!ACC[@]}"; do
  samtools faidx "${name}.fa"
  n_seqs=$(wc -l < "${name}.fa.fai")
  total_bp=$(awk '{s+=$2} END {print s}' "${name}.fa.fai")
  echo "[$(date)] $name: $n_seqs sequences, $total_bp bp"
done

echo "[$(date)] done"
ls -lh "$OUTDIR"

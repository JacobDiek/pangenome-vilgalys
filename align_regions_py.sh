#!/bin/bash
#SBATCH --job-name=align_regions_py
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/home/jad33864/pangenome/apes/logs/%x_%j.out
#SBATCH --error=/home/jad33864/pangenome/apes/logs/%x_%j.err

# Step 3, Python version: bash runs the tools, parse_paf.py does the parsing.

set -euo pipefail

APES=/home/jad33864/pangenome/apes
ALN=$APES/regions/alignments
THREADS=${SLURM_CPUS_PER_TASK:-8}
mkdir -p "$ALN"

module purge
module load minimap2/2.29-GCCcore-13.3.0 SAMtools/1.23.1-GCC-13.3.0 Python/3.12.3-GCCcore-13.3.0

# ---- 1. Combine the four human regions into one query file, named by region
> "$ALN/query.fa"
for name in GYP HOXA LCT MHC; do
  awk -v n="$name" '/^>/ {print ">" n; next} {print}' "$APES/regions/${name}_chm13.fa" >> "$ALN/query.fa"
done

# ---- 2. Align to each chimp haplotype
for hap in ptr_pri ptr_alt; do
  echo "[$(date)] aligning to $hap"
  minimap2 -x asm20 -c --secondary=no -t "$THREADS" \
    "$APES/assemblies/${hap}.fa" "$ALN/query.fa" > "$ALN/query_vs_${hap}.paf"
done

# ---- 3-4. Parse alignments and extract chimp regions (Python)
python3 "$APES/scripts/parse_paf.py" \
  --paf-dir "$ALN" --asm-dir "$APES/assemblies" --out-dir "$APES/regions" \
  --haps ptr_pri ptr_alt --flank 100000 --min-mapq 20

echo "[$(date)] done"

#!/bin/bash
#SBATCH --job-name=index_linear_refs
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=/home/jad33864/pangenome/apes/logs/%x_%j.out
#SBATCH --error=/home/jad33864/pangenome/apes/logs/%x_%j.err

# Make the two linear references bwa maps to, and index them:
#   CHM13 region           -> refs/REGION/chm13.fa    (header >REGION, matches the graph VCF)
#   in-graph chimp region  -> refs/REGION/ptr_*.fa    (header >PTR1_REGION or >PTR2_REGION)
# Each reference is built once; delete its files to force a rebuild.
#
# Usage:  sbatch 2_graphs/index_linear_refs.sh REGION pri|alt   (seconds; fine in interact too)

set -euo pipefail
source /home/jad33864/pangenome/apes/scripts/config.sh
set_region "$@"
module purge
module load $MAP_MODULES

need "$CHM13_FA" "$IN_FA"
mkdir -p "$REFS"

# make_ref SOURCE_FA DEST_FA HEADER: copy with a clean header, then bwa index (only if missing)
make_ref() {
    local src=$1 dest=$2 header=$3
    if [[ ! -s "$dest.bwt" ]]; then
        echo "[$(date)] Writing and indexing $dest"
        awk -v h=">$header" 'NR==1 {$0=h} 1' "$src" > "$dest"
        bwa index "$dest"
    else
        echo "[$(date)] Already indexed: $dest (skipping)"
    fi
}

make_ref "$CHM13_FA" "$CHM13_REF" "$REGION"
make_ref "$IN_FA"    "$IN_REF"    "${IN_TAG}_${REGION}"

echo "[$(date)] Done: $REFS"

#!/bin/bash
#SBATCH --job-name=index_ape_pg
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=24:00:00
#SBATCH --output=/home/jad33864/pangenome/apes/logs/%x_%j.out
#SBATCH --error=/home/jad33864/pangenome/apes/logs/%x_%j.err

# Prepare the downloaded great ape pangenome (0_download/download_ape_pangenome.sh)
# for Giraffe:
#   1. list the graph's reference paths (names needed for surject and region coordinates)
#   2. build the zipcode minimizer index (the downloaded .min predates zipcodes)

set -euo pipefail
source /home/jad33864/pangenome/apes/scripts/config.sh
module purge
need "$GAPG_GBZ" "$GAPG_DIST"

CORES=${SLURM_CPUS_PER_TASK:-16}
mkdir -p "$GAPG_ZIDX"

# ---- 1. Reference path names and lengths
echo "[$(date)] Listing reference paths"
vg paths -M -x "$GAPG_GBZ" > "$GAPG_DIR/paths_metadata.tsv" || echo "WARNING: vg paths -M failed" >&2
vg paths -E -x "$GAPG_GBZ" > "$GAPG_DIR/paths_lengths.tsv"  || echo "WARNING: vg paths -E failed" >&2
echo "Samples with reference paths:"
awk 'NR > 1 {print $3}' "$GAPG_DIR/paths_metadata.tsv" | sort | uniq -c || true
echo "First paths:"
head -5 "$GAPG_DIR/paths_lengths.tsv" || true

# ---- 2. Zipcode minimizer index. If vg rejects the old distance index here, rebuild it
#         first (takes a long time and a lot of memory):
#   vg index -j "$GAPG_DIST" -t "$CORES" "$GAPG_GBZ"
echo "[$(date)] Building zipcode minimizer index"
vg minimizer -d "$GAPG_DIST" -z "$GAPG_ZIP" -o "$GAPG_ZMIN" -t "$CORES" "$GAPG_GBZ"

echo "[$(date)] done"
ls -lh "$GAPG_ZIDX"

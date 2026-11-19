#!/bin/bash
#SBATCH --job-name=index_ape_pg
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=48:00:00
#SBATCH --output=/home/jad33864/pangenome/apes/logs/%x_%j.out
#SBATCH --error=/home/jad33864/pangenome/apes/logs/%x_%j.err

# Prepare the downloaded great ape pangenome (0_download/download_ape_pangenome.sh)
# for Giraffe:
#   1. list the graph's paths (names needed for surject and region coordinates)
#   2. rebuild the distance index with our vg: the downloaded one (vg from early 2024)
#      makes vg 1.66 segfault while building zipcodes (job 50264570)
#   3. build the zipcode minimizer index from the rebuilt distance index
# Each step is skipped if its output already exists, so a rerun after a crash
# picks up where it stopped. Delete a step's output to force it to run again.
#
# Usage:  sbatch 2_graphs/index_ape_pangenome.sh
# Memory and time are guesses for a whole-genome graph: check afterwards with
#   seff JOBID
# If it runs out of memory, resubmit on the high-memory partition:
#   sbatch --partition=highmem_p --mem=400G 2_graphs/index_ape_pangenome.sh

set -euo pipefail
source /home/jad33864/pangenome/apes/scripts/config.sh
module purge
need "$GAPG_GBZ"

CORES=${SLURM_CPUS_PER_TASK:-16}
mkdir -p "$GAPG_ZIDX"

# ---- 1. Path names and lengths (informational: a failure here doesn't stop the job)
if [[ ! -s "$GAPG_DIR/paths_lengths.tsv" ]]; then
    echo "[$(date)] Listing paths"
    vg paths -M -x "$GAPG_GBZ" > "$GAPG_DIR/paths_metadata.tsv" || echo "WARNING: vg paths -M failed" >&2
    vg paths -E -x "$GAPG_GBZ" > "$GAPG_DIR/paths_lengths.tsv"  || echo "WARNING: vg paths -E failed" >&2
else
    echo "[$(date)] Path lists exist, skipping"
fi

# ---- 2. Distance index, rebuilt with the container's vg
if [[ ! -s "$GAPG_DIST" ]]; then
    echo "[$(date)] Building distance index (slow)"
    vg index -j "$GAPG_DIST.tmp" -t "$CORES" "$GAPG_GBZ"
    mv "$GAPG_DIST.tmp" "$GAPG_DIST"     # only a finished index gets the real name
else
    echo "[$(date)] Distance index exists, skipping"
fi
ls -lh "$GAPG_DIST"

# ---- 3. Zipcode minimizer index
echo "[$(date)] Building zipcode minimizer index"
vg minimizer -d "$GAPG_DIST" -z "$GAPG_ZIP" -o "$GAPG_ZMIN" -t "$CORES" "$GAPG_GBZ"

echo "[$(date)] done"
ls -lh "$GAPG_ZIDX"

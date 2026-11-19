#!/bin/bash
#SBATCH --job-name=smoke_ape_pg
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --time=06:00:00
#SBATCH --output=/home/jad33864/pangenome/apes/logs/%x_%j.out
#SBATCH --error=/home/jad33864/pangenome/apes/logs/%x_%j.err

# Smoke test for the great ape pangenome: map the existing HOXA reads (simulated
# from PTR.2 = mPanTro3.2, which IS in this graph, so near-perfect results are
# expected) with the Cactus v2.7.1 vg. Checks that:
#   1. vg_old's giraffe accepts the downloaded .gbz/.min/.dist
#   2. how much memory and time whole-genome mapping takes (check with seff JOBID)
#   3. surject onto hs1 (CHM13) works and which chromosome the reads land on
#
# Usage:  sbatch 3_mapping/smoke_test_ape_pg.sh [COVERAGE]   (default 1: a quick test)

set -euo pipefail
source /home/jad33864/pangenome/apes/scripts/config.sh
set_region HOXA pri          # only to find the existing reads: reads/HOXA_ptr_alt/
set_run "${1:-1}"
module purge

CORES=${SLURM_CPUS_PER_TASK:-8}
OUT=$APES/results/great_ape_smoke/$RUN
HS1_PATHS=$GAPG_DIR/hs1_paths.txt
need "$GAPG_GBZ" "$GAPG_MIN" "$GAPG_DIST" "$SIF_OLD" "$SIM_FQ" "$GAPG_DIR/paths_lengths.tsv"
mkdir -p "$OUT"

echo "vg: $(vg_old version | head -1)"
echo "reads: $SIM_FQ ($(( $(wc -l < "$SIM_FQ") / 4 )) reads)"

# ---- 1. Giraffe to the whole graph
echo "[$(date)] giraffe"
vg_old giraffe -Z "$GAPG_GBZ" -m "$GAPG_MIN" -d "$GAPG_DIST" \
    -f "$SIM_FQ" -t "$CORES" > "$OUT/giraffe.gam"
vg_old stats -a "$OUT/giraffe.gam" > "$OUT/giraffe_stats.txt"
head -8 "$OUT/giraffe_stats.txt"

# ---- 2. Surject onto hs1 only (the graph has two reference genomes, hs1 and hg38)
awk '$1 ~ /^hs1#/ {print $1}' "$GAPG_DIR/paths_lengths.tsv" > "$HS1_PATHS"
echo "[$(date)] surject onto $(wc -l < "$HS1_PATHS") hs1 paths"
vg_old surject -x "$GAPG_GBZ" -F "$HS1_PATHS" -b -t "$CORES" "$OUT/giraffe.gam" \
    > "$OUT/giraffe_hs1.unsorted.bam"

module load SAMtools/1.23.1-GCC-13.3.0
samtools sort -@ "$CORES" -o "$OUT/giraffe_hs1.bam" "$OUT/giraffe_hs1.unsorted.bam"
rm "$OUT/giraffe_hs1.unsorted.bam"
samtools index "$OUT/giraffe_hs1.bam"
samtools flagstat "$OUT/giraffe_hs1.bam" > "$OUT/giraffe_hs1.flagstat"
cat "$OUT/giraffe_hs1.flagstat"

# Which hs1 chromosome did the mapped reads land on? (HOXA is on chr7)
echo "Reads per chromosome:"
samtools idxstats "$OUT/giraffe_hs1.bam" | awk '$3 > 0 {print $1, $3}' | sort -k2,2nr | head

echo "[$(date)] done"
ls -lh "$OUT"

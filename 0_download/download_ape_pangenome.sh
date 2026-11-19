#!/bin/bash
#SBATCH --job-name=dl_ape_pg
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --output=/home/jad33864/pangenome/apes/logs/%x_%j.out
#SBATCH --error=/home/jad33864/pangenome/apes/logs/%x_%j.err

# Download the UCSC 10-way T2T great ape pangenome (Minigraph-Cactus, Cactus v2.7.1):
#   hs1 (CHM13), hg38, HG002 (2 haplotypes), chimp mPanTro3 (2), bonobo mPanPan1 (2),
#   gorilla mGorGor1 (2). Orangutans are in a separate graph.
# Downloads the files Giraffe needs (~49 GB) and the Cactus v2.7.1 container whose vg
# built them (vg_old in config.sh). Skipped on purpose:
#   full.*    og 54 GB, xg 31 GB: Giraffe uses the default (clip) graph
#   MAF/HAL   whole-genome alignments (HAL 13 GB): not needed for mapping
#
# Usage:  sbatch 0_download/download_ape_pangenome.sh
# Safe to rerun: wget -c resumes partial files and skips finished ones.

set -euo pipefail
source /home/jad33864/pangenome/apes/scripts/config.sh

command -v wget >/dev/null || { echo "ERROR: wget not found on this node" >&2; exit 1; }
mkdir -p "$GAPG_DIR"
cd "$GAPG_DIR"

FILES=(
  # Small reference files: which genomes are in the graph and what their paths are called
  "$GAPG_NAME.README.md"
  "$GAPG_NAME.seqfile"        # genome name -> assembly used in the build
  "$GAPG_NAME.stats.tgz"      # graph statistics from Cactus
  chromfile.txt
  rename-hal-from-gca.tsv     # GenBank accession <-> genome name tables
  rename-hal-to-gca.tsv
  # Giraffe inputs (default clip graph)
  "$GAPG_NAME.gbz"            # 7.4 GB: graph + haplotype paths
  "$GAPG_NAME.dist"           # 3.7 GB: distance index
  "$GAPG_NAME.min"            # 37 GB: minimizer index
  # Optional: uncomment when needed
  # "$GAPG_NAME.hapl"         # 2.9 GB: haplotype-sampling index (giraffe --haplotype-name)
  # "$GAPG_NAME.vcf.gz"       # 1.9 GB: graph variants against hs1
  # "$GAPG_NAME.vcf.gz.tbi"
)

# Compare a local file's size with the size the server reports (no md5s exist for these)
check_size() {
    local f=$1 remote local_size
    remote=$(wget --spider -S "$GAPG_URL/$f" 2>&1 \
             | awk 'tolower($1) == "content-length:" {print $2}' | tail -1 | tr -d '\r')
    local_size=$(stat -c %s "$f")
    if [[ -n "$remote" && "$remote" != "$local_size" ]]; then
        echo "ERROR: $f is $local_size bytes, server says $remote (rerun to resume)" >&2
        exit 1
    fi
}

for f in "${FILES[@]}"; do
    echo "[$(date)] $f"
    wget -c -nv "$GAPG_URL/$f"
    check_size "$f"
done

# Container with the vg that built these indexes (apptainer runs on compute nodes only)
if [[ ! -s "$SIF_OLD" ]]; then
    echo "[$(date)] pulling $SIF_OLD_URL"
    apptainer pull "$SIF_OLD" "$SIF_OLD_URL"
fi
echo "vg in $SIF_OLD:"
vg_old version | head -1

# Unpack the stats and print what the graph contains, for the log
mkdir -p stats
tar -xzf "$GAPG_NAME.stats.tgz" -C stats
echo "---- seqfile (genomes in the graph) ----"
cat "$GAPG_NAME.seqfile"
echo "---- chromfile.txt ----"
cat chromfile.txt

echo "[$(date)] done"
ls -lh "$GAPG_DIR"
du -sh "$GAPG_DIR"

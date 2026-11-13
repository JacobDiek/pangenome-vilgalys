#!/bin/bash
#SBATCH --job-name=build_graph
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --output=/home/jad33864/pangenome/apes/logs/%x_%j.out
#SBATCH --error=/home/jad33864/pangenome/apes/logs/%x_%j.err

# Build a Minigraph-Cactus pangenome graph for one region and direction, plus
# the zipcode minimizer index Giraffe needs.
#
# Usage:  sbatch 2_graphs/build_graph.sh REGION pri|alt
#   e.g.  sbatch -J build_graph_HOXA_pri 2_graphs/build_graph.sh HOXA pri
#   pri: CHM13 + PTR.1 (ptr_pri)    alt: CHM13 + PTR.2 (ptr_alt)
# Big regions (MHC) may need more: sbatch --mem=64G --time=08:00:00 ...

set -euo pipefail
source /home/jad33864/pangenome/apes/scripts/config.sh
set_region "$@"
module purge

CORES=${SLURM_CPUS_PER_TASK:-8}
MEM_MB=${SLURM_MEM_PER_NODE:-32768}
MAXMEM="$(( MEM_MB * 9 / 10 ))M"    # 90% of the job's memory, leaving room for Toil's own overhead

JOBSTORE=$GRAPH_DIR/js
SEQFILE=$GRAPH_DIR/seqfile.txt

need "$SIF" "$CHM13_FA" "$IN_FA"
echo "[$(date)] Building $GRAPH_NAME: CHM13 + $IN_SAMPLE ($IN_HAP)"

# ---- 1. Clean copies of the region FASTAs, with a simple ">REGION" header
mkdir -p "$GRAPH_DIR"
awk -v h=">$REGION" 'NR==1 {$0=h} 1' "$CHM13_FA" > "$GRAPH_DIR/${REGION}_chm13.fa"
awk -v h=">$REGION" 'NR==1 {$0=h} 1' "$IN_FA"    > "$GRAPH_DIR/${REGION}_${IN_HAP}.fa"

# ---- 2. Seqfile: sample name, then FASTA. "CHM13" (no dot) = reference haplotype 0;
#         "PTR.1"/"PTR.2" = sample PTR, haplotype 1 or 2.
echo "CHM13" "$GRAPH_DIR/${REGION}_chm13.fa" > "$SEQFILE"
echo "$IN_SAMPLE" "$GRAPH_DIR/${REGION}_${IN_HAP}.fa" >> "$SEQFILE"

# ---- 3. Remove a previous run's jobstore, output and Giraffe index
#         (:? stops the script if a variable is empty, so rm can't hit the wrong folder)
rm -rf "${JOBSTORE:?}" "${GRAPH_OUT:?}" "${ZIDX:?}"

# ---- 4. Minigraph-Cactus (official container)
apptainer exec "$SIF" cactus-pangenome "$JOBSTORE" "$SEQFILE" \
  --outDir "$GRAPH_OUT" --outName "$GRAPH_NAME" \
  --reference CHM13 --gfa --gbz --vcf --giraffe clip \
  --maxCores "$CORES" --maxMemory "$MAXMEM"

# ---- 5. Zipcode minimizer index for Giraffe. Cactus's .min has no zipcodes, so
#         Giraffe would otherwise rebuild one inside out/ on every mapping run.
echo "[$(date)] Building zipcode minimizer index"
mkdir -p "$ZIDX"
vg minimizer -d "$DIST" -z "$ZIP" -o "$ZMIN" -t "$CORES" "$GBZ"

echo "[$(date)] Done: $GRAPH_OUT and $ZIDX"

#!/bin/bash
#SBATCH --job-name=build_hoxa_graph
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --output=/home/jad33864/pangenome/apes/logs/%x_%j.out
#SBATCH --error=/home/jad33864/pangenome/apes/logs/%x_%j.err
set -euo pipefail

module purge

CONTAINER=/home/jad33864/pangenome/containers
SIF=$CONTAINER/cactus_v2.9.9.sif
APES=/home/jad33864/pangenome/apes
GRAPHS=$APES/graphs
REGIONS=$APES/regions
HOXA=$GRAPHS/HOXA_pri
JOBSTOREPATH=$HOXA/js
SEQFILE=$HOXA/seqfile.txt
OUTPUT=$HOXA/out

CORES=${SLURM_CPUS_PER_TASK:-8}
MEM_MB=${SLURM_MEM_PER_NODE:-32768}
MAXMEM="$(( MEM_MB * 9 / 10 ))M"    # 90% of the job's memory, leaving room for Toil's own overhead

# Remove a previous run's jobstore and output
rm -rf "${JOBSTOREPATH:?}" "${OUTPUT:?}"

# ---- Step 1: clean copies of the region FASTAs, with simple ">HOXA" headers
mkdir -p "$HOXA"
for genome in HOXA_chm13 HOXA_ptr_pri; do
    awk 'NR==1 {$0=">HOXA"} 1' "$REGIONS/${genome}.fa" > "$HOXA/${genome}.fa"
done

# ---- Step 2: seqfile listing each sample name and its FASTA
echo "CHM13" "$HOXA/HOXA_chm13.fa" > "$SEQFILE"
echo "PTR.1" "$HOXA/HOXA_ptr_pri.fa" >> "$SEQFILE"

# ---- Step 3: build the graph with Minigraph-Cactus (official container)
# Check to see if cactus_v2.9.9.sif container exists
[[ -s "$SIF" ]] || { echo "ERROR: container not found: $SIF" >&2; exit 1; }

apptainer exec "$SIF" cactus-pangenome "$JOBSTOREPATH" "$SEQFILE" --outDir "$OUTPUT" --outName HOXA_pri \
  --reference CHM13 --gfa --gbz --vcf --giraffe clip \
  --maxCores "$CORES" --maxMemory "$MAXMEM"
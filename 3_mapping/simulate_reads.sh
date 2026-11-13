#!/bin/bash
#SBATCH --job-name=simulate_reads
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/home/jad33864/pangenome/apes/logs/%x_%j.out
#SBATCH --error=/home/jad33864/pangenome/apes/logs/%x_%j.err

# Simulate short reads from the HELD-OUT chimp haplotype (the one not in the graph).
# Read settings (length, seed, error rates, paired) live in config.sh.
#
# Usage:  sbatch 3_mapping/simulate_reads.sh REGION pri|alt [COVERAGE]
#   e.g.  sbatch 3_mapping/simulate_reads.sh HOXA pri 30     (reads from PTR.2)
# Existing reads for the same run tag are kept; delete them to regenerate.

set -euo pipefail
source /home/jad33864/pangenome/apes/scripts/config.sh
set_region "${1:-}" "${2:-}"
set_run "${3:-}"
module purge

need "$SIF" "$HELD_FA"
mkdir -p "$READS"

if [[ -s "$SIM_GAM" && -s "$SIM_FQ" ]]; then
    echo "[$(date)] Reads already exist for $RUN in $READS (skipping; delete them to regenerate)"
    exit 0
fi

# reads = region length x coverage / bases per read (a pair carries 2 x READ_LEN)
HELD_LEN=$(grep -v '>' "$HELD_FA" | tr -d '\n' | wc -c)
if [[ "$PAIRED" == true ]]; then
    N_READS=$(( HELD_LEN * COVERAGE / (2 * READ_LEN) ))
    SIM_PE="-p $FRAG_LEN -v $FRAG_SD"
else
    N_READS=$(( HELD_LEN * COVERAGE / READ_LEN ))
    SIM_PE=""
fi
echo "[$(date)] $REGION $HELD_HAP ($HELD_TAG): $HELD_LEN bp; simulating $N_READS reads (run $RUN)"

# vg sim needs a graph, so make a plain linear graph of the held-out sequence.
# The GAM keeps each read's true origin on it, for scoring placement later.
SRC=$READS/${HELD_TAG}_${REGION}
awk -v h=">${HELD_TAG}_${REGION}" 'NR==1 {$0=h} 1' "$HELD_FA" > "$SRC.fa"
vg construct -r "$SRC.fa" > "$SRC.vg"
vg convert -x "$SRC.vg" > "$SRC.xg"     # xg: gives vg sim path positions

vg sim -x "$SRC.xg" -n "$N_READS" -l "$READ_LEN" -s "$SEED" \
    -e "$SUB_RATE" -i "$INDEL_RATE" $SIM_PE -a > "$SIM_GAM"

# The same reads as FASTQ, so every mapper gets identical input
vg view -X "$SIM_GAM" > "$SIM_FQ"

if [[ "$PAIRED" == true ]]; then
    # bwa -p pairs consecutive reads only if mates share a name (or end /1, /2)
    echo "First two read names (mates should match):"
    awk 'NR % 4 == 1' "$SIM_FQ" | head -n 2
fi

echo "[$(date)] Done: $SIM_GAM and $SIM_FQ"

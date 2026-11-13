#!/bin/bash
#SBATCH --job-name=map_reads
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/home/jad33864/pangenome/apes/logs/%x_%j.out
#SBATCH --error=/home/jad33864/pangenome/apes/logs/%x_%j.err

# Map the simulated held-out reads three ways and summarize:
#   1. vg giraffe -> pangenome graph (CHM13 + in-graph chimp), surjected to CHM13
#   2. bwa mem    -> CHM13 region                 ("just use human")
#   3. bwa mem    -> in-graph chimp region         (same-species best case)
#
# Needs, from earlier steps:
#   2_graphs/build_graph.sh REGION DIR          (graph + zipcode index)
#   2_graphs/index_linear_refs.sh REGION DIR    (bwa references)
#   3_mapping/simulate_reads.sh REGION DIR COV  (reads)
#
# Usage:  sbatch 3_mapping/map_reads.sh REGION pri|alt [COVERAGE]
#   e.g.  sbatch 3_mapping/map_reads.sh HOXA pri 30

set -euo pipefail
source /home/jad33864/pangenome/apes/scripts/config.sh
set_region "${1:-}" "${2:-}"
set_run "${3:-}"
module purge
module load $MAP_MODULES

THREADS=${SLURM_CPUS_PER_TASK:-4}

need "$SIF" "$GBZ" "$DIST" "$ZMIN" "$ZIP" "$SIM_FQ" "$CHM13_REF.bwt" "$IN_REF.bwt"
if [[ "$GBZ" -nt "$ZMIN" ]]; then
    # -nt = "newer than": the graph was rebuilt after its Giraffe index was made
    echo "ERROR: $ZMIN is older than $GBZ; rerun 2_graphs/build_graph.sh $REGION $DIR" >&2
    exit 1
fi

mkdir -p "$RUN_OUT"
echo "[$(date)] Mapping $RUN reads from $HELD_TAG to $GRAPH_NAME, CHM13 and $IN_HAP"

# ---- 1. Graph mapping: giraffe -> GAM -> BAM on CHM13
echo "[$(date)] Mapping with vg giraffe"
VG_PE=""; BWA_PE=""
if [[ "$PAIRED" == true ]]; then VG_PE="-i"; BWA_PE="-p"; fi   # interleaved pairs
vg giraffe -Z "$GBZ" -m "$ZMIN" -z "$ZIP" -d "$DIST" -f "$SIM_FQ" $VG_PE \
    -t "$THREADS" > "$RUN_OUT/giraffe.gam"

# surject projects each graph alignment onto the reference path (CHM13#0#REGION)
vg surject -x "$GBZ" -b $VG_PE -t "$THREADS" "$RUN_OUT/giraffe.gam" \
    | samtools sort -@ "$THREADS" -o "$RUN_OUT/giraffe.sorted.bam"
samtools index "$RUN_OUT/giraffe.sorted.bam"

# ---- 2. Linear mapping: bwa mem to each reference
for ref in chm13 "$IN_HAP"; do
    echo "[$(date)] Mapping with bwa mem to $ref"
    bwa mem -t "$THREADS" $BWA_PE "$REFS/$ref.fa" "$SIM_FQ" \
        | samtools sort -@ "$THREADS" -o "$RUN_OUT/bwa_${ref}.sorted.bam"
    samtools index "$RUN_OUT/bwa_${ref}.sorted.bam"
done

# ---- 3. Summarize: flagstat + MAPQ breakdown
# -F 0x900 drops secondary (0x100) and supplementary (0x800) lines, so each read
# counts once. FLAG bit 4 = unmapped: int($2 / 4) % 2 tests it without gawk.
echo "[$(date)] Summarizing"
SUMMARY=$RUN_OUT/summary.tsv
printf "method\treads\tmapped\tpct_mapped\tmapq0\tmapq1_29\tmapq30_59\tmapq60\n" > "$SUMMARY"
for m in giraffe bwa_chm13 "bwa_${IN_HAP}"; do
    bam="$RUN_OUT/$m.sorted.bam"
    samtools flagstat "$bam" > "$RUN_OUT/$m.flagstat.txt"
    samtools view -F 0x900 "$bam" | awk -v m="$m" '
        { n++
          if (int($2 / 4) % 2 == 1) next
          mapped++
          if ($5 == 0) q0++; else if ($5 < 30) q1++; else if ($5 < 60) q30++; else q60++ }
        END { printf "%s\t%d\t%d\t%.2f\t%d\t%d\t%d\t%d\n",
              m, n, mapped, (n ? 100 * mapped / n : 0), q0, q1, q30, q60 }' >> "$SUMMARY"
done

# Giraffe before surject: reads aligned to the graph itself (GAM). Reads aligned
# only to chimp-specific sequence have no CHM13 position, so surject marks them
# unmapped; the difference between this row and the giraffe row is that loss.
vg stats -a "$RUN_OUT/giraffe.gam" > "$RUN_OUT/giraffe.gam_stats.txt"
awk -F': ' '/^Total alignments/ {n = $2} /^Total aligned/ {a = $2}
    END { printf "giraffe_graph\t%d\t%d\t%.2f\tNA\tNA\tNA\tNA\n", n, a, (n ? 100 * a / n : 0) }' \
    "$RUN_OUT/giraffe.gam_stats.txt" >> "$SUMMARY"

column -t "$SUMMARY"
echo "(giraffe_graph = aligned to the graph before surject; giraffe = after projecting onto CHM13)"
echo "[$(date)] Done. Results in $RUN_OUT"

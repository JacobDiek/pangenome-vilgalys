#!/bin/bash
#SBATCH --job-name=map_hoxa_reads
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/home/jad33864/pangenome/apes/logs/%x_%j.out
#SBATCH --error=/home/jad33864/pangenome/apes/logs/%x_%j.err

# Simulate short reads from the held-out chimp haplotype (PTR.2 = ptr_alt HOXA),
# then map the SAME reads three ways and summarize:
#   1. vg giraffe -> HOXA_pri graph (CHM13 + PTR.1), surjected to a BAM on CHM13
#   2. bwa mem    -> CHM13 HOXA region            (the "just use human" option)
#   3. bwa mem    -> ptr_pri HOXA region (PTR.1)  (same-species best case)
# PTR.2 is in none of the three references, so the reads act like a new chimp.
#
# Re-running: simulation and bwa indexing are skipped if their outputs already
# exist, so changing only the mapping steps reuses the same reads.

set -euo pipefail

module purge
module load BWA/0.7.18-GCCcore-13.3.0 SAMtools/1.23.1-GCC-13.3.0

# ---------- Settings ----------
COVERAGE=30           # target read depth over the PTR.2 region
READ_LEN=150
SEED=42               # fixed seed so the same reads come out every run
SUB_RATE=0.01         # substitution error rate
INDEL_RATE=0.002      # indel error rate
PAIRED=false          # true = paired-end (check mate names first, see step 2)
FRAG_LEN=300          # paired-end only: mean fragment length
FRAG_SD=50            # paired-end only: fragment length std dev

THREADS=${SLURM_CPUS_PER_TASK:-4}

# ---------- Paths ----------
APES=/home/jad33864/pangenome/apes
SIF=/home/jad33864/pangenome/containers/cactus_v2.9.9.sif
GRAPH=$APES/graphs/HOXA_pri/out
GBZ=$GRAPH/HOXA_pri.gbz
DIST=$GRAPH/HOXA_pri.dist

# Giraffe wants a minimizer index with zipcodes; Cactus's HOXA_pri.min has none.
# Build it here, OUTSIDE out/ (build_hoxa_graph.sh wipes out/ on every rebuild).
ZIDX=$APES/graphs/HOXA_pri/giraffe_idx
ZMIN=$ZIDX/HOXA_pri.withzip.min
ZIP=$ZIDX/HOXA_pri.zipcodes

ALT_SRC=$APES/regions/HOXA_ptr_alt.fa                # held-out haplotype (read source)
CHM13_SRC=$APES/graphs/HOXA_pri/HOXA_chm13.fa        # header is already >HOXA
PRI_SRC=$APES/graphs/HOXA_pri/HOXA_ptr_pri.fa        # header >HOXA, renamed below

REFS=$APES/refs/HOXA                                 # linear references + bwa indexes
READS=$APES/reads/HOXA_ptr_alt                       # simulated reads (shared by all mappers)

# vg comes from the Cactus container: the same vg that built the graph indexes.
# This function makes "vg ..." below run inside the container.
# --cleanenv stops the loaded modules' settings from leaking into the container.
vg() { apptainer exec --cleanenv "$SIF" vg "$@"; }

[[ -s "$SIF" ]] || { echo "ERROR: container not found: $SIF" >&2; exit 1; }
for f in "$GBZ" "$DIST" "$ALT_SRC" "$CHM13_SRC" "$PRI_SRC"; do
    [[ -s "$f" ]] || { echo "ERROR: missing input: $f" >&2; exit 1; }
done

# ---------- Read count from coverage ----------
# reads = region length x coverage / bases per read (a pair carries 2 x READ_LEN)
ALT_LEN=$(grep -v '>' "$ALT_SRC" | tr -d '\n' | wc -c)
if [ "$PAIRED" = true ]; then
    N_READS=$(( ALT_LEN * COVERAGE / (2 * READ_LEN) ))   # number of pairs
    SIM_PE="-p $FRAG_LEN -v $FRAG_SD"
    VG_PE="-i"          # giraffe/surject: GAM is interleaved pairs
    BWA_PE="-p"         # bwa mem: FASTQ is interleaved pairs
    RUN=cov${COVERAGE}_s${SEED}_pe
else
    N_READS=$(( ALT_LEN * COVERAGE / READ_LEN ))
    SIM_PE=""
    VG_PE=""
    BWA_PE=""
    RUN=cov${COVERAGE}_s${SEED}
fi

OUT=$APES/results/HOXA_pri/$RUN
mkdir -p "$REFS" "$READS" "$OUT"
echo "[$(date)] PTR.2 region: $ALT_LEN bp; simulating $N_READS reads (run tag $RUN)"

# ---------- 1. Simulate reads from PTR.2 ----------
# vg sim needs a graph, so build a plain linear graph from the PTR.2 sequence.
# The GAM keeps each read's true origin on PTR.2, for scoring placement later.
if [ ! -s "$READS/$RUN.gam" ]; then
    echo "[$(date)] Simulating reads from PTR.2"
    awk 'NR==1 {$0=">PTR2_HOXA"} 1' "$ALT_SRC" > "$READS/PTR2_HOXA.fa"
    vg construct -r "$READS/PTR2_HOXA.fa" > "$READS/PTR2_HOXA.vg"
    vg convert -x "$READS/PTR2_HOXA.vg" > "$READS/PTR2_HOXA.xg"   # xg: path positions for vg sim
    vg sim -x "$READS/PTR2_HOXA.xg" -n "$N_READS" -l "$READ_LEN" -s "$SEED" \
        -e "$SUB_RATE" -i "$INDEL_RATE" $SIM_PE -a \
        > "$READS/$RUN.gam"
else
    echo "[$(date)] Reads already simulated: $READS/$RUN.gam (skipping)"
fi

# ---------- 2. Same reads as FASTQ (for giraffe -f and bwa) ----------
if [ ! -s "$READS/$RUN.fq" ]; then
    vg view -X "$READS/$RUN.gam" > "$READS/$RUN.fq"
fi
if [ "$PAIRED" = true ]; then
    # bwa -p pairs consecutive reads only if mates share a name (or end /1, /2)
    echo "First two read names (mates should match):"
    awk 'NR % 4 == 1' "$READS/$RUN.fq" | head -n 2
fi

# ---------- 3. Linear references + bwa indexes (built once) ----------
cp "$CHM13_SRC" "$REFS/chm13.fa"                                   # keeps >HOXA (matches the VCF)
awk 'NR==1 {$0=">PTR1_HOXA"} 1' "$PRI_SRC" > "$REFS/ptr_pri.fa"    # distinct name in BAMs
for ref in chm13 ptr_pri; do
    if [ ! -s "$REFS/$ref.fa.bwt" ]; then
        echo "[$(date)] bwa index $ref"
        bwa index "$REFS/$ref.fa"
    fi
done

# ---------- 4. Graph mapping: giraffe -> GAM -> BAM on CHM13 ----------
# Zipcode minimizer index: built once, and again whenever the graph is rebuilt
# (-nt = "newer than": the GBZ is newer than the index, so the index is stale).
if [ ! -s "$ZMIN" ] || [ ! -s "$ZIP" ] || [ "$GBZ" -nt "$ZMIN" ]; then
    echo "[$(date)] Building zipcode minimizer index"
    mkdir -p "$ZIDX"
    vg minimizer -d "$DIST" -z "$ZIP" -o "$ZMIN" -t "$THREADS" "$GBZ"
fi

echo "[$(date)] Mapping with vg giraffe"
vg giraffe -Z "$GBZ" -m "$ZMIN" -z "$ZIP" -d "$DIST" -f "$READS/$RUN.fq" $VG_PE \
    -t "$THREADS" > "$OUT/giraffe.gam"

# surject projects each graph alignment onto the reference path (CHM13#0#HOXA)
vg surject -x "$GBZ" -b $VG_PE -t "$THREADS" "$OUT/giraffe.gam" \
    | samtools sort -@ "$THREADS" -o "$OUT/giraffe.sorted.bam"
samtools index "$OUT/giraffe.sorted.bam"

# ---------- 5. Linear mapping: bwa mem to each reference ----------
for ref in chm13 ptr_pri; do
    echo "[$(date)] Mapping with bwa mem to $ref"
    bwa mem -t "$THREADS" $BWA_PE "$REFS/$ref.fa" "$READS/$RUN.fq" \
        | samtools sort -@ "$THREADS" -o "$OUT/bwa_${ref}.sorted.bam"
    samtools index "$OUT/bwa_${ref}.sorted.bam"
done

# ---------- 6. Summarize: flagstat + MAPQ breakdown ----------
# -F 0x900 drops secondary (0x100) and supplementary (0x800) lines, so each read
# counts once. FLAG bit 4 = unmapped: int($2 / 4) % 2 tests it without gawk.
echo "[$(date)] Summarizing"
printf "method\treads\tmapped\tpct_mapped\tmapq0\tmapq1_29\tmapq30_59\tmapq60\n" > "$OUT/summary.tsv"
for m in giraffe bwa_chm13 bwa_ptr_pri; do
    bam="$OUT/$m.sorted.bam"
    samtools flagstat "$bam" > "$OUT/$m.flagstat.txt"
    samtools view -F 0x900 "$bam" | awk -v m="$m" '
        { n++
          if (int($2 / 4) % 2 == 1) next
          mapped++
          if ($5 == 0) q0++; else if ($5 < 30) q1++; else if ($5 < 60) q30++; else q60++ }
        END { printf "%s\t%d\t%d\t%.2f\t%d\t%d\t%d\t%d\n",
              m, n, mapped, (n ? 100 * mapped / n : 0), q0, q1, q30, q60 }' >> "$OUT/summary.tsv"
done

# Giraffe before surject: reads aligned to the graph itself (GAM). Reads aligned
# only to chimp-specific sequence have no CHM13 position, so surject marks them
# unmapped; the difference between this row and the giraffe BAM row is that loss.
vg stats -a "$OUT/giraffe.gam" > "$OUT/giraffe.gam_stats.txt"
awk -F': ' '/^Total alignments/ {n = $2} /^Total aligned/ {a = $2}
    END { printf "giraffe_graph\t%d\t%d\t%.2f\tNA\tNA\tNA\tNA\n", n, a, (n ? 100 * a / n : 0) }' \
    "$OUT/giraffe.gam_stats.txt" >> "$OUT/summary.tsv"

column -t "$OUT/summary.tsv"
echo "(giraffe_graph = aligned to the graph before surject; giraffe = after projecting onto CHM13)"
echo "[$(date)] Done. Results in $OUT"

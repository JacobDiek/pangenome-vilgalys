#!/bin/bash
#SBATCH --job-name=lct_benchmark
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jad33864@uga.edu

# Simulate reads from the LCT enhancer graph, map them with vg giraffe
# (pangenome) and bwa-mem (linear), then compare with samtools flagstat.

set -euo pipefail

ml vg/1.65.0
ml BWA/0.7.18-GCCcore-13.3.0
ml SAMtools/1.23.1-GCC-13.3.0

# ---------- Settings ----------
N_READS=1000          # number of reads (or read pairs, if PAIRED=true)
READ_LEN=150
SEED=42               # fixed seed so the same reads come out every run
SUB_RATE=0.01         # substitution error rate
INDEL_RATE=0.002      # indel error rate
PAIRED=false          # set to true for paired-end
FRAG_LEN=300          # paired-end only: mean fragment length
FRAG_SD=50            # paired-end only: fragment length std dev

THREADS=$SLURM_CPUS_PER_TASK

# ---------- Paths ----------
BASE=~/pangenome
REF=$BASE/ref/lct_enhancer.fa
IDX=$BASE/vg-builds/lct_index
GBZ=$IDX.giraffe.gbz
MIN=$IDX.shortread.withzip.min
ZIP=$IDX.shortread.zipcodes
DIST=$IDX.dist

RUN=sim_n${N_READS}_s${SEED}      # tag so different settings don't overwrite each other

# Paired-end flags (empty when single-end)
if [ "$PAIRED" = true ]; then
    SIM_PE="-p $FRAG_LEN -v $FRAG_SD"
    VG_PE="-i"          # giraffe/surject: input GAM is interleaved pairs
    BWA_PE="-p"         # bwa mem: input FASTQ is interleaved pairs
    RUN=${RUN}_pe
else
    SIM_PE=""
    VG_PE=""
    BWA_PE=""
fi

READS=$BASE/reads
OUT=$BASE/results/$RUN
mkdir -p "$READS" "$OUT"

# ---------- 1. Simulate reads (GAM keeps each read's true origin) ----------
echo "[$(date)] Simulating reads"
vg sim -x "$GBZ" -n "$N_READS" -l "$READ_LEN" -s "$SEED" \
    -e "$SUB_RATE" -i "$INDEL_RATE" $SIM_PE -a \
    > "$READS/$RUN.gam"

# ---------- 2. Convert the same reads to FASTQ for bwa ----------
echo "[$(date)] Converting GAM to FASTQ"
vg view -X "$READS/$RUN.gam" > "$READS/$RUN.fq"

# ---------- 3. Pangenome mapping: giraffe -> GAM -> BAM ----------
echo "[$(date)] Mapping with vg giraffe"
vg giraffe -Z "$GBZ" -m "$MIN" -z "$ZIP" -d "$DIST" \
    -G "$READS/$RUN.gam" $VG_PE -t "$THREADS" \
    > "$OUT/giraffe.gam"

vg surject -x "$GBZ" -b $VG_PE -t "$THREADS" "$OUT/giraffe.gam" \
    | samtools sort -@ "$THREADS" -o "$OUT/giraffe.sorted.bam"
samtools index "$OUT/giraffe.sorted.bam"

# ---------- 4. Linear mapping: bwa index (once) + bwa mem ----------
if [ ! -f "$REF.bwt" ]; then
    echo "[$(date)] Building bwa index"
    bwa index "$REF"
fi

echo "[$(date)] Mapping with bwa mem"
bwa mem -t "$THREADS" $BWA_PE "$REF" "$READS/$RUN.fq" \
    | samtools sort -@ "$THREADS" -o "$OUT/bwa.sorted.bam"
samtools index "$OUT/bwa.sorted.bam"

# ---------- 5. Compare ----------
echo "[$(date)] Running flagstat"
samtools flagstat "$OUT/giraffe.sorted.bam" > "$OUT/giraffe.flagstat.txt"
samtools flagstat "$OUT/bwa.sorted.bam"     > "$OUT/bwa.flagstat.txt"

echo "===== giraffe ====="; cat "$OUT/giraffe.flagstat.txt"
echo "===== bwa mem ====="; cat "$OUT/bwa.flagstat.txt"

echo "[$(date)] Done. Results in $OUT"

#!/bin/bash
#SBATCH --job-name=align_regions
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --output=/home/jad33864/pangenome/apes/logs/%x_%j.out
#SBATCH --error=/home/jad33864/pangenome/apes/logs/%x_%j.err

# Step 3: align the four CHM13 regions to both chimpanzee haplotypes,
# finds where each region sits in chimp, and extract the chimp sequences.

set -euo pipefail

APES=/home/jad33864/pangenome/apes
ASM=$APES/assemblies
REGIONS=$APES/regions
ALN=$APES/regions/alignments
FLANK=100000                          # must match FLANK in define_regions.sh
MINQ=20                               # minimum mapping quality for an anchor alignment
THREADS=${SLURM_CPUS_PER_TASK:-8}     # use the CPUs Slurm gave us (8 if run outside Slurm)
mkdir -p "$ALN"

module purge
module load minimap2/2.29-GCCcore-13.3.0 SAMtools/1.23.1-GCC-13.3.0

# ---- 1. Combine the four human regions into one query file, named by region
> "$ALN/query.fa"
for name in GYP HOXA LCT MHC; do
  awk -v n="$name" '/^>/ {print ">" n; next} {print}' "$REGIONS/${name}_chm13.fa" >> "$ALN/query.fa"
done
echo "[$(date)] query sequences: $(grep -c '>' "$ALN/query.fa")"

# ---- 2. Align the regions to each chimp haplotype (PAF output)
for hap in ptr_pri ptr_alt; do
  echo "[$(date)] aligning to $hap"
  minimap2 -x asm20 -c --secondary=no -t "$THREADS" \
    "$ASM/${hap}.fa" "$ALN/query.fa" > "$ALN/query_vs_${hap}.paf"
  echo "[$(date)] $hap: $(wc -l < "$ALN/query_vs_${hap}.paf") alignment blocks"
done

# ---- 3. Use the unique flanks to pin down each region's chimp coordinates
echo -e "chimp_chrom\tstart\tend\tregion\thap\tstrand\tchimp_len\thuman_len\tanchors" > "$REGIONS/regions_chimp.tsv"
for hap in ptr_pri ptr_alt; do
  awk -F'\t' -v flank="$FLANK" -v minq="$MINQ" -v hap="$hap" '
    $12 >= minq && ($3 < flank || $4 > $2 - flank) {
      key = $1 SUBSEP $6
      bp[key] += $10
      if (!(key in lo) || $8 < lo[key]) lo[key] = $8
      if ($9 > hi[key]) hi[key] = $9
      if ($5 == "+") plus[key] += $10; else minus[key] += $10
      if ($3 < flank)       left[key] = 1
      if ($4 > $2 - flank)  right[key] = 1
      qlen[$1] = $2
    }
    END {
      for (key in bp) {
        split(key, k, SUBSEP)
        if (bp[key] > best[k[1]]) { best[k[1]] = bp[key]; bestkey[k[1]] = key }
      }
      for (q in bestkey) {
        key = bestkey[q]; split(key, k, SUBSEP)
        strand  = (plus[key] >= minus[key]) ? "+" : "-"
        anchors = (left[key] && right[key]) ? "both" : (left[key] ? "left_only" : "right_only")
        print k[2] "\t" lo[key] "\t" hi[key] "\t" q "\t" hap "\t" strand "\t" hi[key] - lo[key] "\t" qlen[q] "\t" anchors
      }
    }
  ' "$ALN/query_vs_${hap}.paf"
done | sort -k4,4 -k5,5 >> "$REGIONS/regions_chimp.tsv"

# ---- 4. Extract each chimp region (reverse-complemented if it aligned to the minus strand)
echo "[$(date)] chimp regions:"
while IFS=$'\t' read -r chr start end name hap strand clen hlen anchors; do
  [[ "$chr" == "chimp_chrom" ]] && continue
  flip=""
  [[ "$strand" == "-" ]] && flip="-i"
  samtools faidx $flip "$ASM/${hap}.fa" "${chr}:$((start + 1))-${end}" > "$REGIONS/${name}_${hap}.fa"
  ratio=$(awk -v c="$clen" -v h="$hlen" 'BEGIN {printf "%.2f", c / h}')
  flag=""
  if [[ "$anchors" != "both" ]] || awk -v r="$ratio" 'BEGIN {exit !(r < 0.5 || r > 2)}'; then
    flag="  <-- CHECK"
  fi
  echo "  $name $hap  ${chr}:$((start + 1))-${end} ($strand)  $((clen / 1000)) kb vs human $((hlen / 1000)) kb  ratio=$ratio  anchors=$anchors$flag"
done < "$REGIONS/regions_chimp.tsv"

echo "[$(date)] done"
ls -lh "$REGIONS"

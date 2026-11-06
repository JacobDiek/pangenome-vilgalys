#!/bin/bash
#SBATCH --job-name=define_regions
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=/home/jad33864/pangenome/apes/logs/%x_%j.out
#SBATCH --error=/home/jad33864/pangenome/apes/logs/%x_%j.err

# Step 2: define the four study regions in T2T-CHM13 coordinates using
# NCBI's gene annotation, then extract each region's sequence.

set -euo pipefail

APES=/home/jad33864/pangenome/apes
ANNOT=$APES/annotation
REGIONS=$APES/regions
FLANK=100000          # extra sequence added on each side of the anchor genes (100 kb)
mkdir -p "$ANNOT" "$REGIONS"

# region name -> anchor genes; each region spans from the leftmost to the
# rightmost of its anchor genes, plus FLANK on both sides
declare -A GENES=(
  [HOXA]="HOXA1 HOXA13"          # whole HOXA cluster (easy control)
  [LCT]="MCM6 LCT"               # LCT gene + MCM6, which holds the lactase enhancer
  [MHC]="MOG COL11A2"            # classical MHC boundaries (class I through class II)
  [GYP]="GYPA GYPB GYPE"         # glycophorin cluster (segmental duplications)
)

# ---- 1. Download CHM13's RefSeq gene annotation (GFF3 format)
module purge
module load NCBI-Datasets-CLI/18.26.0
cd "$ANNOT"
if [[ ! -s chm13.gff ]]; then
  echo "[$(date)] downloading CHM13 annotation"
  datasets download genome accession GCF_009914755.1 --include gff3 --filename chm13_annot.zip
  unzip -o -q chm13_annot.zip -d tmp
  find tmp -name "*.gff" -exec mv {} chm13.gff \;
  rm -rf tmp chm13_annot.zip
fi

# keep only the "gene" rows so the lookups below scan a small file
awk -F'\t' '$3=="gene"' chm13.gff > chm13_genes.gff

# ---- 2. Translate chromosome names: the annotation uses RefSeq names (NC_...),
# FASTA uses GenBank names (CP...). Match them by chromosome length.
awk -F'\t' '$3=="region" && $4==1 {print $1"\t"$5}' chm13.gff | sort -t$'\t' -k2,2 > refseq_lengths.tsv
cut -f1,2 "$APES/assemblies/chm13.fa.fai"                       | sort -t$'\t' -k2,2 > genbank_lengths.tsv
join -t$'\t' -1 2 -2 2 refseq_lengths.tsv genbank_lengths.tsv | cut -f2,3 > nc_to_cp.tsv
echo "[$(date)] matched $(wc -l < nc_to_cp.tsv) chromosome names (expect 25)"

# ---- 3. Look up each anchor gene's position
> gene_hits.tsv
for region in HOXA LCT MHC GYP; do
  for gene in ${GENES[$region]}; do
    hits=$(awk -F'\t' -v g="$gene" '$9 ~ ("(^|;)Name=" g ";") {print $1"\t"$4"\t"$5}' chm13_genes.gff)
    n=$(printf "%s" "$hits" | grep -c . || true)
    if [[ $n -eq 0 ]]; then
      echo "ERROR: gene $gene not found in annotation" >&2
      exit 1
    elif [[ $n -gt 1 ]]; then
      echo "WARNING: gene $gene has $n entries, using the first:" >&2
      echo "$hits" >&2
    fi
    echo -e "$region\t$gene\t$(echo "$hits" | head -1)" >> gene_hits.tsv
  done
done
echo "[$(date)] anchor genes (region, gene, RefSeq chrom, start, end):"
cat gene_hits.tsv

# ---- 4. Combine anchor genes into one interval per region (BED format)
# BED is 0-based and half-open: start = (1-based start - 1)
awk -F'\t' -v flank="$FLANK" '
  NR==FNR { cp[$1]=$2; next }                 # first file: NC -> CP lookup table
  {
    r=$1
    if (!(r in lo) || $4 < lo[r]) lo[r]=$4    # leftmost gene start
    if ($5 > hi[r])               hi[r]=$5    # rightmost gene end
    chrom[r]=cp[$3]
  }
  END {
    for (r in lo) {
      s = lo[r] - 1 - flank; if (s < 0) s = 0
      print chrom[r] "\t" s "\t" hi[r] + flank "\t" r
    }
  }
' nc_to_cp.tsv gene_hits.tsv | sort -t$'\t' -k4,4 > "$REGIONS/regions_chm13.bed"

# ---- 5. Extract each region's sequence from CHM13
module purge
module load SAMtools/1.23.1-GCC-13.3.0
cd "$REGIONS"
echo "[$(date)] regions:"
while read -r chr start end name; do
  samtools faidx "$APES/assemblies/chm13.fa" "${chr}:$((start + 1))-${end}" > "${name}_chm13.fa"
  echo "  $name  ${chr}:$((start + 1))-${end}  $(( (end - start) / 1000 )) kb"
done < regions_chm13.bed

echo "[$(date)] done"
ls -lh "$REGIONS"

#!/bin/bash
#SBATCH --job-name=define_regions
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=/home/jad33864/pangenome/apes/logs/%x_%j.out
#SBATCH --error=/home/jad33864/pangenome/apes/logs/%x_%j.err

# Step 2: define the study regions (listed in regions.tsv) in T2T-CHM13
# coordinates from NCBI's gene annotation, then extract each region's sequence.
# Bash runs the tools (download, samtools); define_regions.py does the parsing.
#
# Usage:  sbatch 1_data_prep/define_regions.sh

set -euo pipefail
source /home/jad33864/pangenome/apes/scripts/config.sh

ANNOT=$APES/annotation
CHM13=$APES/assemblies/chm13.fa
mkdir -p "$ANNOT" "$REGIONS_DIR"
need "$REGIONS_TSV" "$CHM13" "$CHM13.fai"

# ---- 1. Download CHM13's RefSeq gene annotation (GFF3), once
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

# ---- 2-4. Anchor genes -> one BED interval per region (Python)
module purge
module load Python/3.12.3-GCCcore-13.3.0 SAMtools/1.23.1-GCC-13.3.0
python3 "$APES/scripts/1_data_prep/define_regions.py" \
  --gff "$ANNOT/chm13.gff" --fai "$CHM13.fai" \
  --regions "$REGIONS_TSV" --flank "$FLANK" \
  --annot-dir "$ANNOT" --out-bed "$REGIONS_DIR/regions_chm13.bed"

# ---- 5. Extract each region's sequence from CHM13 (BED is 0-based; samtools is 1-based)
cd "$REGIONS_DIR"
echo "[$(date)] regions:"
while read -r chr start end name; do
  samtools faidx "$CHM13" "${chr}:$((start + 1))-${end}" > "${name}_chm13.fa"
  echo "  $name  ${chr}:$((start + 1))-${end}  $(( (end - start) / 1000 )) kb"
done < regions_chm13.bed

echo "[$(date)] done"
ls -lh "$REGIONS_DIR"

#!/bin/bash
# build_pangenome.sh
#
# Builds a vg/Giraffe pangenome graph + indexes from an EXISTING FASTA + VCF
# Requires: Docker Desktop running, quay.io/vgteam/vg image pulled.

set -e   # exit immediately

echo ""
echo "=== Pangenome Graph Builder ==="
echo ""
echo "Select FASTA file"

# Pick the FASTA

FASTA_PATH=$(osascript -e 'POSIX path of (choose file with prompt "Select your reference FASTA file:")' 2>/dev/null) || true

if [ -z "$FASTA_PATH" ]; then
    echo "No FASTA file selected. Exiting."
    exit 1
fi
echo "Reference FASTA: $FASTA_PATH"

# Pick the VCF file via Finder

echo ""
echo "Select VCF file"

VCF_PATH=$(osascript -e 'POSIX path of (choose file with prompt "Select your VCF file (.vcf or .vcf.gz):")' 2>/dev/null) || true

if [ -z "$VCF_PATH" ]; then
    echo "No VCF file selected. Exiting."
    exit 1
fi
echo "VCF file: $VCF_PATH"
echo ""

# Ask for a name/prefix for this graph

read -p "Enter a name/prefix for this graph (e.g. 'lct_index'): " INDEX_PREFIX
if [ -z "$INDEX_PREFIX" ]; then
    INDEX_PREFIX="pangenome_index"
    echo "No prefix entered, defaulting to '$INDEX_PREFIX'"
fi

# Set up a working directory

WORKDIR="$HOME/vg-builds/$INDEX_PREFIX"
mkdir -p "$WORKDIR"

FASTA_NAME=$(basename "$FASTA_PATH")
VCF_NAME=$(basename "$VCF_PATH")

cp "$FASTA_PATH" "$WORKDIR/$FASTA_NAME"
cp "$VCF_PATH" "$WORKDIR/$VCF_NAME"

echo ""
echo "Working directory: $WORKDIR"
cd "$WORKDIR"

# Check VCF is bgzip-compressed

if [[ "$VCF_NAME" != *.gz ]]; then
    echo "VCF is not compressed — compressing with bgzip..."
    docker run --rm -v "$(pwd):/data" quay.io/vgteam/vg bgzip -f "/data/$VCF_NAME"
    VCF_NAME="${VCF_NAME}.gz"
fi

# Build/rebuild the tabix index

echo "Indexing VCF with tabix..."
docker run --rm -v "$(pwd):/data" quay.io/vgteam/vg tabix -p vcf -f "/data/$VCF_NAME"

# Clean up any old index files with this prefix

rm -f ${INDEX_PREFIX}.*

# Build graph + Giraffe indexes

echo "Building pangenome graph with vg autoindex (this may take a while for large regions)..."
docker run --rm -v "$(pwd):/data" quay.io/vgteam/vg vg autoindex --workflow giraffe \
    -r "/data/$FASTA_NAME" \
    -v "/data/$VCF_NAME" \
    -p "/data/$INDEX_PREFIX"

# Show results

echo ""
echo "=== Done! Index files created in: $WORKDIR ==="
ls -la ${INDEX_PREFIX}*
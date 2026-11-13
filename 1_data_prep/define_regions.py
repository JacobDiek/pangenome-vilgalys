#!/usr/bin/env python3
"""Define study regions in T2T-CHM13 coordinates from anchor genes.

For each region in regions.tsv: find its anchor genes in the CHM13 RefSeq
annotation (GFF3), span from the leftmost gene start to the rightmost gene end,
add FLANK bp on both sides, and write one BED line (0-based, end-exclusive)
using the GenBank chromosome names (CP...) that the FASTA uses.

Usage:
  python3 define_regions.py --gff chm13.gff --fai chm13.fa.fai \\
      --regions regions.tsv --flank 100000 \\
      --annot-dir ANNOT_DIR --out-bed regions_chm13.bed [--strict]

Writes:
  out-bed                  one line per region: chrom, start, end, name
  ANNOT_DIR/chm13_genes.gff     the annotation's "gene" rows only
  ANNOT_DIR/nc_to_cp.tsv        RefSeq name -> GenBank name, matched by length
  ANNOT_DIR/gene_hits.tsv       region, gene, RefSeq chrom, start, end (1-based)
"""

import argparse
import sys
from collections import defaultdict
from pathlib import Path


def die(msg):
    sys.exit(f"ERROR: {msg}")


def read_regions(path):
    """regions.tsv -> {region: [gene, ...]}, in file order. '#' lines are comments."""
    regions = {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 2:
                die(f"{path}: expected 'region<TAB>genes', got: {line!r}")
            name = fields[0].strip()
            genes = [g.strip() for g in fields[1].split(",") if g.strip()]
            if not genes:
                die(f"{path}: region {name} has no anchor genes")
            if name in regions:
                die(f"{path}: region {name} is listed twice")
            regions[name] = genes
    if not regions:
        die(f"{path}: no regions found")
    return regions


def read_gff(gff_path, genes_out):
    """Scan the annotation once.

    Returns (refseq_lengths, genes):
      refseq_lengths: {RefSeq chrom: length}, from the 'region' rows that start at 1
      genes:          {gene Name: [(chrom, start, end), ...]}, 1-based coordinates
    Also writes the 'gene' rows to genes_out.
    """
    refseq_lengths = {}
    genes = defaultdict(list)
    with open(gff_path) as gff, open(genes_out, "w") as out:
        for line in gff:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            chrom, ftype, start, end, attributes = f[0], f[2], f[3], f[4], f[8]
            if ftype == "region" and start == "1":
                refseq_lengths[chrom] = int(end)
            elif ftype == "gene":
                out.write(line)
                # attributes look like: ID=gene-HOXA1;Dbxref=...;Name=HOXA1;...
                attrs = dict(
                    kv.split("=", 1) for kv in attributes.split(";") if "=" in kv
                )
                name = attrs.get("Name")
                if name:
                    genes[name].append((chrom, int(start), int(end)))
    return refseq_lengths, genes


def match_names(refseq_lengths, fai_path):
    """Match RefSeq names (NC_...) to GenBank names (CP...) by chromosome length.

    Returns (nc_to_cp, genbank_lengths). Stops if a length matches more than one
    chromosome, so a coincidence can never pair the wrong names silently.
    """
    genbank_lengths = {}
    with open(fai_path) as fh:
        for line in fh:
            name, length = line.split("\t")[:2]
            genbank_lengths[name] = int(length)

    by_length = defaultdict(list)
    for name, length in genbank_lengths.items():
        by_length[length].append(name)

    nc_to_cp = {}
    for nc, length in refseq_lengths.items():
        matches = by_length.get(length, [])
        if len(matches) > 1:
            die(
                f"{nc} ({length} bp) matches several FASTA sequences by length: {matches}"
            )
        if matches:
            nc_to_cp[nc] = matches[0]
    return nc_to_cp, genbank_lengths


def find_anchors(regions, genes, strict):
    """Look up every anchor gene. Returns [(region, gene, chrom, start, end), ...]."""
    hits = []
    for region, names in regions.items():
        for gene in names:
            entries = genes.get(gene, [])
            if not entries:
                die(f"gene {gene} (region {region}) not found in the annotation")
            if len(entries) > 1:
                listing = "; ".join(f"{c}:{s}-{e}" for c, s, e in entries)
                if strict:
                    die(f"gene {gene} has {len(entries)} entries: {listing}")
                print(
                    f"WARNING: gene {gene} has {len(entries)} entries, using the first: {listing}",
                    file=sys.stderr,
                )
            hits.append((region, gene, *entries[0]))
    return hits


def build_bed(hits, nc_to_cp, genbank_lengths, flank):
    """One BED row per region: leftmost gene start to rightmost gene end, plus flank.

    BED is 0-based and end-exclusive, so start = (1-based start - 1).
    """
    by_region = defaultdict(list)
    for region, gene, chrom, start, end in hits:
        by_region[region].append((gene, chrom, start, end))

    rows = []
    for region, anchors in by_region.items():
        chroms = {chrom for _, chrom, _, _ in anchors}
        if len(chroms) > 1:
            die(
                f"anchor genes for {region} are on different chromosomes: {sorted(chroms)}"
            )
        nc = chroms.pop()
        if nc not in nc_to_cp:
            die(f"no FASTA chromosome matches {nc} (needed for {region})")
        cp = nc_to_cp[nc]
        lo = min(start for _, _, start, _ in anchors)
        hi = max(end for _, _, _, end in anchors)
        bed_start = max(lo - 1 - flank, 0)
        bed_end = min(hi + flank, genbank_lengths[cp])  # never past the chromosome end
        rows.append((cp, bed_start, bed_end, region))
    rows.sort(key=lambda r: r[3])  # by region name, like the old version
    return rows


def main():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument(
        "--gff", required=True, type=Path, help="CHM13 RefSeq annotation (GFF3)"
    )
    p.add_argument(
        "--fai",
        required=True,
        type=Path,
        help="chm13.fa.fai (GenBank names and lengths)",
    )
    p.add_argument("--regions", required=True, type=Path, help="regions.tsv")
    p.add_argument(
        "--flank",
        type=int,
        default=100000,
        help="bp added on each side (default 100000)",
    )
    p.add_argument(
        "--annot-dir", required=True, type=Path, help="where the annotation tables go"
    )
    p.add_argument("--out-bed", required=True, type=Path, help="output BED file")
    p.add_argument(
        "--strict",
        action="store_true",
        help="stop if an anchor gene has several entries, instead of using the first",
    )
    args = p.parse_args()

    args.annot_dir.mkdir(parents=True, exist_ok=True)
    regions = read_regions(args.regions)

    refseq_lengths, genes = read_gff(args.gff, args.annot_dir / "chm13_genes.gff")
    nc_to_cp, genbank_lengths = match_names(refseq_lengths, args.fai)
    # 24 for CHM13: the RefSeq annotation has no row for the mitochondrion (CP068254.1)
    print(f"matched {len(nc_to_cp)} chromosome names (expect 24: 22 autosomes + X + Y)")
    with open(args.annot_dir / "nc_to_cp.tsv", "w") as out:
        out.writelines(f"{nc}\t{cp}\n" for nc, cp in nc_to_cp.items())

    hits = find_anchors(regions, genes, args.strict)
    with open(args.annot_dir / "gene_hits.tsv", "w") as out:
        out.writelines("\t".join(map(str, row)) + "\n" for row in hits)
    print("anchor genes (region, gene, RefSeq chrom, start, end):")
    for row in hits:
        print("  " + "\t".join(map(str, row)))

    rows = build_bed(hits, nc_to_cp, genbank_lengths, args.flank)
    args.out_bed.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out_bed, "w") as out:
        for row in rows:
            out.write("\t".join(map(str, row)) + "\n")
    print(f"wrote {len(rows)} regions to {args.out_bed}")


if __name__ == "__main__":
    main()

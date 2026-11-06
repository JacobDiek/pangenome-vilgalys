#!/usr/bin/env python3
"""
Find where each human (CHM13) region sits in a chimp haplotype, using only the
alignment blocks that touch the region's unique 100 kb flanks, then extract
the chimp sequences with samtools.

Python rewrite of steps 3-4 of align_regions.sh. Uses only the standard library.

Usage:
  python parse_paf.py --paf-dir ALN_DIR --asm-dir ASM_DIR --out-dir REGIONS_DIR \
                      --haps ptr_pri ptr_alt [--flank 100000] [--min-mapq 20]
"""

import argparse
import subprocess
from collections import defaultdict
from pathlib import Path


def parse_paf(paf_path, flank, min_mapq):
    """Tally anchor alignments per (region, chimp chromosome) pair."""
    # one record per (region, chimp chromosome); defaultdict creates it on first use
    stats = defaultdict(
        lambda: {
            "bp": 0,  # total matching bases (the "support")
            "lo": None,
            "hi": 0,  # leftmost start / rightmost end on the chimp chromosome
            "plus": 0,
            "minus": 0,  # matching bases on each strand
            "left": False,
            "right": False,  # did any block touch the left / right flank?
        }
    )
    region_len = {}

    with open(paf_path) as paf:
        for line in paf:
            col = line.rstrip("\n").split("\t")
            region, qlen, qstart, qend = col[0], int(col[1]), int(col[2]), int(col[3])
            strand, chrom = col[4], col[5]
            tstart, tend = int(col[7]), int(col[8])
            matches, mapq = int(col[9]), int(col[11])

            touches_left = qstart < flank
            touches_right = qend > qlen - flank
            if mapq < min_mapq or not (touches_left or touches_right):
                continue  # not a confident flank anchor

            s = stats[(region, chrom)]
            s["bp"] += matches
            s["lo"] = tstart if s["lo"] is None else min(s["lo"], tstart)
            s["hi"] = max(s["hi"], tend)
            s["plus" if strand == "+" else "minus"] += matches
            s["left"] = s["left"] or touches_left
            s["right"] = s["right"] or touches_right
            region_len[region] = qlen

    # for each region, keep the chimp chromosome with the most supporting bases
    best = {}
    for (region, chrom), s in stats.items():
        if region not in best or s["bp"] > best[region][1]["bp"]:
            best[region] = (chrom, s)

    results = []
    for region, (chrom, s) in best.items():
        if s["left"] and s["right"]:
            anchors = "both"
        elif s["left"]:
            anchors = "left_only"
        else:
            anchors = "right_only"
        results.append(
            {
                "chrom": chrom,
                "start": s["lo"],
                "end": s["hi"],
                "region": region,
                "strand": "+" if s["plus"] >= s["minus"] else "-",
                "chimp_len": s["hi"] - s["lo"],
                "human_len": region_len[region],
                "anchors": anchors,
            }
        )
    return results


def extract(fasta, r, out_path):
    """Cut the region out with samtools; reverse-complement if it's on the minus strand."""
    cmd = ["samtools", "faidx"]
    if r["strand"] == "-":
        cmd.append("-i")
    cmd += [
        str(fasta),
        f"{r['chrom']}:{r['start'] + 1}-{r['end']}",
    ]  # PAF 0-based -> samtools 1-based
    with open(out_path, "w") as out:
        subprocess.run(
            cmd, stdout=out, check=True
        )  # check=True: stop if samtools fails


def main():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--paf-dir", type=Path, required=True)
    p.add_argument("--asm-dir", type=Path, required=True)
    p.add_argument("--out-dir", type=Path, required=True)
    p.add_argument("--haps", nargs="+", required=True)
    p.add_argument("--flank", type=int, default=100_000)
    p.add_argument("--min-mapq", type=int, default=20)
    p.add_argument("--no-extract", action="store_true", help="only write the table")
    args = p.parse_args()

    rows = []
    for hap in args.haps:
        for r in parse_paf(
            args.paf_dir / f"query_vs_{hap}.paf", args.flank, args.min_mapq
        ):
            r["hap"] = hap
            rows.append(r)
    rows.sort(key=lambda r: (r["region"], r["hap"]))

    cols = [
        "chrom",
        "start",
        "end",
        "region",
        "hap",
        "strand",
        "chimp_len",
        "human_len",
        "anchors",
    ]
    table = args.out_dir / "regions_chimp.tsv"
    with open(table, "w") as out:
        out.write("chimp_chrom\t" + "\t".join(cols[1:]) + "\n")
        for r in rows:
            out.write("\t".join(str(r[c]) for c in cols) + "\n")

    print("chimp regions:")
    for r in rows:
        if not args.no_extract:
            extract(
                args.asm_dir / f"{r['hap']}.fa",
                r,
                args.out_dir / f"{r['region']}_{r['hap']}.fa",
            )
        ratio = r["chimp_len"] / r["human_len"]
        flag = "  <-- CHECK" if r["anchors"] != "both" or not 0.5 <= ratio <= 2 else ""
        print(
            f"  {r['region']} {r['hap']}  {r['chrom']}:{r['start'] + 1}-{r['end']} ({r['strand']})  "
            f"{r['chimp_len'] // 1000} kb vs human {r['human_len'] // 1000} kb  "
            f"ratio={ratio:.2f}  anchors={r['anchors']}{flag}"
        )


if __name__ == "__main__":
    main()

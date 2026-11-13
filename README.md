# pangenome-vilgalys

Scripts for the cross-species ape pangenome project (Vilgalys lab, UGA).
Measures reference bias when mapping non-human ape short reads, comparing
linear references (bwa) with a Minigraph-Cactus pangenome (vg giraffe).

## Layout (run in order)

Study regions are listed in `regions.tsv` (region name + anchor genes); add a line there to add a region. Shared paths and settings live in `config.sh`, which the scripts load. Scripts in `2_graphs/` and `3_mapping/` take a region and a direction:

- `REGION`: GYP, HOXA, LCT or MHC (any region with `{REGION}_{chm13,ptr_pri,ptr_alt}.fa` in `apes/regions/`)
- `pri`: graph = CHM13 + PTR.1 (ptr_pri), reads simulated from PTR.2 (ptr_alt)
- `alt`: graph = CHM13 + PTR.2 (ptr_alt), reads simulated from PTR.1 (ptr_pri)

| Folder | Script | Does |
|---|---|---|
| `1_data_prep/` | `download_pilot.sh` | Download CHM13 + chimp mPanTro3 (pri/alt) from NCBI, `samtools faidx` |
| | `define_regions.sh` + `define_regions.py` | Anchor genes from `regions.tsv` -> CHM13 annotation -> BED -> extract human regions (bash runs the tools, Python does the parsing) |
| | `align_regions_py.sh` | minimap2 human regions -> chimp (regions from `regions.tsv`), then calls `parse_paf.py` |
| | `parse_paf.py` | Parse PAF, pick best chimp locus per region, extract chimp FASTAs |
| `2_graphs/` | `build_graph.sh REGION DIR` | Minigraph-Cactus graph (CHM13 + in-graph chimp) + Giraffe zipcode index |
| | `index_linear_refs.sh REGION DIR` | CHM13 and in-graph chimp region FASTAs for bwa, + `bwa index` |
| `3_mapping/` | `simulate_reads.sh REGION DIR [COV]` | Simulate short reads from the held-out chimp haplotype (`vg sim`) |
| | `map_reads.sh REGION DIR [COV]` | giraffe (graph, surjected to CHM13) + bwa (CHM13, in-graph chimp); mapping-rate and MAPQ summary |
| `archive/` | | Superseded scripts: HOXA-only `build_hoxa_graph.sh` and `map_hoxa_reads.sh`, awk PAF parser, FASTA+VCF builder, Mac/Docker Graphviz renderer, Part 1 LCT scripts in `archive/lct/` |

Example, one region and direction end to end:

```bash
cd ~/pangenome/apes/scripts
sbatch -J build_graph_HOXA_pri 2_graphs/build_graph.sh HOXA pri
sbatch 2_graphs/index_linear_refs.sh HOXA pri
sbatch 3_mapping/simulate_reads.sh HOXA pri 30
sbatch 3_mapping/map_reads.sh HOXA pri 30     # after the three above finish
```

## Reference docs
- [`docs/FILE_TYPES.md`](docs/FILE_TYPES.md): every file type the pipeline makes or reads so far (FASTA to the final graph), with real examples, which program makes each, and how to look inside it.

# pangenome-vilgalys

Scripts for the cross-species ape pangenome project (Vilgalys lab, UGA).
Measures reference bias when mapping non-human ape short reads, comparing
linear references (bwa) with a Minigraph-Cactus pangenome (vg giraffe).

## Layout (run in order)

| Folder | Script | Does |
|---|---|---|
| `1_data_prep/` | `download_pilot.sh` | Download CHM13 + chimp mPanTro3 (pri/alt) from NCBI, `samtools faidx` |
| | `define_regions.sh` | Gene anchors from the CHM13 annotation -> BED -> extract human regions |
| | `align_regions_py.sh` | minimap2 human regions -> chimp, then calls `parse_paf.py` |
| | `parse_paf.py` | Parse PAF, pick best chimp locus per region, extract chimp FASTAs |
| `2_graphs/` | `build_hoxa_graph.sh` | Minigraph-Cactus (`cactus-pangenome`, container) graph of HOXA: CHM13 + PTR.1 |
| `3_mapping/` | `map_hoxa_reads.sh` | Simulate reads from held-out PTR.2, map with giraffe (graph) and bwa (CHM13, PTR.1), summarize mapping rate and MAPQ |
| `archive/` | | Superseded scripts: awk PAF parser, FASTA+VCF builder, Mac/Docker Graphviz renderer (`visualize_pangenome.sh`, template for a future container-based viz script), Part 1 LCT scripts in `archive/lct/` |

## Reference docs
- [`docs/FILE_TYPES.md`](docs/FILE_TYPES.md): every file type the pipeline makes or reads so far (FASTA to the final graph), with real examples, which program makes each, and how to look inside it.

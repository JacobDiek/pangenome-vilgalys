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
| `archive/` | | Superseded scripts: awk PAF parser, FASTA+VCF builder, Mac/Docker Graphviz renderer (`visualize_pangenome.sh`, template for a future container-based viz script), Part 1 LCT scripts in `archive/lct/` |

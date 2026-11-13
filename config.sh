# config.sh: settings shared by the pipeline scripts (1_data_prep/, 2_graphs/, 3_mapping/).
# Not run on its own. Each script loads it with:
#   source /home/jad33864/pangenome/apes/scripts/config.sh
# (A full path, because sbatch runs a COPY of the script from Slurm's own folder,
#  so a script can't find files "next to itself".)

# ---------- Locations ----------
PANGENOME=/home/jad33864/pangenome
APES=$PANGENOME/apes
SIF=$PANGENOME/containers/cactus_v2.9.9.sif
REGIONS_DIR=$APES/regions            # {REGION}_{chm13,ptr_pri,ptr_alt}.fa from 1_data_prep

# ---------- Study regions ----------
REGIONS_TSV=$APES/scripts/regions.tsv   # region name + anchor genes; add a line to add a region
FLANK=100000                            # bp added on each side of a region's anchor genes

# Cluster modules for bwa and samtools (vg comes from the container instead)
MAP_MODULES="BWA/0.7.18-GCCcore-13.3.0 SAMtools/1.23.1-GCC-13.3.0"

# ---------- Read simulation settings (shared so simulate and map agree on the run tag) ----------
READ_LEN=150
SEED=42               # fixed seed: the same reads every run
SUB_RATE=0.01         # substitution error rate
INDEL_RATE=0.002      # indel error rate
PAIRED=false          # true = paired-end (check that mate names match first)
FRAG_LEN=300          # paired-end only: mean fragment length
FRAG_SD=50            # paired-end only: fragment length std dev
DEFAULT_COVERAGE=30

# ---------- Helpers ----------
# Run vg inside the Cactus container: the same vg that builds the graph indexes.
# --cleanenv keeps loaded modules' settings from leaking into the container.
vg() { apptainer exec --cleanenv "$SIF" vg "$@"; }

# Stop with a clear message if any listed file is missing or empty.
need() {
    local f
    for f in "$@"; do
        [[ -s "$f" ]] || { echo "ERROR: missing input: $f" >&2; exit 1; }
    done
}

# set_region REGION DIRECTION
# Turns the two arguments into every name and path the scripts use, so all
# scripts agree on where things live. DIRECTION says which chimp haplotype is
# in the graph and which one is held out to simulate reads from:
#   pri: graph = CHM13 + PTR.1 (ptr_pri); reads from PTR.2 (ptr_alt)
#   alt: graph = CHM13 + PTR.2 (ptr_alt); reads from PTR.1 (ptr_pri)
set_region() {
    REGION=${1:-}
    DIR=${2:-}
    if [[ -z "$REGION" || -z "$DIR" ]]; then
        echo "Usage: <script> REGION pri|alt   (e.g. HOXA pri)" >&2
        exit 1
    fi
    case "$DIR" in
        pri) IN_HAP=ptr_pri; IN_SAMPLE=PTR.1; IN_TAG=PTR1; HELD_HAP=ptr_alt; HELD_TAG=PTR2 ;;
        alt) IN_HAP=ptr_alt; IN_SAMPLE=PTR.2; IN_TAG=PTR2; HELD_HAP=ptr_pri; HELD_TAG=PTR1 ;;
        *)   echo "ERROR: direction must be pri or alt, got '$DIR'" >&2; exit 1 ;;
    esac

    # Region FASTAs from 1_data_prep
    CHM13_FA=$REGIONS_DIR/${REGION}_chm13.fa
    IN_FA=$REGIONS_DIR/${REGION}_${IN_HAP}.fa         # chimp haplotype in the graph
    HELD_FA=$REGIONS_DIR/${REGION}_${HELD_HAP}.fa     # held-out chimp haplotype (read source)

    # Pangenome graph (2_graphs/build_graph.sh)
    GRAPH_NAME=${REGION}_${DIR}                       # e.g. HOXA_pri
    GRAPH_DIR=$APES/graphs/$GRAPH_NAME
    GRAPH_OUT=$GRAPH_DIR/out                          # Cactus output only: wiped on rebuild
    GBZ=$GRAPH_OUT/$GRAPH_NAME.gbz
    DIST=$GRAPH_OUT/$GRAPH_NAME.dist
    ZIDX=$GRAPH_DIR/giraffe_idx                       # zipcode minimizer index (outside out/)
    ZMIN=$ZIDX/$GRAPH_NAME.withzip.min
    ZIP=$ZIDX/$GRAPH_NAME.zipcodes

    # Linear references for bwa (2_graphs/index_linear_refs.sh)
    REFS=$APES/refs/$REGION
    CHM13_REF=$REFS/chm13.fa                          # header >REGION (matches the graph VCF)
    IN_REF=$REFS/${IN_HAP}.fa                         # header >PTR1_REGION or >PTR2_REGION

    # Simulated reads (3_mapping/simulate_reads.sh), one folder per held-out haplotype
    READS=$APES/reads/${REGION}_${HELD_HAP}

    # Mapping results (3_mapping/map_reads.sh)
    RESULTS=$APES/results/$GRAPH_NAME
}

# set_run [COVERAGE]
# Run tag shared by simulate_reads.sh and map_reads.sh, e.g. cov30_s42 (or cov30_s42_pe).
set_run() {
    COVERAGE=${1:-$DEFAULT_COVERAGE}
    RUN=cov${COVERAGE}_s${SEED}
    if [[ "$PAIRED" == true ]]; then RUN=${RUN}_pe; fi
    SIM_GAM=$READS/$RUN.gam
    SIM_FQ=$READS/$RUN.fq
    RUN_OUT=$RESULTS/$RUN
}

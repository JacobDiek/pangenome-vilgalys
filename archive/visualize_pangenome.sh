#!/bin/bash
# visualize_pangenome.sh
#
# Takes an EXISTING pangenome graph and renders it as a PNG image using 
# Graphviz
#
# Requires: Docker Desktop running, quay.io/vgteam/vg image pulled,
# Graphviz installed locally

set -e

echo "=== Pangenome Graph Visualizer ==="
echo ""

# Make sure Graphviz's `dot` is installed locally

if ! command -v dot &> /dev/null; then
    echo "Graphviz's 'dot' command was not found."
    exit 1
fi

# Pick the graph file

GRAPH_PATH=$(osascript -e 'POSIX path of (choose file with prompt "Select your pangenome graph (.gbz):")' 2>/dev/null) || true

if [ -z "$GRAPH_PATH" ]; then
    echo "No graph file selected. Exiting."
    exit 1
fi
echo "Graph file: $GRAPH_PATH"
echo ""

# Split into directory + filename

GRAPH_DIR=$(dirname "$GRAPH_PATH")
GRAPH_NAME=$(basename "$GRAPH_PATH")
cd "$GRAPH_DIR"

# Full or subgraph

echo "Enter a region to visualize, as  contig:start-end  (e.g. chr2:1-2000)"
read -p "Leave blank to render the FULL graph" REGION

if [ -n "$REGION" ]; then
    # Subregion
    read -p "Context to include around the region, in bp (default 10): " CONTEXT
    CONTEXT="${CONTEXT:-10}"

    echo "Extracting subgraph for $REGION (+/- $CONTEXT bp context)..."
    docker run --rm -v "$(pwd):/data" quay.io/vgteam/vg vg find \
        -x "/data/$GRAPH_NAME" -p "$REGION" -c "$CONTEXT" > subgraph.vg

    DOT_FILE="subgraph_${REGION//[:\-]/_}.dot"
    PNG_FILE="subgraph_${REGION//[:\-]/_}.png"

    docker run --rm -v "$(pwd):/data" quay.io/vgteam/vg vg view -d "/data/subgraph.vg" > "$DOT_FILE"
    rm -f subgraph.vg
else
    # Full graph mode
    echo "Rendering the FULL graph — this can be slow or unreadable for large graphs."
    DOT_FILE="${GRAPH_NAME%.giraffe.gbz}_full.dot"
    PNG_FILE="${GRAPH_NAME%.giraffe.gbz}_full.png"

    docker run --rm -v "$(pwd):/data" quay.io/vgteam/vg vg view -d "/data/$GRAPH_NAME" > "$DOT_FILE"
fi

# Warn if the DOT file is large before attempting to render
LINE_COUNT=$(wc -l < "$DOT_FILE" | tr -d ' ')
echo "Graph description has $LINE_COUNT lines."

if [ "$LINE_COUNT" -gt 3000 ]; then
    echo ""
    echo "WARNING: This graph is large ($LINE_COUNT lines) — rendering may take"
    echo "a long time and the resulting image may be too dense to read."
    read -p "Continue rendering anyway? (y/n): " CONFIRM
    if [[ "$CONFIRM" != "y" && "$CONFIRM" != "Y" ]]; then
        echo "Skipped rendering. DOT file was still saved at: $GRAPH_DIR/$DOT_FILE"
        exit 0
    fi
fi

# Render and open the image

echo "Rendering image with Graphviz..."
dot -Tpng "$DOT_FILE" -o "$PNG_FILE"

echo ""
echo "=== Done! Image saved to: $GRAPH_DIR/$PNG_FILE ==="
open "$PNG_FILE"
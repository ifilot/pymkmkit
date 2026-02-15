#!/usr/bin/env bash

set -e

STATES_BASE="/mnt/d/data/Ru1121_dataset/states"
TS_BASE="/mnt/d/data/Ru1121_dataset/ts"

ISFS_OUT="ISFS"
TS_OUT="TS"

mkdir -p "$ISFS_OUT"
mkdir -p "$TS_OUT"

echo "================================"
echo "Processing state minima"
echo "================================"

for dir in "$STATES_BASE"/*; do
    name=$(basename "$dir")

    # skip empty surface for now
    if [[ "$name" == "empty" ]]; then
        echo "Skipping empty state (handled separately)"
        continue
    fi

    outcar="$dir/OUTCAR"

    if [[ ! -f "$outcar" ]]; then
        echo "No OUTCAR found in $name — skipping"
        continue
    fi

    lname=$(echo "$name" | tr '[:upper:]' '[:lower:]')

    outfile="$ISFS_OUT/${lname}.yaml"

    echo "State: $name → $outfile"

    pymkmkit freq2yaml "$outcar" \
        -o "$outfile" \
        --average-pairs
done


echo ""
echo "================================"
echo "Processing transition states"
echo "================================"

for dir in "$TS_BASE"/*; do
    name=$(basename "$dir")

    outcar="$dir/OUTCAR"

    if [[ ! -f "$outcar" ]]; then
        echo "No OUTCAR found in TS $name — skipping"
        continue
    fi

    lname=$(echo "$name" | tr '[:upper:]' '[:lower:]')

    outfile="$TS_OUT/${lname}.yaml"

    echo "TS: $name → $outfile"

    pymkmkit freq2yaml "$outcar" \
        -o "$outfile" \
        --average-pairs
done


echo ""
echo "Dataset build complete."

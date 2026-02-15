#!/usr/bin/env bash

set -e

STATES_BASE="/mnt/d/data/Ru1121_dataset/states"
TS_BASE="/mnt/d/data/Ru1121_dataset/ts"
GAS_BASE="/mnt/d/data/Ru1121_dataset/gas"

ISFS_OUT="ISFS"
TS_OUT="TS"
GAS_OUT="GAS"

mkdir -p "$ISFS_OUT"
mkdir -p "$TS_OUT"
mkdir -p "$GAS_OUT"

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

    if [[ -f "$outfile" ]]; then
        echo "State: $name → output exists, skipping"
        continue
    fi

    echo "State: $name → $outfile"

    pymkmkit freq2yaml "$outcar" \
        -o "$outfile" \
        --average-pairs
done


echo ""
echo "================================"
echo "Processing empty surface"
echo "================================"

EMPTY_DIR="$STATES_BASE/empty"
EMPTY_OUT="$ISFS_OUT/empty.yaml"

if [[ -d "$EMPTY_DIR" ]]; then
    outcar="$EMPTY_DIR/OUTCAR"

    if [[ ! -f "$outcar" ]]; then
        echo "No OUTCAR found in empty — skipping"
    elif [[ -f "$EMPTY_OUT" ]]; then
        echo "Empty surface → output exists, skipping"
    else
        echo "Empty surface → $EMPTY_OUT"

        pymkmkit opt2yaml "$outcar" \
            -o "$EMPTY_OUT"
    fi
else
    echo "Empty surface directory not found — skipping"
fi


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

    if [[ -f "$outfile" ]]; then
        echo "TS: $name → output exists, skipping"
        continue
    fi

    echo "TS: $name → $outfile"

    pymkmkit freq2yaml "$outcar" \
        -o "$outfile" \
        --average-pairs
done


echo ""
echo "================================"
echo "Processing gas phase states"
echo "================================"

for dir in "$GAS_BASE"/*; do
    name=$(basename "$dir")

    outcar="$dir/OUTCAR"

    if [[ ! -f "$outcar" ]]; then
        echo "No OUTCAR found in GAS $name — skipping"
        continue
    fi

    lname=$(echo "$name" | tr '[:upper:]' '[:lower:]')
    outfile="$GAS_OUT/${lname}.yaml"

    if [[ -f "$outfile" ]]; then
        echo "GAS: $name → output exists, skipping"
        continue
    fi

    echo "GAS: $name → $outfile"

    pymkmkit freq2yaml "$outcar" \
        -o "$outfile"
done


echo ""
echo "Dataset build complete."

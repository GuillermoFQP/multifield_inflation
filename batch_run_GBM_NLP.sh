#!/usr/bin/env bash

set -euo pipefail

PROG="ps_full"
OUTDIR="ps_full_GBM_NLP_amp"

# Build
make "$PROG"

# Prepare output directory
mkdir -p "$OUTDIR"

# Run jobs
i=1
for val in $(seq 0.05 0.05 1.00); do
    printf -v idx "%d" "$i"
    printf -v val_fmt "%.2f" "$val"

    "./$PROG" "$val_fmt" > "$OUTDIR/${OUTDIR}${idx}.txt"

    ((i++))
done

echo "All jobs finished."

#!/bin/bash

RUNME="./runme.sh"
PROG="./efold_cmap.out"

# 1. Build efold_cmap.out using runme.sh
$RUNME build3 || { echo "Build failed"; exit 1; }
$RUNME clean

# 2. Prepare output directory
mkdir -p Frames

# 3. Job settings
i=1

# 4. Loop over parameter range
for exp in $(seq -6.0 0.1 6.0); do
	echo "Computing e-fold colormap for M = exp($exp)"
    M=$(echo "e($exp)" | bc -l)
    $PROG "$M" > "Frames/frame_$i"
    i=$((i + 1))
done

wait  # Wait for remaining background jobs
$RUNME cleanall
echo "All jobs finished."

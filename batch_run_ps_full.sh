#!/bin/bash

RUNME="./runme.sh"
PROG="./ps_full.out"

# 1. Build ps_mode.out using runme.sh
$RUNME build1 || { echo "Build failed"; exit 1; }
$RUNME clean

# 2. Prepare output directory
mkdir -p Frames

# 3. Job settings
i=1

# 4. Loop over parameter range
for exp in $(seq 5.0 0.01 9.0); do
	echo "Computing power spectrum for M = exp($exp)"
    M=$(echo "e($exp)" | bc -l)
    $PROG "$M" > "Frames/frame_$i"
    i=$((i + 1))
done

wait  # Wait for remaining background jobs
$RUNME cleanall
echo "All jobs finished."

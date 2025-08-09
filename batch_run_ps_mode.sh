#!/bin/bash

RUNME="./runme.sh"
PROG="./ps_mode.out"

# 1. Build ps_mode.out using runme.sh
$RUNME build2 || { echo "Build failed"; exit 1; }
$RUNME clean

# 2. Prepare output directory
mkdir -p Frames

# 3. Job settings
MAX_JOBS=8  # Number of concurrent jobs
i=1

# 4. Loop over parameter range
for exp in $(seq -2.0 0.1 8.0); do
    echo "Computing power spectrum for M = exp($exp)"
    M=$(echo "e($exp)" | bc -l)   # Calculate M
    $PROG "$M" > "Frames/frame_$i" &

    # Limit concurrent jobs
    if (( i % MAX_JOBS == 0 )); then
        wait
    fi

    i=$((i + 1))
done

wait  # Wait for remaining background jobs
$RUNME cleanall
echo "All jobs finished."

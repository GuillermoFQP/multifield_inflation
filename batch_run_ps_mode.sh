#!/bin/bash

PROG="ps_mode"

# 1. Build ps_mode.out using runme.sh
make $PROG || { echo "Build failed"; exit 1; }
make clean

# 2. Prepare output directory
mkdir -p Frames

# 3. Job settings
MAX_JOBS=8  # Number of concurrent jobs
i=1

# 4. Loop over parameter range
for exp in $(seq -2.0 0.1 8.0); do
    echo "Computing power spectrum for M = exp($exp)"
    M=$(echo "e($exp)" | bc -l)
    time ./$PROG "$M" > "Frames/frame_$i" &
    echo
    # Limit concurrent jobs
    if (( i % MAX_JOBS == 0 )); then
        wait
    fi
    i=$((i + 1))
done

make cleanall
echo "All jobs finished."

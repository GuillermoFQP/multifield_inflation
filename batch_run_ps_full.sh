#!/bin/bash

PROG="ps_full"

# 1. Build ps_mode.out using runme.sh
make $PROG || { echo "Build failed"; exit 1; }
make clean

# 2. Prepare output directory
mkdir -p Frames

# 3. Job settings
i=1

# 4. Loop over parameter range
for exp in $(seq 5.0 0.01 9.0); do
	echo "Computing power spectrum for M = exp($exp)"
    M=$(echo "e($exp)" | bc -l)
    time ./$PROG "$M" > "Frames/frame_$i"
    echo
    i=$((i + 1))
done

make cleanall
echo "All jobs finished."

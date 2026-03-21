#!/bin/bash

PROG="ps_full"

# 1. Build ps_mode.out using runme.sh
make $PROG || { echo "Build failed"; exit 1; }
make clean

# 2. Prepare output directory
mkdir -p ps_full_RPM_ELP15_frames

# 3. Job settings
i=1

# 4. Loop over parameter range
for exp in $(seq -3.0 0.5 3.0); do
	echo "Computing power spectrum for M = exp($exp)"
    M=$(echo "e($exp)" | bc -l)
    time ./$PROG "$M" > "ps_full_RPM_ELP15_frames/frames_$i.txt"
    echo
    i=$((i + 1))
done

make cleanall
echo "All jobs finished."

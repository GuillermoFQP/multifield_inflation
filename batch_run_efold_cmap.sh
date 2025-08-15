#!/bin/bash

PROG="efold_cmap"

# 1. Build efold_cmap.out using runme.sh
make $PROG || { echo "Build failed"; exit 1; }
make clean

# 2. Prepare output directory
mkdir -p Frames

# 3. Job settings
i=1

# 4. Loop over parameter range
for exp in $(seq -6.00 0.01 3.00); do
	echo "Computing e-fold colormap for M = exp($exp)"
    M=$(echo "e($exp)" | bc -l)
    time ./$PROG "$M" > "Frames/frame_$i"
    i=$((i + 1))
done

make cleanall
echo "All jobs finished."

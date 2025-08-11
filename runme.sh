#!/bin/bash

# ===== Configuration =====
FC="gfortran"
FLAGS="-fopenmp -O3 -fdefault-real-8"
FLAGS2="-fopenmp -O3"
MODULES="multifield_globals.f90 multifield_utils.f90"

# ===== Functions =====
help_menu() {
    echo "Usage: ./runme.sh [command] [argument]"
    echo
    echo "Commands:"
    echo "  build           Compile *.f90 and  -> *.out"
    echo "  build1          Compile ps_full.f90 and  -> ps_full.out"
    echo "  build2          Compile ps_mode.f90 and  -> ps_mode.out"
    echo "  build3          Compile efold_cmap.f90 and  -> efold_cmap.out"
    echo "  run1 [argument] Run ps_full.out with optional arguments"
    echo "                  Arguments must be real numbers"
    echo "  run2 [argument] Run ps_mode.out with optional arguments"
	echo "                  Arguments must be real numbers"
    echo "  run3 [argument] Run efold_cmap.out with optional arguments"
	echo "                  Arguments must be real numbers"
    echo "  clean           Remove .o and .mod files"
    echo "  cleanout        Remove .out executables"
    echo "  cleanall        Remove .o, .mod, and .out files"
    echo "  help            Show this menu"
    echo
}

build_prog() {
	build_prog1
	build_prog2
	build_prog3
}

build_prog1() {
    echo "Compiling ps_full.f90..."
    $FC $FLAGS $MODULES ps_full.f90 -o ps_full.out
}

build_prog2() {
    echo "Compiling ps_mode.f90..."
    $FC $FLAGS $MODULES ps_mode.f90 -o ps_mode.out
}

build_prog3() {
    echo "Compiling efold_cmap.f90..."
    $FC $FLAGS2 $MODULES efold_cmap.f90 -o efold_cmap.out
}

run_prog1() {
    if [[ ! -f ps_full.out ]]; then
        build_prog1
    fi
    echo "Running ps_full.out $@"
    ./ps_full.out "$@"
}

run_prog2() {
    if [[ ! -f ps_mode.out ]]; then
        build_prog2
    fi
    echo "Running ps_mode.out $@"
    ./ps_mode.out "$@"
}

run_prog3() {
    if [[ ! -f efold_cmap.out ]]; then
        build_prog3
    fi
    echo "Running efold_cmap.out $@"
    ./efold_cmap.out "$@"
}

clean() {
    echo "Removing .o and .mod files..."
    rm -f *.o *.mod
}

cleanout() {
    echo "Removing executables..."
    rm -f *.out
}

cleanall() {
    clean
    cleanout
}

# ===== Main =====
case "$1" in
    build)    build_prog ;;
    build1)   build_prog1 ;;
    build2)   build_prog2 ;;
    build3)   build_prog3 ;;
    run1)     shift; run_prog1 "$@" ;;
    run2)     shift; run_prog2 "$@" ;;
    run3)     shift; run_prog3 "$@" ;;
    clean)    clean ;;
    cleanout) cleanout ;;
    cleanall) cleanall ;;
    help|"")  help_menu ;;
    *)        echo "Unknown command: $1"; help_menu ;;
esac

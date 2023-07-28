#!/bin/bash

# This script runs the whole workflow for the DDEC charges calculation
source /dccstor/nanopore-2945/cp2k/cp2k-v2023.1/tools/toolchain/install/setup

CP2K_DIR=/dccstor/nanopore-2945/cp2k/cp2k-v2023.1/exe/local

export OMP_NUM_THREADS=1

echo -e "Running CP2K simulation..."
mpirun -np 36 $CP2K_DIR/cp2k.psmp -i simulation_SCF.inp -o simulation_SCF.out

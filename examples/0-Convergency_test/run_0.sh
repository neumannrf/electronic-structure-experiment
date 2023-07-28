#!/bin/bash

# Define environmenta variables
export PATH=$PATH:/home/felipelopes/Simulations/PRs/electronic-density-experiment/bin/
export EPATH=/home/felipelopes/Simulations/PRs/electronic-density-experiment/bin/
export CP2K_DIR=/home/felipelopes/Programas/cp2k
export CP2K_DATA_DIR=${CP2K_DIR}/data
export CHARGEMOL_DIR=/home/felipelopes/Programas/Chargemol/
export CHARGEMOL_DATA_FOLDER=${CHARGEMOL_DIR}/atomic_densities/

# Define specific variables
FrameworkName='MgMOF-74'
OutputFolder=$PWD
NumProcs=4

# Load the CP2K setup script
source ${CP2K_DIR}/tools/toolchain/install/setup

# First ensure that the CIF file has all atoms explicitly defined in a P1 cell
echo -e "\nCreating RASPA P1 cell input file..."
create_supercell.py --FrameworkName ${FrameworkName} \
                    --FrameworkFolder ${OutputFolder} \
                    --UnitCells '1,1,1' \
                    ${OutputFolder}

echo -e "\nRunning RASPA P1 cell simulation..."
cd ${OutputFolder} && $RASPA_DIR/bin/simulate -i simulation-CreateSupercell.input

echo -e "\nCopying P1-symmetry CIF file..."
mv -v Movies/System_0/Framework_0_final_*_P1.cif ${FrameworkName}.cif

rm -r Movies/ Output/ Restart/ VTK/



for PW in 200 250 300 350 400 450 500 550 600 650 700 750 800 850 900 950 1000 1050 1100 1150 1200 1250 1300 1350 1400 1500 1600 1700 1800 1900 2000 2500 3000 3500 4000; do

mkdir $PW

cd $PW

cp ../$FrameworkName.cif .
cp ../run.sh .

echo -e "\nCreating CP2K input file..."
charge_density.py --FrameworkName ${FrameworkName} \
		  --SCFGuess 'atomic' \
		  --CP2KDataDir '/dccstor/nanopore-2945/cp2k/cp2k-v2023.1/data' \
		  --UseSmearing  \
		  --AddedMOs  80 \
		  --PWCutoff $PW \
                  ${OutputFolder}/$PW
cd ..
done

# Electronic Structure Experiment

Automated virtual experiment that uses the CP2K software package to run DFT simulations on nanoporous materials.

It expects two environmental variables to be set:

- `CP2K_DATA_DIR`: Path to the CP2K directory containing basis set and pseudopotential files.
- `CHARGEMOL_DATA_FOLDER`: Path to the chargemol folder containing the atomic density files.

## Examples

### DDEC partial charges calculation

The script below shows a simple example of how to use the scripts provided in this repository to calculate the DDEC charges of a given cif file. To use it is necessary to set the `FrameworkName`, and `NProcs` variables. The script also expect to save all the output files on the same folder where the script is located. You can change it by setting the `OutputFolder` variable accordingly.

```bash
#!/bin/bash

# Define environment variables
export PATH=$PATH:/home/felipelopes/Simulations/electronic-density-experiment/bin/
export EPATH=/home/felipelopes/Simulations/electronic-density-experiment/bin/
export CP2K_DIR=/home/felipelopes/Programs/cp2k
export CP2K_DATA_DIR=${CP2K_DIR}/data
export CHARGEMOL_DIR=/home/felipelopes/Programs/Chargemol/
export CHARGEMOL_DATA_FOLDER=${CHARGEMOL_DIR}/atomic_densities/

# Define specific variables
FrameworkName='CIP-Me'
OutputFolder=$PWD
NProcs=4

# Load the CP2K setup script
source ${CP2K_DIR}/tools/toolchain/install/setup

# First ensure that the CIF file has all atoms explicitly defined in a P1 cell
echo -e "\nCreating primitive P1 cell..."
pmg structure --convert --filename ${OutputFolder}/${FrameworkName}.cif ${OutputFolder}/${FrameworkName}_prim.cif
mv -v ${OutputFolder}/${FrameworkName}_prim.cif ${OutputFolder}/${FrameworkName}.cif

echo -e "\nCreating CP2K input file..."
charge_density.py --FrameworkName ${FrameworkName} \
                  --SCFGuess 'restart' \
                  ${OutputFolder}

export OMP_NUM_THREADS=1
echo -e "\nRunning CP2K simulation with ${NProcs} MPI and ${OMP_NUM_THREADS} OMP process..."
mpirun -np ${NProcs} $CP2K_DIR/exe/local/cp2k.psmp -i simulate_SCF.inp -o simulate_SCF.out

echo -e "\nRenaming the cube file..."
mv *-valence_density-ELECTRON_DENSITY-1_0.cube valence_density.cube

# Make the number of OMP processes to ${NProcs}
export OMP_NUM_THREADS=${NProcs}
echo -e "\nCreating Chargemol input file..."
chargemol.py --FrameworkName ${FrameworkName} ${OutputFolder}

echo -e "\nRunning chargemol on ${OMP_NUM_THREADS} OMP processes..."
${CHARGEMOL_DIR}/Chargemol_09_26_2017_linux_parallel

echo -e "\nParsing the chargemol output file and creating the cif file with the DDEC charges..."
parse_charges.py --FrameworkName ${FrameworkName} ${OutputFolder}
```

### Geometry optimization

Below there is an example of a geometry optimization script. It is very similar to the previous one, but it uses the `structure_optimization.py` script to create the CP2K input file. It also uses the `parse_optimization.py` script to parse the CP2K output file and create a new cif file with the optimized structure.

```bash
#!/bin/bash

# Define environment variables
export PATH=$PATH:/home/felipelopes/Simulations/electronic-density-experiment/bin/
export EPATH=/home/felipelopes/Simulations/electronic-density-experiment/bin/
export CP2K_DIR=/home/felipelopes/Programs/cp2k
export CP2K_DATA_DIR=${CP2K_DIR}/data
export CHARGEMOL_DIR=/home/felipelopes/Programs/Chargemol/
export CHARGEMOL_DATA_FOLDER=${CHARGEMOL_DIR}/atomic_densities/

# Define specific variables
FrameworkName='CIP-Me'
OutputFolder=$PWD
NumProcs=4

# Load the CP2K setup script
source ${CP2K_DIR}/tools/toolchain/install/setup

# First ensure that the CIF file has all atoms explicitly defined in a P1 cell
echo -e "\nCreating primitive P1 cell..."
pmg structure --convert --filename ${OutputFolder}/${FrameworkName}.cif ${OutputFolder}/${FrameworkName}_prim.cif
mv -v ${OutputFolder}/${FrameworkName}_prim.cif ${OutputFolder}/${FrameworkName}.cif

echo -e "\nCreating CP2K input file..."
structure_optimization.py --FrameworkName ${FrameworkName} \
                          ${OutputFolder}

export OMP_NUM_THREADS=1
echo -e "\nRunning CP2K simulation with ${NumProcs} MPI and ${OMP_NUM_THREADS} OMP process..."
mpirun -np ${NumProcs} $CP2K_DIR/exe/local/cp2k.psmp -i simulation_Optimization.inp -o simulation_Optimization.out

parse_optimization.py --FrameworkName ${FrameworkName} \
                      --SaveHistory  \
                      ${OutputFolder}
```

### Vibrational frequencies calculation

Below there is an example of a vibrational frequencies script. It is very similar to the previous one, but it uses the `raman_ir.py` script to create the CP2K input file. It also uses the `parse_vibrations.py` script to parse the CP2K output file and create a new files with the results.

```bash
#!/bin/bash

# Define environment variables
export PATH=$PATH:/home/felipelopes/Simulations/electronic-density-experiment/bin/
export CP2K_DIR=/home/felipelopes/Programs/cp2k
export CP2K_DATA_DIR=${CP2K_DIR}/data

# Define specific variables
FrameworkName='MgMOF74'
OutputFolder=$PWD
NumProcs=32
export OMP_NUM_THREADS=1

# Load the CP2K setup script
source ${CP2K_DIR}/tools/toolchain/install/setup

# First ensure that the CIF file has all atoms explicitly defined in a P1 cell
echo -e "\nCreating primitive P1 cell..."
pmg structure --convert --filename ${OutputFolder}/${FrameworkName}.cif ${OutputFolder}/${FrameworkName}_prim.cif
mv -v ${OutputFolder}/${FrameworkName}_prim.cif ${OutputFolder}/${FrameworkName}.cif

echo -e "\nCreating input file for CP2K Vibrations calculation..."
raman_ir.py --FrameworkName ${FrameworkName} \
            --MaxSCFycles 50 \
            --ProcsPerReplica 16 \
            --dX 0.001 \
            --CalculateRaman \
            --CalculateIR \
            ${OutputFolder}

echo -e "\nRunning CP2K simulation with ${NumProcs} MPI and ${OMP_NUM_THREADS} OMP process..."
mpirun -np ${NumProcs} $CP2K_DIR/exe/local/cp2k.psmp -i simulation_Vibrations.inp -o simulation_Vibrations.out

echo -e "\nCreating RASPA P1 cell input file..."
parse_vibrations.py --FrameworkName ${FrameworkName} \
                    --SaveVibrations  \
                    --HalfWidth  5.0 \
                    ${OutputFolder}
```

The `parse_vibrations.py` creates two csv files with the results. The `{FrameworkName}_RAMAN_IR_Curve.csv` that holds the Raman and/or the IR curve created by adjusting a Lorentzian function to the calculated frequencies. The `{FrameworkName}_VibrationsTable.csv` that holds the calculated frequencies and intensities. The `--HalfWidth` option defines the half width of the Lorentzian function used to create the Raman and/or IR curve.

The `--SaveVibrations` option creates a folder called `Vibrations` with `AXSF` files containing the displacements of the atoms for each vibration mode. It is possible to visualize the vibrations using the [VESTA](https://jp-minerals.org/vesta/en/) software.

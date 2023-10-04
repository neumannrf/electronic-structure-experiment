# Electronic Structure Experiment Tutorials

The Electronic Structure Experiment (ESE) package currently supports the following simulations:

- Charge density calculation:
  - Band gap calculation
  - HOMO and LUMO orbital plots
- Partial charges calculation:
  - CM5 charges
  - DDEC charges
- Geometry optimization
- Vibrational frequencies calculation
  - FTIR spectra
  - Raman spectra

## Environment variables

To run the scripts properly, you must set the following environmental variables:

- `PATH`: The path to the bin folder of this repository.
- `CP2K_DATA_DIR`: The path to the CP2K directory containing basis set and pseudopotential files.
- `CHARGEMOL_DATA_FOLDER`: The path to the chargemol folder containing the atomic density files.
- `FrameworkName`: The name of the framework to be used in the calculations.
- `OutputFolder`: The path where the output files will be saved.
- `NProcs`: The number of MPI processes to be used in the calculations.

By default, the number of OpenMP processes will be set to 1. You can change this by setting the `OMP_NUM_THREADS` variable, but it is not recommended to use more than one OpenMP process per MPI process, as hybrid parallelization may not work properly.

Here is an example of how to set the environmental variables and run the scripts:

```bash
#!/bin/bash

# Define global environment variables
export PATH=$PATH:/dccstor/nanopore-2945/electronic-density-experiment/bin/

export CP2K_DIR=/dccstor/nanopore-2945/cp2k/cp2k-v2023.1
export CP2K_DATA_DIR=${CP2K_DIR}/data

export CHARGEMOL_DIR=/dccstor/nanopore-2945/Chargemol/
export CHARGEMOL_DATA_FOLDER=${CHARGEMOL_DIR}/atomic_densities/

export OMP_NUM_THREADS=1

# Define specific variables
FrameworkName='MgMOF-74'
OutputFolder=$PWD
NProcs=36

# Load the CP2K setup script
source ${CP2K_DIR}/tools/toolchain/install/setup
```

This header can be used in almost all scripts provided in this repository. You should change the specific variables accordingly. If you are not interested in calculating the partial charges, the `CHARGEMOL_DIR` and `CHARGEMOL_DATA_FOLDER` variables are not necessary.

## Charge density calculation

### Partial Charges

The calculation of partial charges can be divided into two main steps:

- Generating a charge density file (`.cube`): The `charge_density.py` script is used to generate the input for CP2K simulation. This can be done using the `--WriteDensity` option.

  - Only the `PBE` functional is supported for the charge density calculation.
  - The `DFTD3(BJ)` and `DFTD3` dispersion corrections are supported.
  - The `DZVP` and `TZV2P` basis sets are supported.

- Calculating the partial charges using Chargemol: The [Chargemol](https://sourceforge.net/projects/ddec/files/?source=navbar) software is used to calculate the partial charges.
  - `CM5` and `DDEC` schemes are supported. The `DDEC` is recomended.

The parse_charges.py script can be used to parse the Chargemol output file and create a new cif file with the partial charges.

**Note:** The `OMP_NUM_THREADS` variable is set to `NProcs` because Chargemol supports OpenMP parallelization but not MPI parallelization.

```bash
#!/bin/bash

# Define global environment variables
export PATH=$PATH:/dccstor/nanopore-2945/electronic-density-experiment/bin/

export CP2K_DIR=/dccstor/nanopore-2945/cp2k/cp2k-v2023.1
export CP2K_DATA_DIR=${CP2K_DIR}/data

export CHARGEMOL_DIR=/dccstor/nanopore-2945/Chargemol/
export CHARGEMOL_DATA_FOLDER=${CHARGEMOL_DIR}/atomic_densities/

export OMP_NUM_THREADS=1

# Define specific variables
FrameworkName='MgMOF-74'
OutputFolder=$PWD
NProcs=36

# Load the CP2K setup script
source ${CP2K_DIR}/tools/toolchain/install/setup

# First ensure that the CIF file has all atoms explicitly defined in a P1 cell
echo -e "\nCreating primitive P1 cell..."
pmg structure --convert --filename ${OutputFolder}/${FrameworkName}.cif ${OutputFolder}/${FrameworkName}_prim.cif
mv -v ${OutputFolder}/${FrameworkName}_prim.cif ${OutputFolder}/${FrameworkName}.cif

echo -e "\nCreating CP2K input file..."
charge_density.py --FrameworkName ${FrameworkName} \
                  --Functional "PBE" \
                  --DispersionCorrection "DFTD3(BJ)" \
                  --PWCutoff 1200 \
                  --BasisSet "TZV2P" \
                  --UseScalapack  \
                  --WriteDensity \
                  ${OutputFolder}

echo -e "\nRunning CP2K simulation with ${NProcs} MPI and ${OMP_NUM_THREADS} OMP process..."
mpirun -np ${NProcs} $CP2K_DIR/exe/local/cp2k.psmp -i simulation_SCF.inp -o simulation_SCF.out

echo -e "\nRenaming the cube file..."
mv *-valence_density-ELECTRON_DENSITY-1_0.cube valence_density.cube

echo -e "\nCreating Chargemol input file..."
chargemol.py --FrameworkName ${FrameworkName} ${OutputFolder}

# Make the number of OMP processes to ${NProcs}
export OMP_NUM_THREADS=${NProcs}

echo -e "\nRunning chargemol on ${OMP_NUM_THREADS} OMP processes..."
${CHARGEMOL_DIR}/Chargemol_09_26_2017_linux_parallel

echo -e "\nParsing the chargemol output file and creating the cif file with the DDEC charges..."
parse_charges.py --FrameworkName ${FrameworkName} \
                 --CM5 \
                 ${OutputFolder}
```

By default the DDEC charges are always calculated. The `--CM5` option can be used to calculate the CM5 charges. At the end the files `{FrameworkName}_DDEC.cif` and `{FrameworkName}_CM5.cif` will be created. These files contain the partial charges in the `DDEC` and `CM5` schemes, respectively.

For details of the DDEC method please see [J. Chem. Theory Comput. 2012, 8, 8, 2844–2867](https://pubs.acs.org/doi/10.1021/ct3002199). For details of the CM5 method please see [J. Chem. Theory Comput. 2012, 8, 2, 527–541](https://pubs.acs.org/doi/full/10.1021/ct200866d).

### Band gap calculation

The band gap is calculated directly by CP2K with the input created by the `charge_density.py` script. The only requirement is to add some unnoccupied states to the simulation. The band gap accuracy is extremally dependent of the level of theory used in the simulation. Currently only the `PBE` functional is implemented, although it is not the best choice for band gap calculations. *In the future the `PBE0` and `HSE06` functionals will be implemented.*

- The `DFTD3(BJ)` and `DFTD3` dispersion corrections are supported.
- The `DZVP` and `TZV2P` basis sets are supported.

```bash
#!/bin/bash

# Define global environment variables
export PATH=$PATH:/dccstor/nanopore-2945/electronic-density-experiment/bin/

export CP2K_DIR=/dccstor/nanopore-2945/cp2k/cp2k-v2023.1
export CP2K_DATA_DIR=${CP2K_DIR}/data

export OMP_NUM_THREADS=1

# Define specific variables
FrameworkName='MgMOF-74'
OutputFolder=$PWD
NProcs=36

# Load the CP2K setup script
source ${CP2K_DIR}/tools/toolchain/install/setup

# First ensure that the CIF file has all atoms explicitly defined in a P1 cell
echo -e "\nCreating primitive P1 cell..."
pmg structure --convert --filename ${OutputFolder}/${FrameworkName}.cif ${OutputFolder}/${FrameworkName}_prim.cif
mv -v ${OutputFolder}/${FrameworkName}_prim.cif ${OutputFolder}/${FrameworkName}.cif

echo -e "\nCreating CP2K input file..."
charge_density.py --FrameworkName ${FrameworkName} \
                  --Functional "PBE" \
                  --DispersionCorrection "DFTD3(BJ)" \
                  --PWCutoff 1200 \
                  --BasisSet "TZV2P" \
                  --UseScalapack  \
                  --NLUMO 10 \
                  ${OutputFolder}

echo -e "\nRunning CP2K simulation with ${NProcs} MPI and ${OMP_NUM_THREADS} OMP process..."
mpirun -np ${NProcs} $CP2K_DIR/exe/local/cp2k.psmp -i simulation_SCF.inp -o simulation_SCF.out

```

Here 10 unnocupied orbitals were added to the calculations in order to estimate the band gap. You can get the HOMO-LUMO gap from the `simulation_SCF.out` file searching for the line with `HOMO - LUMO gap [eV] :`.

### Molecular Orbital Plots

The plot of the molecular orbitals can be activate using the `--WriteMO` option. You also need to provide the number of orbitals that you want to be plotted. This can be done with the options `--NHOMO X` for ploting the `X` occupied orbitals below the Fermi level and `--NLUMO Y` for plotting the `Y` unnocupied orbitals above the Fermi level.

Currently only the `PBE` functional is implemented, although it is not the best choice for band gap calculations. *In the future the `PBE0` and `HSE06` functionals will be implemented.*

- The `DFTD3(BJ)` and `DFTD3` dispersion corrections are supported.
- The `DZVP` and `TZV2P` basis sets are supported.

```bash
#!/bin/bash

# Define global environment variables
export PATH=$PATH:/dccstor/nanopore-2945/electronic-density-experiment/bin/

export CP2K_DIR=/dccstor/nanopore-2945/cp2k/cp2k-v2023.1
export CP2K_DATA_DIR=${CP2K_DIR}/data

export OMP_NUM_THREADS=1

# Define specific variables
FrameworkName='MgMOF-74'
OutputFolder=$PWD
NProcs=36

# Load the CP2K setup script
source ${CP2K_DIR}/tools/toolchain/install/setup

# First ensure that the CIF file has all atoms explicitly defined in a P1 cell
echo -e "\nCreating primitive P1 cell..."
pmg structure --convert --filename ${OutputFolder}/${FrameworkName}.cif ${OutputFolder}/${FrameworkName}_prim.cif
mv -v ${OutputFolder}/${FrameworkName}_prim.cif ${OutputFolder}/${FrameworkName}.cif

echo -e "\nCreating CP2K input file..."
charge_density.py --FrameworkName ${FrameworkName} \
                  --Functional "PBE" \
                  --DispersionCorrection "DFTD3(BJ)" \
                  --PWCutoff 1200 \
                  --BasisSet "TZV2P" \
                  --UseScalapack  \
                  --WriteMO \
                  --NHOMO 10 \
                  --NLUMO 10 \
                  ${OutputFolder}

echo -e "\nRunning CP2K simulation with ${NProcs} MPI and ${OMP_NUM_THREADS} OMP process..."
mpirun -np ${NProcs} $CP2K_DIR/exe/local/cp2k.psmp -i simulation_SCF.inp -o simulation_SCF.out

```

Here the 10 highest energy occupied orbitals and 10 lowest unnocupied orbitals will be saved as `.cube` files. It is possible to use softwares like [VESTA](https://jp-minerals.org/vesta/en/) or [Avogadro](https://two.avogadro.cc/) to visualize the orbitals.

## Geometry optimization

The geometry optimization can be divided into two main steps:

- The `structure_optimization.py` script is used to generate the input for CP2K simulation. Here you can define the functional, the basis set, the dispersion correction, the plane wave cutoff, etc.
  - The `PBE` and `XTB` functionals are supported.
  - The `DFTD3(BJ)` and `DFTD3` dispersion corrections are supported.
  - The `DZVP` and `TZV2P` basis sets are supported.
  - The `--KeepSymmetry` can be used to force that the inicial symmetry of the structure is kept during the optimization.

- The `parse_optimization.py` script is used to parse the CP2K output file and create a new cif file with the optimized structure.

```bash
#!/bin/bash

# Define global environment variables
export PATH=$PATH:/dccstor/nanopore-2945/electronic-density-experiment/bin/

export CP2K_DIR=/dccstor/nanopore-2945/cp2k/cp2k-v2023.1
export CP2K_DATA_DIR=${CP2K_DIR}/data

export OMP_NUM_THREADS=1

# Define specific variables
FrameworkName='MgMOF-74'
OutputFolder=$PWD
NProcs=36

# Load the CP2K setup script
source ${CP2K_DIR}/tools/toolchain/install/setup

# First ensure that the CIF file has all atoms explicitly defined in a P1 cell
echo -e "\nCreating primitive P1 cell..."
pmg structure --convert --filename ${OutputFolder}/${FrameworkName}.cif ${OutputFolder}/${FrameworkName}_prim.cif
mv -v ${OutputFolder}/${FrameworkName}_prim.cif ${OutputFolder}/${FrameworkName}.cif

echo -e "\nCreating CP2K input file..."
structure_optimization.py --FrameworkName ${FrameworkName} \
                          --Functional "PBE" \
                          --DispersionCorrection "DFTD3(BJ)" \
                          --PWCutoff 1200 \
                          --KeepSymmetry \
                          --BasisSet "TZV2P" \
                          ${OutputFolder}

echo -e "\nRunning CP2K simulation with ${NProcs} MPI and ${OMP_NUM_THREADS} OMP process..."
mpirun -np ${NProcs} $CP2K_DIR/exe/local/cp2k.psmp -i simulation_Optimization.inp -o simulation_Optimization.out

parse_optimization.py --FrameworkName ${FrameworkName} \
                      --SaveHistory  \
                      ${OutputFolder}
```

The `--SaveHistory` option will save all the intermediary structures of the optimization as independet `.cif` files. The final optimized structure will be saved as `{FrameworkName}_optimized.cif`.

## Vibrational frequencies calculation

The vibrational frequencies are calculated numerically using the finite difference method. The `raman_ir.py` script is used to generate the input for CP2K simulation.

- The `PBE` and `XTB` functionals are supported.
- The `DFTD3(BJ)` and `DFTD3` dispersion corrections are supported.
- The `DZVP` and `TZV2P` basis sets are supported.
- The `--CalculateRaman` option can be used to calculate the Raman spectrum.
- The `--CalculateIR` option can be used to calculate the IR spectrum.
- The `--ProcsPerReplica` option can be used to define the number of MPI processes per replica. The total number of MPI processes will be `NProcs = NReplicas * ProcsPerReplica`.
- The `--dX` option can be used to define the step size used in the finite difference method. Default value is `0.001`.

The `parse_vibrations.py` script is used to parse the CP2K output file and create a new cif file with the vibrational frequencies and intensities.

- The `--SaveVibrations` option can be used to save the vibrations as `.axsf` files.
- The `--HalfWidth` option can be used to define the half width of the Lorentzian function used to create the Raman and/or IR curve.

```bash
#!/bin/bash

# Define global environment variables
export PATH=$PATH:/dccstor/nanopore-2945/electronic-density-experiment/bin/

export CP2K_DIR=/dccstor/nanopore-2945/cp2k/cp2k-v2023.1
export CP2K_DATA_DIR=${CP2K_DIR}/data

export OMP_NUM_THREADS=1

# Define specific variables
FrameworkName='MgMOF-74'
OutputFolder=$PWD
NProcs=36

# Load the CP2K setup script
source ${CP2K_DIR}/tools/toolchain/install/setup

# First ensure that the CIF file has all atoms explicitly defined in a P1 cell
echo -e "\nCreating primitive P1 cell..."
pmg structure --convert --filename ${OutputFolder}/${FrameworkName}.cif ${OutputFolder}/${FrameworkName}_prim.cif
mv -v ${OutputFolder}/${FrameworkName}_prim.cif ${OutputFolder}/${FrameworkName}.cif

echo -e "\nCreating CP2K input file..."
raman_ir.py --FrameworkName ${FrameworkName} \
            --Functional "PBE" \
            --DispersionCorrection "DFTD3(BJ)" \
            --PWCutoff 1200 \
            --BasisSet "TZV2P" \
            --ProcsPerReplica 9 \
            --dX 0.001 \
            --CalculateRaman \
            --CalculateIR \
            --UseScalapack \
            ${OutputFolder}

echo -e "\nRunning CP2K simulation with ${NProcs} MPI and ${OMP_NUM_THREADS} OMP process..."
mpirun -np ${NProcs} $CP2K_DIR/exe/local/cp2k.psmp -i simulation_Vibrations.inp -o simulation_Vibrations.out

echo -e "\nCreating RASPA P1 cell input file..."
parse_vibrations.py --FrameworkName ${FrameworkName} \
                    --SaveVibrations  \
                    --HalfWidth  5.0 \
                    ${OutputFolder}
```

The `parse_vibrations.py` creates two csv files with the results. The `{FrameworkName}_RAMAN_IR_Curve.csv` that holds the Raman and/or the IR curve created by adjusting a Lorentzian function to the calculated frequencies. The `{FrameworkName}_VibrationsTable.csv` that holds the calculated frequencies and intensities. The `--HalfWidth` option defines the half width of the Lorentzian function used to create the Raman and/or IR curve.

The `--SaveVibrations` option creates a folder called `Vibrations` with `AXSF` files containing the displacements of the atoms for each vibration mode. It is possible to visualize the vibrations using the [VESTA](https://jp-minerals.org/vesta/en/) software.

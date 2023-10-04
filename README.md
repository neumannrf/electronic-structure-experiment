# Electronic Structure Experiment

Automated virtual experiment that uses the CP2K software package to run DFT simulations of materials.

## Prerequisites

### General Prerequisites

- [CP2K](https://www.cp2k.org/download)
- [Chargemol](https://sourceforge.net/projects/ddec/files/)
- [Conda](https://www.anaconda.com/)

### Python Prerequisites

- jupyter
- jupyterlab
- numpy
- pandas
- pip
- python>=3.10
- ase
- gemmi
- pymatgen
- cp2k-input-tools

*Currently the Electronic Structure Experiment is only supported on Linux operating systems.*

## Installation

1. Install CP2K, Chargemol, and Conda.
2. Clone this repository.
3. Create a Conda environment using the `environment.yml` file in the root directory of this repository. The environment can be created using the following command: `conda env create -f environment.yml`.
4. Activate the Conda environment.

## Usage

Currently the Electronic Structure Experiment package can be used to perform the following tasks:

- Band gap calculation
- HOMO and LUMO orbital plots
- Partial charges calculation
- Geometry optimization
- Vibrational frequencies calculation and FTIR/Raman spectra simulation.

To a list of examples of how to use this repository, see the `docs` directory.

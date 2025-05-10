# FMOkit: A Python toolkit for the preparation and analysis of FMO calculations

## Introduction
**FMOkit** is a command-line toolkit written in Python for both preprocessing and postprocessing of FMO (Fragment Molecular Orbital) calculations.  
It consists of three main commands:

- **`mmcifprep`**  
  Adds hydrogen atoms and assigns partial charges to PDB structures. It also repairs missing residues or loops to produce mmCIF files suitable for FMO calculations.

- **`mmcif2fmoinp`**  
  Generates GAMESS input files for FMO calculations from the mmCIF files processed by `mmcifprep`.

- **`gamout2tsv`**  
  Processes GAMESS output files and extracts data into analysis-ready TSV format.

---

## Installation

To install FMOkit, clone the repository and install it in editable mode:

```bash
git clone https://github.com/kzfm/FMOkit.git
cd FMOkit
pip install -e .
```

This will register the commands mmcifprep, mmcif2fmoinp, and gamout2tsv for use in your terminal.

## Usage

### 1. Preprocess a structure for FMO calculations

```bash
mmcifprep INPUT_FILE LIGAND_ID SMILES_STRING
```

### 2. Generate GAMESS input files

```bash
mmcif2fmoinp INPUT
```

### 3. Extract results into a TSV file

```bash
gamout2tsv GAMOUT
```

## Example

## License

Code released under the [BSD license](LICENSE).

## Contributions

Pull requests and feature suggestions are welcome.
Please open an issue to discuss what you’d like to change.
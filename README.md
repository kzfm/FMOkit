# FMOkit: A Python toolkit for the preparation and analysis of FMO calculations

## Introduction
**FMOkit** is a command-line toolkit written in Python for both preprocessing and postprocessing of FMO (Fragment Molecular Orbital) calculations.  
It consists of three main commands:

- **`mmcifprep`**  
  Adds hydrogen atoms and assigns partial charges to PDB structures. It also repairs missing residues or loops to produce mmCIF files suitable for FMO calculations.

- **`mmcif2fmoinp`**  
  Generates GAMESS input files for FMO calculations from the mmCIF files processed by `mmcifprep`.

- **`gamoutparser`**  
  Processes GAMESS output files and extracts data into analysis-ready CSV/TSV format.

## Installation

To install FMOkit, clone the repository and install it in editable mode:

```bash
git clone https://github.com/kzfm/FMOkit.git
cd FMOkit
pip install -e .
```

This will register the commands mmcifprep, mmcif2fmoinp, and gamoutparser for use in your terminal.

## Usage

### 1. Preprocess a structure for FMO calculations

```bash
mmcifprep INPUT_FILE LIGAND_ID SMILES_STRING
```

As this is the most important step, we recommend referring to [the accompanying Jupyter notebook](https://github.com/kzfm/FMOkit/blob/main/examples/prepare_complex-5law.ipynb), which demonstrates the process interactively and step by step.

> **Note:**
> If you are using molecular modeling tools such as Maestro or MOE, you can skip this step.
> Instead of using mmcifprep, please perform preprocessing steps such as hydrogen addition within those tools and export the structure as an  mmCIF format file.

### 2. Generate GAMESS input files

```bash
mmcif2fmoinp INPUT
```

Generates a GAMESS input file for FMO calculations from a preprocessed mmCIF file.
Optionally, you can specify GAMESS computational resources such as the number of nodes, CPU cores, and memory size.
The system is automatically divided into fragments based on residue IDs.
Each fragment’s formal charge is computed by summing the partial charges of the atoms it contains, so partial charge information is required.

> **Note**  
> If you are converting a CIF file exported from **Maestro**, there are two important points to keep in mind:  
> 1. Atomic charge information is stored in the `pdbx_formal_charge` attribute.  
> 2. The chain ID corresponds to `auth_asym_id` instead of `label_asym_id`.
> 
> Therefore, you need to specify the following two options explicitly:
> 
> ```bash
> mmcif2fmoinp INPUT -–asym_id=auth_asym_id -–charge=pdbx_formal_charge
> ```
> 
> Alternatively, you can use the predefined shortcut option:
> 
> ```bash
> mmcif2fmoinp -M INPUT
> ```


### 3. Extract results into a CSV/TSV file

```bash
gamoutparser GAMOUT
```

Parses GAMESS output files and extracts PIEDA analysis results into CSV/TSV format, making the data compatible with common data science tools.

## Example

## License

Code released under the [BSD license](LICENSE).

## Contributions

Pull requests and feature suggestions are welcome.
Please open an issue to discuss what you’d like to change.

## Version history

- 0.1: initial version
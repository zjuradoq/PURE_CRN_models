# A chemical reaction network model of PURE
[![bioRxiv](https://img.shields.io/badge/PDF-bioRxiv-red)](https://www.biorxiv.org/content/10.1101/2023.08.14.553301v1.full.pdf)
[![License](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/license/mit/)

Data and models repository for the paper "A chemical reaction network model of PURE," along with additional simulations and validation experiments. The folder contains all Python scripts used to simulate protein expression and plot figures seen in the paper "A chemical reaction network model of PURE" by Zoila Jurado, Ayush Pandey, and Richard M. Murray. Additionally, all simulation runs and experimental results are included.

# Installation

1. Run `pip install -r requirements.txt` from your terminal.

Or, if you prefer to install everything manually, you can:
1. Install biocrnpyler: `pip install biocrnpyler>=1.1.1`
2. Install bioscrape: `pip install bioscrape>=1.2.0`
3. If Step 2 fails, you can clone the [Bioscrape](https://github.com/biocircuits/bioscrape/) repository and manually install the package. Run `python setup.py install` to install bioscrape once you meet all the [installation requirements](https://github.com/biocircuits/bioscrape/wiki/Installation).
4. Install seaborn, holoviz, corner, svglib, scikit-learn, bokeh-catplot, openpyxl, bokeh, jupyter_bokeh ==2.0.0, and panel==0.13.1

# How to navigate this repository

## PURE CRN models
We use BioCRNpyler to generate chemical reaction network models to building a transcription (TX) and translation (TL) model for cell-free protein synthesis using PURE. The generation of all models are in the `modeling/` directory in this repository.

## Inferencing on TX-only data 
Scripts and results for the parameter inferencing for the TX-only model on MGapt expression is under the `Inferencing_MGapt` directory.

## Experimental data
All calibration data used of this repository is available under the `Data_files` directory. Additional calibration data use can be found under the 'Calibration_MGapt','Calibration_MGapt' directory. The plasmids used in this paper are available from [myTXTL T7 Expression Kit](https://arborbiosci.com/products/cell-free-protein-synthesis/mytxtl-cell-free-expression-kits/mytxtl-t7-expression-kit/).

## Simulation Experiments
All simulations from the model are available under the `Simulation_results` directory in folders corresponding with model run. Subfolder `Nuclues_Simulations` contains, specific simulations around the separation of transcription and translation under different nucleus conditions.

## SBML models
All generated SBML models are available in the `Models` directory.

## Scripts
All scripts used to build model and plotting are located in the main folder. 

## References
The PURE CRN model is based on [Matsuura 2017](https://www.pnas.org/doi/full/10.1073/pnas.1615351114). For more details, see [PURE simulatator](https://sites.google.com/view/puresimulator).

# Contact

Contact Zoila Jurado if you have any questions.

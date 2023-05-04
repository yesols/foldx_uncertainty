# foldx_uncertainty

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7897628.svg)](https://doi.org/10.5281/zenodo.7897628)

## Overview

The code and datasets here accompany the manuscript, "Statistical Modeling to Quantify the Uncertainty of FoldX-Predicted Protein Folding and Binding Stability" by Sapozhnikov et al.

All original datasets as well as intermediate outputs are already present in appropriate folders, and will be simply replaced if you re-run the codes. Some of the model search, especially the best subset selection, may take awhile. Code files are numbered according to the order of analysis.

For applying our model to your own data, skip to the file 07_application.Rmd.


## Description of `data/`

`dssp/`: output and processed files from running DSSP

`fold_xstal_data/`: ddG of folding, from experimental structures

`bind_xtal_data/`: ddG of binding, from experimental structures

`binding_exp/`: corrected experimental ddG of binding for MD+FoldX dataset (the experimentally measured ddG for all others are contained in their respective folders, together with FoldX ddG)

`fold_data/`: ddG of folding, from MD snapshots (MD+FoldX workflow)

`Bind_data/`: ddG of binding, from MD snapshots (MD+FoldX workflow)

## Description of `outputs/`

`best_mods_all.rds`: Final models 1-5 for folding and binding. A list containing ten `lm` objects.

`models.rds`: Models 1 and 5 produced from stepwise selection and from best subset selection, for folding and binding. A list containing eight `lm` objects. Better models from the two methods are contained in `best_mods_all.rds`.

`csv` files: Dataframes containing all mutations and variables used in the model. Separate dataframes were created for model training by removing columns that are not variables.

## Description of `scripts/`

`01_generate-organized-tables.Rmd`: Reads in the data in `data/` folder for experimental and FoldX ddG, and adds biochemical property info. Outputs `csv` files to `outputs/` folder.

`02_model-selection.R`: Reads in the `csv` files in `outputs/` folder and performs model search for Models 1 (simplest) and 5 (full). Outputs models as a list of `lm` objects to `outputs/models.rds`.

`03_cross-validation.Rmd`: Reads in the `csv` files and `models.rds` in `outputs/` folder and performs leave-one-system-out cross validation to pick the better model between stepwise selection and best subset selection. 

`04_other-mods.Rmd`: Reads in the `csv` files and `models.rds` in `outputs/` folder and performs model search and cross validation for Models 2-4. Outputs final models as a list of `lm` objects to `outputs/best_mods_all.rds`.

`05_figures.Rmd`: Reads in the `csv` files and `models.rds` in `outputs/` folder and generates figures that are in the manuscript. Outputs to `figures/` folder.

`06_summary.Rmd`: Reads in the `csv` files and `best_mods_all.rds` in `outputs/` folder and prints model details for all ten models.

`07_application.Rmd`: Reads in `outputs/best_mods_all.rds` to make the uncertainty prediction on user-provided data. No output file.



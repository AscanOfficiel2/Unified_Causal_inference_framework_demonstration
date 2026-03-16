# Unified Causal Inference Framework – CRC Microbiome Demonstration
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19047082.svg)](https://doi.org/10.5281/zenodo.19047082)

This repository provides the code used to demonstrate the unified causal inference framework in our review paper (Briefings in Bioinformatics, 2026). Using eight public CRC metagenomic datasets available from NCBI SRA Read Archive (PRJDB4176, PRJEB27928, PRJEB6070, PRJEB72525, PRJEB72526, PRJEB7774, PRJNA1167935, PRJNA389927), we apply:

- Data harmonization and CLR transformation
- DML estimation with cross-fitting (EconML)
- Bootstrap resampling (300 iterations)
- Negative control validation (*Acidaminococcus intestini*)
- E-value sensitivity analysis

Results: The detailed results can be found in the DML_bootstrap_negative_control_results.xlsx. 

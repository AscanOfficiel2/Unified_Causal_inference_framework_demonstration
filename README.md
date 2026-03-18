# Unified Causal Inference Framework – CRC Microbiome Demonstration
This repository provides the code used to demonstrate the unified causal inference framework for our manuscript currently under review at Briefings in Bioinformatics (2026). It includes detailed results and figures generated from eight public colorectal cancer (CRC) metagenomic datasets available from the NCBI Sequence Read Archive (SRA) under accessions: PRJDB4176, PRJEB27928, PRJEB6070, PRJEB72525, PRJEB72526, PRJEB7774, PRJNA1167935, and PRJNA389927. The analysis applies:

- Data harmonization and CLR transformation
- DML estimation with cross-fitting (EconML)
- Bootstrap resampling (300 iterations)
- Negative control validation (*Acidaminococcus intestini*)
- E-value sensitivity analysis

## Outputs: 
- The code for the DAG construction in dagitty
- The detailed statistical results can be found in the DML_bootstrap_negative_control_results.xlsx. 
- The figure of the DAG (saved as model.png)
- The forest plot with the negative control
- Bootstrap distributions

## Citation
When applying the code or using the results in this repository, please cite:

Ascandari A, Aminu S, Benhida R, Daoud R. (2026). From Association to Causation: A Decision-Aware Framework for Reproducible Biomarker Discovery and Precision Intervention Design in the Human Gut Microbiome. Briefings in Bioinformatics.

This repository is permanently archived at Zenodo with DOI : 
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19047082.svg)](https://doi.org/10.5281/zenodo.19047082)

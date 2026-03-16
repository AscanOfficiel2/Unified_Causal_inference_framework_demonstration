# Unified Causal Inference Framework – CRC Microbiome Demonstration
This repository provides the code used to demonstrate the unified causal inference framework in our review paper (Briefings in Bioinformatics, 2026). Using eight public CRC metagenomic datasets available from NCBI SRA Read Archive (PRJDB4176, PRJEB27928, PRJEB6070, PRJEB72525, PRJEB72526, PRJEB7774, PRJNA1167935, PRJNA389927), we apply:

- Data harmonization and CLR transformation
- DML estimation with cross-fitting (EconML)
- Bootstrap resampling (300 iterations)
- Negative control validation (Acidaminococcus intestini)
- E-value sensitivity analysis

Results: F. nucleatum shows a positive association with CRC (β = 0.026, 95% CI: 0.021–0.032); the negative control shows no effect. All code is provided for reproducibility

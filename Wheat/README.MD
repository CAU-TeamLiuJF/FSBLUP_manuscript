# Benchmarking analyses using data from Krause et al 2019

## Data download and processing
The following files were downloaded from: `https://data.cimmyt.org/dataset.xhtml?persistentId=hdl:11529/10548109`

- `Krause_et_al_2018_Yield_iid_BLUPs.csv`
- `Krause_et_al_2018_Yield_BLUEs.csv`
- `Krause_et_al_2018_Hyper_BLUEs_Individual_Time_Points.csv`
- `Krause_et_al_2018_Genotypes.csv`

Calculation of the kinship matrices was done with the script: `data_prep_Krause.R` (note the commented lines must be run one time to generate the K matrix)

## Analyses

The following scripts run the analyses:
- `univariate_GP.R`
- `MegaLMM_Krause.R`
- `GOBULP.R`

These can be parallelized on a cluster using the `*.sh` files.

Note that the file path should change to your workpath.

Results can be collected with `collect_results.R`, and then figures made with `Figures.R`.

The manuscript results are save in `results` directory.


[note] The scripts `MegaLMM_Krause.R` and `Estimate_gcor_prediction.R` were download from `https://github.com/deruncie/MegaLMM_analyses/blob/v0.9.2` , which supplied in article *Runcie DE, Qu J, Cheng H, Crawford L. MegaLMM: Mega-scale linear mixed models for genomic predictions with thousands of traits. Genome Biol. 2021 Jul 23;22(1):213. doi: 10.1186/s13059-021-02416-w.*
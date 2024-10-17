# Benchmarking analyses using pig data 

raw data are available from the corresponding authors upon reasonable request.
 
All matrices and phenotypes that replicate the results in the manuscript are stored in `Pig\data\data.rdata`.

## Data processing

Calculation of the kinship matrices was done with the script: `data_prep.R` (note the commented lines must be run one time to generate the K matrix).

## Analyses

The following scripts run the analyses:
-`univariate.R`
-`GOBLUP.R`

These can be parallelized on a cluster using the `*.sh` files.

Note that the file path should change to your workpath.

The manuscript results are save in `results` directory.
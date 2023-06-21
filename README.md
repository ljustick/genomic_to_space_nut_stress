# Genomic-to-space measurements reveal global ocean nutrient stress

Code to accompany the manuscript: Genomic-to-space measurements reveal global ocean nutrient stress, and description of data files.

## Data files
The following data files can be downloaded here: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8064616.svg)](https://doi.org/10.5281/zenodo.8064616)


in_situ_data.csv
- csv file containing *in situ* measurements. Collection of metadata from Ustick *et al.* 2021 https://doi.org/10.1126/science.abe6301.

sat_matchup_data.RData
- Theta prime observations averaged across different time and spatial scales matched to the time and sampling location of metagenomic samples. Stored as RData file.

thPrime_1degree_032223.mat
- Global theta prime observations, averaged in 1 degree across 8 day windows.

modisSST_1degree8Day.mat
- Sea surface temperature, averaged in 1 degree across 8 day windows.


## Code
in_situ_analysis.R
- R script for selecting best spatial and temporal matchup between remote sensing and in situ data, and for Random Forest Analysis.
- Written by: Lucas J. Ustick


temporal_analysis_scripts
- Matlab scripts used for temporal of theta prime.
- Written by: Adam C. Martiny
# Genomic-to-space measurements reveal global ocean nutrient stress
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8065231.svg)](https://doi.org/10.5281/zenodo.8065231)

Code to accompany the manuscript: Genomic-to-space measurements reveal global ocean nutrient stress, and description of data files.

## Data files
The following data files can be downloaded here: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8064615.svg)](https://doi.org/10.5281/zenodo.8064615)

data_insitu_v2_25-04.csv
- csv file containing *in situ* measurements. Collection of metadata from Ustick *et al.* 2021 https://doi.org/10.1126/science.abe6301.

sat_data_v2_25-04
- Theta prime observations averaged across different time and spatial scales matched to the time and sampling location of metagenomic samples. Saved as a merged xlsx file and separate csv files .

theta_v2_25-04.csv
- Global theta prime observations, averaged in 1 degree across 8 day windows. Saved as a csv file with the first two rows and columns representing the Year, Day, Degree Latitude, and Degree Longitude, with the cells containing the corresponding values across these four variables.

SST_v2_25-04.csv
- Sea surface temperature, averaged in 1 degree across 8 day windows. Saved as a csv file with the first two rows and columns representing the Year, Day, Degree Latitude, and Degree Longitude, with the cells containing the corresponding values across these four variables.

## Code
calc_thetaprime_v2.m
- Matalab script for calculating theta prime from remote sensing products.
- Written by: Toby K. Westberry

random_forest_analysis_v2.R
- R script for selecting best spatial and temporal matchup between remote sensing and in situ data, and for Random Forest Analysis.
- Written by: Lucas J. Ustick

temporal_analysis_scripts_v2
- Matlab scripts used for temporal analysis of theta prime.
- Written by: Adam C. Martiny
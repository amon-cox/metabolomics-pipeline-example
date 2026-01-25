## Script for performing the processing, analysis, and visualization.
# Calls all scripts in R/
# This script must inhabit and be run from the main folder, metabolomics-pipeline-example/
source("R/00_setup.r")
source("R/01_preprocess.r")
source("R/02_sample_analysis.r")
source("R/03_feature_analysis.r")

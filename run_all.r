## Script for performing the processing, analysis, and visualization.
# Calls all scripts in R/
# This script must inhabit and be run from the main folder, metabolomics-pipeline-example/
source(file.path("R", "00_setup.r")) # file.path() used to ensure funtionality across macOS, Linux, Windows
source(file.path("R", "01_preprocess.r"))
source(file.path("R", "02_sample_analysis.r"))
source(file.path("R", "03_feature_analysis.r"))
source(file.path("R", "04_supervised.r"))
source(file.path("R", "05_significant_features.r"))

if (!file.exists(file.path("data_processed", "06_kegg_compounds.tsv"))) { # checks if KEGG info is already present
  source(file.path("R", "06a_download_kegg_cpds.R")) # time-consuming script to access full compound info from KEGG
}

#source(file.path("R", "06b_annotate_with_kegg.r"))

# execute the Quarto overview file
quarto::quarto_render(input = "overview.qmd", output_format = "pdf")

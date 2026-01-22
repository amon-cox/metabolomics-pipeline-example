## 01) Organizing both metabolomics datasets and establishing groups.
# load raw LC-MS/MS peaks data
peaks_negative <- read.csv("data_raw/peaks_negative.csv", header = TRUE)
peaks_positive <- read.csv("data_raw/peaks_positive.csv", header = TRUE)

## log2 transformation and standardization, if needed (may hold on standardization)
## prepare feature_metadata for both sets
## export cleaned features with mz_rt_min IDs, no feature metadata
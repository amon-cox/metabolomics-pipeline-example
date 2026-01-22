## 02) Overview analysis of sample behavior
# load sample metadata
sample_metadata <- read.csv("data_raw/sample_metadata.csv", header = TRUE) |>
    mutate(treatment = factor(treatment, levels = c("control", "tr1", "tr2", "tr1+tr2")))

## remember, the peaks_X_log2 objects are already in the environment if running these scripts in sequence.
set.seed(123)

## reorient data structure, apply PCA, extract scores and plot, extract loadings and plot
# 2026-01-22
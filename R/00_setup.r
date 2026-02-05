## 00) Set up R environment
# load necessary packages
library(vegan)
library(pairwiseAdonis)
library(limma)
library(mixOmics)
library(ggplot2)
library(viridis)
library(cowplot)
library(rmarkdown)
library(quarto)
library(tidyr)
library(tibble)
library(dplyr) # loaded last to give dplyr priority in function name conflicts

# load sample metadata
sample_metadata <- read.csv("data_raw/sample_metadata.csv", header = TRUE) |>
    mutate(treatment = factor(treatment, levels = c("control", "tr1", "tr2", "tr1+tr2")))

# establish theme for most plots
ggplot2::theme_set(
    cowplot::theme_half_open(12) +
        ggplot2::theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.margin = margin(0.2, 0.2, 0.2, 0.2, unit = "cm")
        )
)
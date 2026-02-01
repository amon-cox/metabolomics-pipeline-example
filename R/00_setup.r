## 00) Set up R environment
# load necessary packages
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(vegan)
library(pairwiseAdonis)
library(limma)
library(viridis)
library(cowplot)
library(rmarkdown)
library(quarto)

# establish theme for most plots
ggplot2::theme_set(
    cowplot::theme_half_open(12) +
        ggplot2::theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.margin = margin(0.2, 0.2, 0.2, 0.2, unit = "cm")
        )
)
## 00) Set up R environment
# load necessary packages
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(viridis)
library(cowplot)
library(limma)

# establish theme for plots
ggplot2::theme_set(
    cowplot::theme_half_open(12) +
        ggplot2::theme(
            panel.grid.major = element_line(color = "grey80"),
            panel.grid.minor = element_blank(),
            plot.margin = margin(6, 0, 6, 0)
        )
)
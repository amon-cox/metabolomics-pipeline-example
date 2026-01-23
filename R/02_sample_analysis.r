## 02) Overview analysis of sample behavior
if(!interactive()) pdf(NULL) # prevents pdf artifact from being generated from plots

# Principal Components Analysis
set.seed(123)

do_pca <- function(.peaks, .mode) {

    # conduct PCA and prepare its outputs
    pca <- .peaks |>
        tibble::column_to_rownames(var = "mz_rt_min") |> # convert column to rownames before transpose
        t() |> # rotate matrix
        as.data.frame() |>
        prcomp(center = TRUE, scale. = TRUE) # apply PCA after standardizing data
    
    scores_matrix <- pca$x |> # extract principal components
        as.data.frame() |>
        tibble::rownames_to_column(var = "sample")
    
    loadings_matrix <- pca$rotation |> # extract feature contributions from PC space
        as.data.frame() |>
        tibble::rownames_to_column(var = "peak")
    
    variance <- round(pca$sdev^2 / sum(pca$sdev^2), 4) * 100 # convert proportions of variance into percentage

    # export data objects
    write.csv(scores_matrix, file = paste0("data_processed/pca_scores_", .mode, "_mode.csv"), row.names = FALSE, quote = FALSE)
    write.csv(loadings_matrix, file = paste0("data_processed/pca_loadings_", .mode, "_mode.csv"), row.names = FALSE, quote = FALSE)

    # prepare scores and loadings as a list output
    pca_object <- list(
        scores = scores_matrix,
        loadings = loadings_matrix,
        prop_var = variance
    )

}

pca_negative <- do_pca(.peaks = peaks_negative_log2, .mode = "negative")
pca_positive <- do_pca(peaks_positive_log2, "positive")

# load sample metadata
sample_metadata <- read.csv("data_raw/sample_metadata.csv", header = TRUE) |>
    mutate(treatment = factor(treatment, levels = c("control", "tr1", "tr2", "tr1+tr2")))

# plotting the principal components
plot_scores <- function(.pca, .mode, .md = sample_metadata) {

    # set up 
    scores_plot <- .pca[['scores']] |>
        select(sample, PC1, PC2) |>
        left_join(y = .md) |>
        ggplot(aes(x = PC1, y = PC2, color = treatment)) +
            geom_point(aes(shape = cohort), alpha = 0.75, size = 2) +
            coord_fixed(ratio = 1) + # set ratio so that PC scales are apparent
            scale_color_viridis(discrete = TRUE, option = "viridis", direction = -1, drop = FALSE) +
            labs(
                x = paste0("PC1 (", .pca[['prop_var']][1], "%)"),
                y = paste0("PC2 (", .pca[['prop_var']][2], "%)"),
                title = paste0(.mode, "-mode"),
                subtitle = "LC-MS/MS peaks"
            )

}

scores_plot_negative <- plot_scores(.pca = pca_negative, .mode = "negative")
scores_plot_positive <- plot_scores(.pca = pca_positive, .mode = "positive")

# arrange plots for export
plot_row <- cowplot::plot_grid(
    scores_plot_negative + theme(legend.position = "none"),
    scores_plot_positive + theme(legend.position = "none"),
    align = "vh",
    labels = "AUTO",
    hjust = -1,
    nrow = 1
)

legend <- cowplot::get_legend(
    scores_plot_negative + theme(legend.box.margin = margin(0, 0, 0, 12))
)

final_scores_plot <- cowplot::plot_grid(
    plot_row,
    legend,
    rel_widths = c(2, .4)
)

cowplot::save_plot(
    filename = "output/pca_scores_plot.png",
    plot = final_scores_plot,
    bg = "white"
)
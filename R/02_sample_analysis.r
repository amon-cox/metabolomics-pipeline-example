## 02) Overview analysis of sample behavior
if(!interactive()) pdf(NULL) # prevents pdf artifact from being generated from plots

# Principal Components Analysis
set.seed(123)

do_pca <- function(.intensities, .mode) {

    # conduct PCA and prepare its outputs
    pca <- .intensities |>
        tibble::column_to_rownames(var = "mz_rt_min") |> # convert column to rownames before transpose
        t() |> # rotate matrix
        as.data.frame() |>
        prcomp(center = TRUE, scale. = TRUE) # apply PCA after standardizing data
    
    scores_matrix <- pca$x |> # extract principal components
        as.data.frame() |>
        tibble::rownames_to_column(var = "sample")
    
    loadings_matrix <- pca$rotation |> # extract feature contributions from PC space
        as.data.frame() |>
        tibble::rownames_to_column(var = "mz_rt_min")
    
    variance <- round(pca$sdev^2 / sum(pca$sdev^2), 4) * 100 # convert proportions of variance into percentage

    # export data objects
    write.csv(scores_matrix, file = file.path("data_processed", paste0("02_pca_scores_", .mode, ".csv")), row.names = FALSE, quote = FALSE)
    write.csv(loadings_matrix, file = file.path("data_processed", paste0("02_pca_loadings_", .mode, ".csv")), row.names = FALSE, quote = FALSE)
    write.csv(variance, file = file.path("data_processed", paste0("02_pca_variance_", .mode, ".csv")), row.names = FALSE, quote = FALSE)

    # prepare scores and loadings as a list output
    pca_object <- list(
        scores = scores_matrix,
        loadings = loadings_matrix,
        prop_var = variance
    )

    return(pca_object)

}

pca_negative <- do_pca(.intensities = intensity_norm_log2_neg, .mode = "neg")
pca_positive <- do_pca(intensity_norm_log2_pos, "pos")

# detect outliers in PCs 1 & 2 space per treatment group
detect_outliers <- function(.pca, .mode, .md = sample_metadata, .group = "treatment", k_iqr = 1.5) {

    scores <- .pca[['scores']] |>
        select(sample, PC1, PC2) |> # grab sample IDs and the first two principal components
        left_join(y = .md, by = "sample")
    
    centroids <- scores |>
        group_by(.data[[.group]]) |> # set grouping to treatment variable
        summarize( # calculate centroid of each PC per treatment group
            centroid_PC1 = mean(PC1),
            centroid_PC2 = mean(PC2),
            .groups = "drop"
        )
    
    distances <- left_join( # copies the centroid info to each sample by group
        x = scores,
        y = centroids,
        by = .group
    ) |>
        mutate(group_dist = sqrt((PC1 - centroid_PC1)^2 + (PC2 - centroid_PC2)^2)) # Pythagorean theorem; samples' distance from centroid
    
    outliers <- distances |>
        group_by(.data[[.group]]) |> # set grouping to treatment variable
        mutate(
            Q1 = quantile(group_dist, 0.25), # calculate first quartile
            Q3 = quantile(group_dist, 0.75), # calculate third quartile
            IQR = Q3 - Q1, # calculate interquartile range
            upper_bound = Q3 + k_iqr * IQR, # establish maximum distance from group centroid
            is_outlier = group_dist > upper_bound # establish outlier criteria
        ) |>
        ungroup() |>
        mutate(mode = .mode) |>
        select(sample, is_outlier, mode, all_of(.group), group_dist, Q1, Q3, IQR, upper_bound, everything())
    
    return(outliers)

}

pca_outliers <- bind_rows( # merge function output for both modes
    detect_outliers(.pca = pca_negative, .mode = "neg"),
    detect_outliers(pca_positive, "pos")
)

write.csv(pca_outliers, file = file.path("output", "tables", "02_pca_outlier_flags.csv"), row.names = FALSE)

# plotting the principal components
plot_scores <- function(.pca, .mode, .md = sample_metadata) {

    # set up 
    p_scores <- .pca[['scores']] |>
        select(sample, PC1, PC2) |>
        left_join(y = .md, by = "sample") |>
        ggplot(aes(x = PC1, y = PC2, color = treatment)) +
            geom_point(aes(shape = cohort), alpha = 0.75, size = 2) +
            coord_fixed() + # set ratio so that axes are on the same scale
            scale_color_viridis(discrete = TRUE, option = "viridis", direction = -1, drop = FALSE) +
            theme(plot.subtitle = element_text(face = "bold")) +
            labs(
                x = paste0("PC1 (", .pca[['prop_var']][1], "%)"),
                y = paste0("PC2 (", .pca[['prop_var']][2], "%)"),
                subtitle = paste0("PCA on ", .mode, "-mode intensities")
            ) +
            geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
            geom_vline(xintercept = 0, linetype = "dashed", color = "black")

}

p_scores_negative <- plot_scores(.pca = pca_negative, .mode = "neg")
p_scores_positive <- plot_scores(.pca = pca_positive, .mode = "pos")

# arrange plots for export
p_scores_row <- cowplot::plot_grid(
    p_scores_negative + theme(legend.position = "none"),
    p_scores_positive + theme(legend.position = "none"),
    align = "vh",
    nrow = 1
)

legend <- cowplot::get_legend(
    p_scores_negative + theme(legend.box.margin = margin(0, 0, 0, 12))
)

p_scores_final <- cowplot::plot_grid(
    p_scores_row,
    legend,
    rel_widths = c(2, .4)
)

cowplot::save_plot(
    filename = file.path("output", "plots", "02_pca_scores_neg_pos.png"),
    plot = p_scores_final,
    bg = "white",
    base_height = 4
)

# Permutational Analysis of Variance on treatment groups. Hn is that all group centroids are the same
do_permanova <- function(.intensities, .mode, .md = sample_metadata) {

    distances <- .intensities |>
        tibble::column_to_rownames(var = "mz_rt_min") |> # convert column to rownames before transpose
        t() |> # rotate matrix
        as.data.frame() |>
        dist(method = "euclidean") # compute euclidean distances between samples

    permanova <- vegan::adonis2(
        formula = distances ~ treatment, # set up model for PERMANOVA
        data = .md # specify where the info for "treatment" comes from
    )

    pairwise_permanova <- pairwiseAdonis::pairwise.adonis( # conduct several paired analyses ...
        x = distances,
        factors = .md$treatment,
        p.adjust.m = "bonferroni" # ... then correct p values for multiple comparisons
    )

    write.csv( # export global PERMANOVA
        as.data.frame(permanova),
        file = file.path("output", "tables", paste0("02_permanova_global_", .mode, ".csv")),
        row.names = TRUE
    )

    write.csv( # export pairwise comparisons
        as.data.frame(pairwise_permanova),
        file = file.path("output", "tables", paste0("02_permanova_pairwise_", .mode, ".csv")),
        row.names = FALSE
    )

    group_stats <- list(
        distance = distances,
        global = permanova,
        pairwise = pairwise_permanova,
        mode = .mode
    )

    return(group_stats)

}

permanova_neg <- do_permanova(.intensities = intensity_norm_log2_neg, .mode = "neg")
permanova_pos <- do_permanova(intensity_norm_log2_pos, "pos")

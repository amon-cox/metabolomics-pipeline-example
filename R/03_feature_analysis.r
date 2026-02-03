## 03) Analysis of feature behavior
# analyze features using linear models for microarray data
design_matrix <- model.matrix(~ treatment, data = sample_metadata) # expands metadata into dummy variables for a design matrix
colnames(design_matrix) <- c("Intercept", "tr1", "tr2", "tr1_tr2") # replace "+" with "_" to satisfy R naming conventions

design_contrasts <- makeContrasts( # establish the contrast matrix of the treatment groups
    tr1_vs_ctrl = tr1,
    tr2_vs_ctrl = tr2,
    tr2_tr1_vs_ctrl = tr1_tr2,
    tr2_vs_tr1 = tr2 - tr1,
    levels = design_matrix
)

apply_limma <- function(.intensities, .mode, .md = sample_metadata, .design = design_matrix, .contrasts = design_contrasts) {

    correlations <- duplicateCorrelation( # estimate intra-block correlations
        object = .intensities,
        design = .design,
        block = .md$cohort
    )

    model <- lmFit( # fit linear model for each LC-MS/MS peak in the series
        object = .intensities,
        design = .design,
        correlation = correlations$consensus.correlation
    )

    limma_output <- contrasts.fit( # computes estimated coefficients and standard errors for the linear model contrasts
        fit = model,
        contrasts = .contrasts
    ) |>
        eBayes() # borrows information from across features
    
    return(limma_output)

}

limma_output_neg <- apply_limma(.intensities = column_to_rownames(intensity_norm_log2_neg, "mz_rt_min"), .mode = "neg")
limma_output_pos <- apply_limma(column_to_rownames(intensity_norm_log2_pos, "mz_rt_min"), "pos")

limma_res_neg <- topTable( # extracts features from a linear model fit for viewing in a table
    fit = limma_output_neg,
    adjust = "BH", # adjusts p values for multiple comparisons
    number = Inf # grabs all features, not just the top X features
) |>
    rownames_to_column("mz_rt_min")

limma_res_pos <- topTable(
    fit = limma_output_pos,
    adjust = "BH",
    number = Inf
) |>
    rownames_to_column("mz_rt_min")

write.csv(limma_res_neg, file = "output/tables/03_limma_global_neg.csv", row.names = FALSE)
write.csv(limma_res_pos, file = "output/tables/03_limma_global_pos.csv", row.names = FALSE)

export_limma_res <- function(.limma_output, .mode) { # function to export each contrast

    res_list <- list()

    for (coef_name in colnames(.limma_output$coefficients)) { # prepare to run through each contrast

        res <- topTable(
            fit = .limma_output, # take original limma output
            adjust = "BH", # adjust p values for multiple comparisons
            coef = coef_name, # specify which column/contrast
            number = Inf, # grab all features, for now
            sort.by = "p" # sort by p value
        ) |>
            rownames_to_column("mz_rt_min") |>
            rename( # simplifying some names
                log2FC = logFC,
                p = P.Value,
                p_adj_BH = adj.P.Val,
                log_odds = B
            )
        
        write.csv(
            filter(res, p_adj_BH < .05 & abs(log2FC) > 1), # export only significant results
            file = paste0("output/tables/03_limma_", coef_name, "_", .mode, ".csv"),
            row.names = FALSE
        )

        res_list[[coef_name]] <- res

    }

    # build volcano plots in a panel, per mode
    p_volcano <- res_list |>
        bind_rows(.id = "contrast") |> # collapse list into large dataframe
        mutate(sig_feature = p_adj_BH < .05 & abs(log2FC) > 1) |>
        ggplot(aes(x = log2FC, y = -log10(p_adj_BH), color = sig_feature)) +
            geom_point(alpha = 0.75) +
            scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
            facet_wrap(~ contrast, ncol = 2, nrow = 2) +
            scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "#48bda6")) +
            theme(panel.grid.major = element_blank()) +
            geom_vline(xintercept = 1, linetype = "dotted", color = "grey60") +
            geom_vline(xintercept = -1, linetype = "dotted", color = "grey60") +
            geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "grey80") +
            labs(
                x = "log2 fold change",
                y = "-log10 adj. p value",
                title = paste0("limma contrasts ", .mode, "-mode"),
                color = "significant"
            )

    cowplot::save_plot(
        filename = paste0("output/plots/03_limma_contrasts_", .mode, ".png"),
        plot = p_volcano,
        bg = "white"
    )

    invisible(res_list)
}

export_limma_res(.limma_output = limma_output_neg, .mode = "neg")
export_limma_res(limma_output_pos, "pos")

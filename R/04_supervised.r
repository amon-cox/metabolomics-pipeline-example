## 04) Supervised method for sample analysis
if(!interactive()) pdf(NULL) # prevents pdf artifact from being generated from plots

# conduct a sparse Partial Least Squares Discriminant Analysis
do_splsda <- function(.intensities, .mode, .md = sample_metadata) {

    # reformat data objects
    X_intensities <- .intensities |>
        column_to_rownames(var = "mz_rt_min") |>
        t() |>
        as.data.frame()
    
    Y_treatment <- .md$treatment
    study_cohort <- .md$cohort

    # establish sparsity parameters and tune for sparse variables
    test_keepX <- c(10, 25, 50, 100, 200, 400, 800)

    set.seed(123)

    tune_res <- tune.splsda(
        X = X_intensities,
        Y = Y_treatment,
        ncomp = 2,
        test.keepX = test_keepX,
        scale = TRUE,
        nrepeat = 3
    )

    keepX_optimal <- tune_res$choice.keepX

    # conduct sparse PLS-DA
    set.seed(123)
    
    s_model <- mixOmics(
        X = X_intensities,
        Y = Y_treatment,
        study = study_cohort,
        ncomp = 2,
        keepX = keepX_optimal,
        mode = "regression",
        scale = TRUE,
        max.iter = 100
    )

    return(s_model)

}

splsda_neg <- do_splsda(.intensities = intensity_norm_log2_neg, .mode = "neg")
splsda_pos <- do_splsda(intensity_norm_log2_pos, "pos")

# function to export tables and plots per mode
export_splsda <- function(.s_model, .mode) {

    .s_model$variates$X |> # export full sample information; variates
        as.data.frame() |>
        rownames_to_column(var = "sample") |>
        write.csv(file = file.path("data_processed", paste0("04_plsda_variates_", .mode, ".csv")), row.names = FALSE)
    
    .s_model$loadings$X |> # export full feature information; loadings
        as.data.frame() |>
        rownames_to_column(var = "mz_rt_min") |>
        write.csv(file = file.path("data_processed", paste0("04_plsda_sparse_loadings_", .mode, ".csv")), row.names = FALSE)
    
    .s_model$loadings$Y |> # export loadings associations with treatment groups
        as.data.frame() |>
        rownames_to_column(var = "treatment") |>
        write.csv(file = file.path("output", "tables", paste0("04_plsda_loadings_association_", .mode, ".csv")), row.names = FALSE)
    
    # prepare retained features for export
    select_features <- full_join(
        x = selectVar(.s_model, comp = 1)$value |> # selectVar can only act on one component at a time
            as.data.frame() |> rownames_to_column(var = "mz_rt_min") |> rename(comp1_value = value.var),
        y = selectVar(.s_model, comp = 2)$value |> 
            as.data.frame() |> rownames_to_column(var = "mz_rt_min") |> rename(comp2_value = value.var),
        by = "mz_rt_min"
    )

    write.csv(
        select_features,
        file = file.path("output", "tables", paste0("04_plsda_select_features_", .mode, ".csv")),
        row.names = FALSE
    )

    # export mixOmics' default plots. Need to use traditional R object setup and dev.off routine
    png( # plot the samples in PLS-DA space
        filename = file.path("output", "plots", paste0("04_plsda_samples_", .mode, ".png")),
        width = 6, height = 6, units = "in", res = 300
    )

    plotIndiv(
        object = .s_model,
        comp = c(1,2),
        col = viridis::viridis(n = 4, option = "viridis", direction = -1),
        ellipse = TRUE,
        centroid = TRUE,
        title = paste0("sparse PLS-DA ", .mode, "-mode"),
        legend = TRUE,
        legend.title = "treatment"
    )

    dev.off()

    png( # plot features in PLS-DA space
        filename = file.path("output", "plots", paste0("04_plsda_features_", .mode, ".png")),
        width = 6, height = 6, units = "in", res = 300
    )
    
    plotVar(
        object = .s_model,
        title = paste0("sparse PLS-DA features ", .mode, "-mode")
    )

    dev.off()


    png( # plot the association between features and component 1
        filename = file.path("output", "plots", paste0("04_plsda_comp1_", .mode, ".png"))
    )

    plotLoadings(
        object = .s_model,
        comp = 1, # plotLoadings can only act on one component at a time
        title = paste0("sparse PLS-DA comp1 | ", .mode, "-mode"),
        size.title = 1
    )

    dev.off()

    png( # plot the association between features and component 2
        filename = file.path("output", "plots", paste0("04_plsda_comp2_", .mode, ".png"))
    )

    plotLoadings(
        object = .s_model,
        comp = 2,
        title = paste0("sparse PLS-DA comp2 | ", .mode, "-mode"),
        size.title = 1
    )

    dev.off()

}

export_splsda(.s_model = splsda_neg, .mode = "neg")
export_splsda(splsda_pos, "pos")

## 01) Organizing both metabolomics datasets and establishing groups.
if(!interactive()) pdf(NULL) # prevents pdf artifact from being generated from plots

# load raw LC-MS/MS peaks data
intensity_negative_raw <- read.csv(file.path("data_raw", "intensity_negative.csv"), header = TRUE)
intensity_positive_raw <- read.csv(file.path("data_raw", "intensity_positive.csv"), header = TRUE)

# merge duplicate peaks (same mz_ratio and rt_min) by summing intensities
intensity_negative <- intensity_negative_raw |>
    group_by(mz_ratio, rt_min) |>
    summarize(
        across(
          .cols = where(is.numeric),
          .fns = \(x) sum(x, na.rm = TRUE)
        ),
        across(
          .cols = -where(is.numeric),
          .fns = ~ first(.x)
        ),
        .groups = "drop"
    ) |>
    mutate(mz_rt_min = paste0(mz_ratio, "_", rt_min)) |>
    relocate(mz_rt_min, calculated_mw, mz_ratio, rt_min, reference_ion)

intensity_positive <- intensity_positive_raw |>
    group_by(mz_ratio, rt_min) |>
    summarize(
        across(
          .cols = where(is.numeric),
          .fns = \(x) sum(x, na.rm = TRUE)
        ),
        across(
          .cols = -where(is.numeric),
          .fns = ~ first(.x)
        ),
        .groups = "drop"
    ) |>
    mutate(mz_rt_min = paste0(mz_ratio, "_", rt_min)) |>
    relocate(mz_rt_min, calculated_mw, mz_ratio, rt_min, reference_ion)

# export information describing the LC-MS/MS peaks
feature_info <- c("calculated_mw", "mz_ratio", "rt_min", "reference_ion") # identify feature metadata columns

feature_metadata <- bind_rows( # stacks the feature info for both sets
    select(intensity_negative, all_of(c("mz_rt_min", feature_info))) |> # grabs feature info from negative set
        mutate(mode = "negative"), # adds an identifier for the LC-MS/MS negative-mode set
    select(intensity_positive, all_of(c("mz_rt_min", feature_info))) |> 
        mutate(mode = "positive")
)

write.csv(feature_metadata, file = file.path("data_raw", "feature_metadata.csv"), quote = TRUE, row.names = FALSE) # export feature metadata

# perform median normalization followed by log2 transformation
sample_medians_neg <- intensity_negative |>
    select(!any_of(feature_info)) |> # ignore feature metadata
    summarize(across(where(is.numeric), \(x) median(x, na.rm = TRUE))) # grab median peak from each sample

sample_medians_pos <- intensity_positive |>
    select(!any_of(feature_info)) |> 
    summarize(across(where(is.numeric), \(x) median(x, na.rm = TRUE)))

target_median_neg <- sample_medians_neg |> 
    summarize(across(everything(), median)) |> # find the median across all (negative) LC-MS/MS peaks
    pull(1) # grab the value, not the column

target_median_pos <- sample_medians_pos |> 
    summarize(across(everything(), median)) |>
    pull(1)

intensity_norm_log2_neg <- intensity_negative |>
    select(!any_of(feature_info)) |> # remove feature metadata
    mutate(across(
        .cols = where(is.numeric), # only target sample columns, not mz_rt_min
        .fns = ~ .x / sample_medians_neg[[cur_column()]] * target_median_neg # applies median normalization
    )) |>
    mutate(across(
        .cols = where(is.numeric), # agian, only target sample columns
        .fns = ~ log2(.x + 1e-6) # add small offset to avoid instances of log2(0)
    ))

intensity_norm_log2_pos <- intensity_positive |>
    select(!any_of(feature_info)) |>
    mutate(across(
        .cols = where(is.numeric),
        .fns = ~ .x / sample_medians_pos[[cur_column()]] * target_median_pos
    )) |>
    mutate(across(
        .cols = where(is.numeric),
        .fns = ~ log2(.x + 1e-6)
    ))

# export for inspection
write.csv(intensity_norm_log2_neg, file = file.path("data_processed", "01_intensity_norm_log2_neg.csv"), row.names = FALSE)
write.csv(intensity_norm_log2_pos, file = file.path("data_processed", "01_intensity_norm_log2_pos.csv"), row.names = FALSE)

# plot comparisons of sample-wise data before and after normalization & transformation
if(!interactive()) pdf(NULL) # prevents pdf artifact from being generated from plots

intensity_long_log2 <- bind_rows( # gather the raw and normalized-log2 data into a single object
    intensity_negative |> select(!any_of(feature_info)) |> mutate(mode = "negative", stage = "raw"), # specify mode and raw vs normalized
    intensity_positive |> select(!any_of(feature_info)) |> mutate(mode = "positive", stage = "raw"),
    intensity_norm_log2_neg |> mutate(mode = "negative", stage = "median-norm"),
    intensity_norm_log2_pos |> mutate(mode = "positive", stage = "median-norm")
) |>
    pivot_longer( # collapse the various sample columns into a long format
        cols = !c(stage, mode, mz_rt_min), # new col.names are stage, mode, mz_rt_min, sample, and intensity
        names_to = "sample",
        values_to = "intensity"
    ) |>
    mutate(
        stage = factor(stage, levels = c("raw", "median-norm")),
        intensity_log2 = ifelse( # convert raw data to log2 scale
            stage == "raw",
            yes = log2(intensity), # convert if raw
            no = intensity # leave unmodified if already log2 ("median-norm")
        )
    )

write.csv( # export the long-data formatted table for easy re-plotting in the Quarto report
    intensity_long_log2,
    file = file.path("data_processed", "01_intensity_long_log2.csv"),
    row.names = FALSE
)

p_intensity_distribution <- intensity_long_log2 |> # create boxplots of the data
    ggplot(aes(x = intensity_log2, y = sample)) +
        geom_boxplot(outlier.size = 0.5) +
        facet_wrap(~stage + mode) + # calculate medians and IQR by stage and mode
        expand_limits(x = 0) +
        scale_x_continuous(expand = c(0,0)) +
        theme(
            strip.text = element_text(face = "bold"),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title = element_text(face = "bold")
        ) +
        labs(x = "log2 intensity")

cowplot::save_plot(
    filename = file.path("output", "plots", "01_sample_distribution.png"),
    plot = p_intensity_distribution,
    bg = "white",
    base_height = 8,
    base_width = 8 / 1.618
)

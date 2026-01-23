## 01) Organizing both metabolomics datasets and establishing groups.
# load raw LC-MS/MS peaks data
peaks_negative_raw <- read.csv("data_raw/peaks_negative.csv", header = TRUE)
peaks_positive_raw <- read.csv("data_raw/peaks_positive.csv", header = TRUE)

# merge duplicate peaks (same mz_ratio and rt_min) by summing intensities
peaks_negative <- peaks_negative_raw |>
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
    relocate(mz_rt_min)

peaks_positive <- peaks_positive_raw |>
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
    relocate(mz_rt_min)

# export feature metadata
feature_info <- c("calculated_mw", "mz_ratio", "rt_min", "reference_ion") # identify feature metadata columns

feature_metadata <- bind_rows( # stacks the feature info for both sets
    select(peaks_negative, all_of(c("mz_rt_min", feature_info))) %>% # grabs feature info from negative set
        mutate(mode = "negative"), # adds an identifier for the LC-MS/MS negative-mode set
    select(peaks_positive, all_of(c("mz_rt_min", feature_info))) %>% 
        mutate(mode = "positive"),
)

write.csv(feature_metadata, file = "data_raw/feature_metadata.csv", quote = TRUE, row.names = FALSE) # export feature metadata

# apply log2 transformation to peak area, then export
peaks_negative_log2 <- peaks_negative |>
    select(!any_of(feature_info)) |> # remove feature metadata
    mutate(across(
        .cols = where(is.numeric), # only target remaining numeric columns
        .fns = ~ log2(.x + 1e-6)) # add small offset to avoid instances of log2(0)
    )
 
peaks_positive_log2 <- peaks_positive |>
    select(!any_of(feature_info)) |>
    mutate(across(
        .cols = where(is.numeric),
        .fns = ~ log2(.x + 1e-6))
    )

write.csv(peaks_negative_log2, file = "data_processed/peaks_negative_log2.csv", row.names = FALSE)
write.csv(peaks_positive_log2, file = "data_processed/peaks_positive_log2.csv", row.names = FALSE)
## 01) Organizing both metabolomics datasets and establishing groups.
# load raw LC-MS/MS peaks data
peaks_negative <- read.csv("data_raw/peaks_negative.csv", header = TRUE)
peaks_positive <- read.csv("data_raw/peaks_positive.csv", header = TRUE)

# prepare feature metadata
feature_info <- c("calculated_mw", "mz_ratio", "rt_min", "reference_ion") # identify feature metadata columns

feature_metadata <- bind_rows( # stacks the feature info for both sets
    select(peaks_negative, all_of(c(mz_rt_min, feature_info))) %>% # grabs feature info from negative set
        mutate(mode = "negative"), # adds an identifier for the LC-MS/MS negative-mode set
    select(peaks_positive, all_of(c(mz_rt_min, feature_info))) %>% 
        mutate(mode = "positive"),
)

write.csv(feature_metadata, file = "data_processed/feature_metadata.csv", quote = TRUE, row.names = FALSE) # export feature metadata

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
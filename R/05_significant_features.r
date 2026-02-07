## 05) Organize features selected by limma and PLS-DA for putative KEGG annotation
feature_metadata <- read.csv(file.path("data_raw", "feature_metadata.csv"), header = TRUE) # call in manually in case running script separately

# function to pull feature IDs from each limma contrast
feature_files <- list.files(
    path = file.path("output", "tables"), # set file path to find the features
    pattern = "^(03_limma_tr|04_plsda_select_features).*\\.csv$", # grab features from limma and plsda outputs
    full.names = TRUE
)

read_path <- function(.path) {
    read.csv(.path, header = TRUE) |>
    mutate(source_file = basename(.path), mz_rt_min = as.character(mz_rt_min)) |>
    select(source_file, mz_rt_min)
}

features_to_annotate <- lapply(
    X = feature_files,
    FUN = read_path
) |>
    bind_rows() |>
    left_join(y = select(feature_metadata, mz_rt_min, calculated_mw, mode), by = "mz_rt_min") |> # add in feature information
    select(-source_file) |>
    distinct() |>
    arrange(calculated_mw, mode)

write.csv(features_to_annotate, file = file.path("data_processed", "05_features_to_annotate.csv"), row.names = FALSE, quote = FALSE)

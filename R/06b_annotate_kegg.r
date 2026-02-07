## 06b) Use locally-stored KEGG compound info to putatively annotate significant features from the data
# grab the features to be annotated and the KEGG info
features_to_annotate <- read.csv(file.path("data_processed", "05_features_to_annotate.csv"), header = TRUE)

kegg_compounds <- read.table(file.path("data_processed", "06_kegg_compounds.tsv"), sep = "\t", header = TRUE) |>
    as_tibble()

# establish criteria for matching compound info to calculated_mw from the LC-MS features
ppm_tolerance <- 10 # establish a 10 parts per million tolerance between feature calcualted molecular weight and KEGG exact mass

features_to_annotate <- features_to_annotate |>
    distinct(mz_rt_min, .keep_all = TRUE) |>
    dplyr::filter(!is.na(calculated_mw)) |>
    mutate(
        mw_min = calculated_mw * (1 - ppm_tolerance / 1e6),
        mw_max = calculated_mw * (1 + ppm_tolerance / 1e6)
    )

features_annotated <- tidyr::crossing( # generates massive combinations of cases & variables from each object (3,808,310 rows)
    features_to_annotate,
    kegg_compounds
) |>
    filter(exact_mass >= mw_min, exact_mass <= mw_max) |> # retain only matches that fall within mass/mw tolerance range
    group_by(mz_rt_min) |> # group to evalute features with multiple matches
    mutate(
        n_matches = n(), # tracks number of times a feature appears
        annotation = if_else(n_matches == 1, "putative_single", "putative_multi") # labels feature on their number of annotations
    ) |>
    ungroup() |>
    select(
        mz_rt_min,
        calculated_mw,
        exact_mass,
        kegg_id,
        name,
        formula,
        pathways,
        annotation
    ) # 307 putative annotations

write.table(
    features_annotated,
    file = file.path("output", "tables", "06_putative_annotations.tsv"),
    sep = "\t",
    row.names = FALSE
)
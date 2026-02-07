## 06a) Accessing KEGG for compound information to be stored locally
## This script is usually skipped in run_all.r after being run once, as it can be time-consuming to access KEGG
compounds <- read.table( # download and format KEGG's full compound list
    "https://rest.kegg.jp/list/compound",
    sep = "\t", quote = "", stringsAsFactors = FALSE
) |>
    rename(kegg_id = V1, name = V2) |>
    as_tibble()

id_chunks <- split( # need to break up KEGG submissions into chunks to respect online access limits
    compounds$kegg_id, # grab compound IDs
    ceiling(seq_along(compounds$kegg_id) / 10L) # divide compounds into groups; keggGet only accepts 10 entries at a time
)

get_cpd_details <- function(ids) {
    Sys.sleep(0.5) # create delays between KEGG requests
    keggGet(ids) # query KEGG for ID info
}

cat("\nBeginning KEGG query for compound masses. This will take a while\n")

details_list <- lapply(id_chunks, get_cpd_details) # apply keggGet() in chunks

details_flat <- unlist(details_list, recursive = FALSE) # flatten one level, but keep a list of entry objects

kegg_compounds_list <- lapply( # parse each entry into a data.frame, within a list
    details_flat,
    function(x) {
        data.frame(
            kegg_id = x$ENTRY, # the compound ID, e.g., cpd:C00001
            name = paste(x$NAME, collapse = "; "), # compound names, compressed format
            formula = if (!is.null(x$FORMULA)) x$FORMULA else NA_character_, # extract formula if present
            exact_mass = if (!is.null(x$EXACT_MASS)) as.numeric(x$EXACT_MASS) else NA_real_, # extract exact mass if present
            pathways = if (!is.null(x$PATHWAY)) paste(names(x$PATHWAY), collapse = ";") else NA_character_, # extract pathways list if present
            stringsAsFactors = FALSE
        )
  }
)

kegg_compounds <- bind_rows(kegg_compounds_list) # collapse list into a data.frame/tibbble. The info is now easier to work with

kegg_compounds <- kegg_compounds |> # convoluted extra step to clean up names and pathways lists. Not strictly necessary
    mutate(
        name = gsub(";+", ";", name),
        name = gsub("^;\\s*|\\s*;$", "", name),
        name = trimws(name),
        pathways = ifelse(is.na(pathways), pathways,
            gsub(";+", ";", pathways)),
        pathways = ifelse(is.na(pathways), pathways,
            gsub("^;\\s*|\\s*;$", "", pathways)),
        pathways = ifelse(is.na(pathways), pathways, trimws(pathways))
    ) |>
    distinct(kegg_id, .keep_all = TRUE)

write.table( # export the KEGG compound information to compare against the data in next script
    kegg_compounds,
    file = file.path("data_processed", "06_kegg_compounds.tsv"), # saved as .tsv to avoid conflicts with commas in chemical names
    sep = "\t",
    row.names = FALSE
)

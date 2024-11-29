library(tidyverse)
library(yaml)


config <- yaml::read_yaml("_pipelines/validation_dataset/config/config-fragmentation_benchmark.yaml")

samples <- config$samples
samples <- samples %>%
    names() %>%
    as_tibble() %>%
    rename(sample = value)

df <- samples %>%
    mutate(
        motifs = map(.x = sample, ~ read_delim(
            paste(
                "data/monocultures", .x, "nanomotif_0.4.13", "bin-motifs.tsv",
                sep = "/"
            ),
            delim = "\t",
            col_types = cols(.default = col_character())
        ))
    ) %>%
    unnest(motifs) %>%
    select(-sample)

write_delim(df, "bin-motifs.tsv", delim = "\t")

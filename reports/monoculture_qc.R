library(data.table)
library(tidyverse)
library(ggtext)
library(seqinr)
library(yaml)


config <- yaml::read_yaml("_pipelines/validation_dataset/config/config-fragmentation_benchmark.yaml")

assembly <- read.fasta("output/monocultures_qc/concatenated_assembly.fasta", seqtype = "DNA")
genome_cov <- fread("output/monocultures_qc/concatenated_coverage.txt", header = FALSE)


assembly_info <- data.frame(
    contig = names(assembly),
    length = sapply(assembly, function(seq) {
        length(seq)
    }),
    gc = sapply(assembly, function(seq) {
        seqinr::GC(seq)
    })
)
rm(assembly)

assembly_info <- assembly_info %>%
    as_tibble() %>%
    mutate(
        bin = str_remove(contig, "_contig.*")
    )


genome_cov <- genome_cov %>%
    rename(
        contig = V1,
        cov = V4
    ) %>%
    select(contig, cov) %>%
    group_by(contig) %>%
    summarise(
        mean_cov = mean(cov),
        sd_cov = sd(cov)
    ) %>%
    mutate(
        bin = str_remove(contig, "_contig.*")
    )



genome_cov %>%
    ggplot(aes(x = bin, y = mean_cov, color = bin)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(alpha = 0.5) +
    theme_minimal() +
    scale_y_log10() +
    # geom_hline(yintercept = 40, linetype = "dashed", color = "red") +
    coord_flip() +
    labs(
        x = "Bin",
        y = "Mean coverage",
        title = "Mean coverage per bin"
    ) +
    theme(
        legend.position = "none",
        axis.text.y = element_markdown()
    )




genome_cov_contig %>%
    filter(mean_cov < 100) %>%
    ggplot(aes(x = bin_label, y = mean_cov, color = bin)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(alpha = 0.5) +
    theme_minimal() +
    geom_hline(yintercept = 40, linetype = "dashed", color = "red") +
    coord_flip() +
    labs(
        x = "Bin",
        y = "Mean coverage",
        title = "Mean coverage per bin"
    ) +
    theme(
        legend.position = "none",
        axis.text.y = element_markdown()
    )




genome_cov_gc <- genome_cov %>%
    left_join(assembly_info, by = "contig")

genome_cov_gc %>%
    ggplot(aes(x = gc, y = mean_cov, color = bin, size = length)) +
    geom_point(alpha = 0.5) +
    theme_minimal() +
    scale_y_log10() +
    labs(
        x = "GC content",
        y = "Mean coverage",
        title = "Mean coverage per GC content"
    ) +
    theme(
        legend.position = "none"
    )

genome_cov_gc$bin %>% unique()
genome_cov_gc %>%
    filter(str_detect(bin, "CVM73")) %>%
    mutate(
        contig_nr = str_extract(contig, "contig_\\d+")
    ) %>%
    group_by(contig_nr) %>%
    summarise(
        mean_cov = mean(mean_cov),
        length = sum(length)
    ) %>%
    arrange(mean_cov)





genome_cov_gc %>%
    mutate(
        group = case_when(
            str_detect(bin, "CVM77") ~ "CVM77",
            TRUE ~ "other"
        )
    ) %>%
    ggplot(aes(x = gc, y = mean_cov, size = length, color = group)) +
    geom_point(alpha = 0.5) +
    theme_minimal() +
    scale_y_log10() +
    labs(
        x = "GC content",
        y = "Mean coverage",
        title = "Mean coverage per GC content"
    ) +
    theme(
        legend.position = "none"
    )

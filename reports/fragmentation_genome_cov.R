library(data.table)
library(tidyverse)
library(ggtext)
library(seqinr)

genome_cov <- fread("data/datasets/fragmentation_benchmark_20kb/genomecov.txt", header = FALSE)
contig_bin <- fread("data/datasets/fragmentation_benchmark_20kb/contig_bin.tsv", header = FALSE)
assembly <- read.fasta("data/datasets/fragmentation_benchmark_20kb/concatenated_assembly.fasta", seqtype = "DNA")

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
    as_tibble()

contig_bin <- contig_bin %>%
    rename(
        contig = V1,
        bin = V2
    )

genome_cov <- genome_cov %>%
    rename(
        contig = V1,
        cov = V4
    ) %>%
    select(contig, cov)


genome_cov_contig <- genome_cov %>%
    group_by(contig) %>%
    summarise(
        mean_cov = mean(cov),
        sd_cov = sd(cov)
    ) %>%
    mutate(
        bin = str_remove(contig, "_contig.*"),
        bin_label = case_when(
            !bin %in% contig_bin$bin ~ paste0("<span style='color: red;'>", bin, "</span>"),
            TRUE ~ bin
        )
    )


genome_cov_contig %>%
    ggplot(aes(x = bin_label, y = mean_cov, color = bin)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(alpha = 0.5) +
    theme_minimal() +
    scale_y_log10() +
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




genome_cov_gc <- genome_cov_contig %>%
    left_join(assembly_info, by = "contig") %>%
    left_join(contig_bin, by = c("bin", "contig")) %>%
    mutate(
        col = case_when(
            bin == "DSMZ11109-Desulfobacca_acetoxidans" ~ "DSMZ11109-Desulfobacca_acetoxidans",
            TRUE ~ "other"
        )
    )


genome_cov_gc %>%
    filter(bin == "DSMZ11109-Desulfobacca_acetoxidans") %>%
    ggplot(aes(x = gc, y = mean_cov, color = col)) +
    geom_point(alpha = 0.5) +
    theme_minimal() +
    scale_y_log10()


genome_cov_gc %>%
    filter(bin == "DSMZ11109-Desulfobacca_acetoxidans") %>%
    mutate(
        contig_nr = str_extract(contig, "contig_\\d+")
    ) %>%
    group_by(contig_nr) %>%
    summarise(
        mean_cov = mean(mean_cov),
        length = sum(length)
    ) %>%
    arrange(mean_cov)

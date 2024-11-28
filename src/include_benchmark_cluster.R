#!/bin/env Rscript
if (!require("argparse")) install.packages("argparse")
library(argparse)

# Create argument parser
parser <- ArgumentParser(description = "Create benchmark report")

parser$add_argument("--baselinedir", help = "Path to the semibin2 folder [default: output/baseline]", default = "output/baseline/")
parser$add_argument("--benchmarkdir", help = "Path to the semibin3 folder [default: output/benchmarks]", default = "output/benchmarks")
parser$add_argument("--benchmarks", help = "Comma-separated list of benchmarks", default = "", type = "character")
parser$add_argument("--samples", help = "Comma-separated list of samples. ['.' specifies all samples]", default = "", type = "character")
parser$add_argument("--output", required = TRUE, help = "Path to the output folder")

# Parse arguments

# Mock args
args <- list(
    baselinedir = "output/baseline/detect_contamination/",
    benchmarkdir = "output/benchmarks/include_contigs/2024-11-27_read-level_clustering_init",
    benchmarks = ".",
    samples = ".",
    output = "analysis/2024-11-27_clustering_include_init"
)

# Create output folder if not exists
if (!dir.exists(args$output)) {
    dir.create(args$output, recursive = TRUE)
}

# Load libraries
if (!require("data.table")) install.packages("data.table")
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("ggtext")) install.packages("ggtext")
if (!require("janitor")) install.packages("janitor")
if (!require("ggpubr")) install.packages("ggpubr")
library(data.table)
library(tidyverse)
library(ggtext)
library(janitor)
library(ggpubr)

custom_theme <- theme_bw() +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "right"
    )


# split arguments
benchmarks <- strsplit(args$benchmarks, ",")[[1]]
samples <- strsplit(args$samples, ",")[[1]]

# if benchmark or samples is . find all
if ("." %in% samples) {
    for (s in samples) {
        samples <- c(samples, list.files(file.path(args$benchmarkdir, s)))
    }
    # remove .
    samples <- samples[samples != "."]
}

if ("." %in% benchmarks) {
    for (s in samples) {
        benchmarks <- c(benchmarks, list.files(file.path(args$benchmarkdir, s)))
    }

    benchmarks <- unique(benchmarks)

    # remove .
    benchmarks <- benchmarks[benchmarks != "."]
}

benchmarks <- benchmarks[benchmarks != "include.done"]
# load file if exists
load_file <- function(file, has_header = TRUE) {
    if (file.exists(file)) {
        fread(file, header = has_header)
    } else {
        tibble()
    }
}


# load all benchmarks
df_benchmark <- crossing(
    sample = samples,
    benchmark = benchmarks
) %>%
    mutate(
        contig_bin = map2(sample, benchmark, ~ load_file(file.path(args$benchmarkdir, .x, .y, "contig_bin.tsv"), has_header = FALSE)),
        include = map2(sample, benchmark, ~ load_file(file.path(args$benchmarkdir, .x, .y, "include_contigs.tsv")) %>% mutate(bin = as.character(bin))),
    )


# Create bin truth
bin_truth <- tibble(sample = samples) %>%
    mutate(
        bin_truth = map(sample, ~ load_file(file.path("data/datasets/", .x, "contig_bin.tsv"), has_header = FALSE))
    ) %>%
    unnest(bin_truth) %>%
    rename(
        contig = V1,
        bin_truth = V2
    )

contig_bin <- df_benchmark %>%
    select(sample, benchmark, contig_bin) %>%
    unnest(contig_bin) %>%
    rename(
        contig = V1,
        bin = V2
    )

include <- df_benchmark %>%
    select(sample, benchmark, include) %>%
    mutate(
        tib_len = map_int(include, ~ nrow(.x))
    ) %>%
    filter(tib_len > 0) %>%
    select(-tib_len) %>%
    unnest(include) %>%
    select(sample, benchmark, contig, cluster, bin, assigned_bin, assignment_is_unique)

# Combine data
unbinned_contigs <- bin_truth %>%
    crossing(tibble(benchmark = benchmarks)) %>%
    left_join(contig_bin) %>%
    mutate(bin = ifelse(is.na(bin), "unbinned", bin)) %>%
    filter(bin == "unbinned")


d_unique <- unbinned_contigs %>% # only unbinned
    left_join(include %>% filter(assignment_is_unique) %>% select(sample:contig, assigned_bin)) %>%
    mutate(
        measure = case_when(
            bin_truth == assigned_bin ~ "correct",
            bin_truth != assigned_bin ~ "incorrect",
            is.na(assigned_bin) ~ "unassigned",
            TRUE ~ "unknown"
        ),
        assignment_type = "unique"
    )

d_candidate <- unbinned_contigs %>%
    left_join(include %>% filter(!assignment_is_unique) %>% select(sample:contig, assigned_bin)) %>%
    filter(!is.na(assigned_bin)) %>%
    mutate(
        measure = case_when(
            bin_truth == assigned_bin ~ "correct",
            bin_truth != assigned_bin ~ "incorrect",
            TRUE ~ "unknown"
        )
    ) %>%
    group_by(sample, benchmark, contig, measure) %>%
    summarise(
        n_bins = n()
    )



# Plot confusion matrix
confusion_df <- d_unique %>%
    group_by(measure, sample, benchmark) %>%
    summarise(
        count = n()
    )

metrics <- confusion_df %>%
    select(sample, benchmark, measure, count) %>%
    pivot_wider(names_from = measure, values_from = count) %>%
    # fill NA with 0
    mutate_all(~ replace(., is.na(.), 0)) %>%
    clean_names() %>%
    mutate(
        misclassification_rate = incorrect / (correct + incorrect),
        precision = correct / (correct + incorrect),
        recall = correct / (correct + unassigned),
    ) %>%
    pivot_longer(
        cols = correct:recall,
        names_to = "metric",
        values_to = "value"
    )


METRICS <- c("misclassification_rate", "precision", "recall")

metric_labels <- metrics %>%
    filter(metric %in% METRICS) %>%
    group_by(sample, metric) %>%
    summarise(
        value = mean(value, na.rm = TRUE)
    ) %>%
    mutate(
        label = paste0(round(value, 2))
    )

m1 <- metrics %>%
    filter(metric %in% METRICS) %>%
    mutate(
        metric = factor(metric, levels = METRICS)
    ) %>%
    ggplot(aes(x = metric, y = value, color = sample)) +
    geom_boxplot(outlier.alpha = 0) +
    geom_point(position = position_jitterdodge(jitter.width = 0.15), alpha = 0.6) +
    geom_text(data = metric_labels, aes(label = label, y = 1.05), position = position_dodge(width = 1)) +
    # facet_wrap(~sample, scales = "free") +
    # limit y axis
    scale_y_continuous(limits = c(0, 1.05)) +
    # scale_color_manual(values = c("nanomotif" = "#3ab7ff", "developement benchmark" = "#196900")) +
    custom_theme +
    labs(
        y = "Value"
    )

m2 <- metrics %>%
    filter(!metric %in% METRICS) %>%
    ggplot(aes(x = metric, y = value, color = sample)) +
    geom_boxplot(outlier.alpha = 0) +
    geom_point(position = position_jitterdodge(jitter.width = 0.15), alpha = 0.6) +
    # facet_wrap(~sample, scales = "free") +
    # scale_color_manual(values = c("nanomotif" = "#3ab7ff", "developement benchmark" = "#196900")) +
    # limit y axis
    custom_theme +
    labs(
        y = "Count"
    )

ggarrange(m1, m2, ncol = 2, nrow = 1, common.legend = TRUE, legend = "bottom")

ggsave(file.path(args$output, "metrics_boxplot.png"), width = 10, height = 5)


for (s in samples) {
    motifs_scored <- read_delim(file.path("output/baseline/include_files", s, "motifs-scored-read-methylation-1.tsv"), "\t")


    motifs_scored_f <- motifs_scored %>%
        filter(mean_read_cov > 2) %>%
        filter((N_motif_obs * mean_read_cov) >= 40) %>%
        left_join(bin_truth %>% filter(sample == s)) %>%
        rename(bin = bin_truth) %>%
        mutate(
            motif_mod = paste0(motif, "_", mod_type, "_", mod_position)
        ) %>%
        filter(!is.na(bin)) %>%
        group_by(bin, motif_mod) %>%
        summarise(
            median = mean(median)
        ) %>%
        mutate(
            group = "methylation"
        )

    df_wide <- motifs_scored_f %>%
        select(bin, motif_mod, median) %>%
        pivot_wider(names_from = motif_mod, values_from = median, values_fill = list(median = 0))

    # Convert to matrix
    data_matrix <- as.matrix(df_wide[, -1])
    rownames(data_matrix) <- df_wide$bin

    # Perform clustering
    dist_rows <- dist(data_matrix, method = "euclidean")
    dist_cols <- dist(t(data_matrix), method = "euclidean")
    hc_rows <- hclust(dist_rows, method = "average")
    hc_cols <- hclust(dist_cols, method = "average")

    mod_cols <- colnames(data_matrix)[hc_cols$order]
    bin_rows <- rownames(data_matrix)[hc_rows$order]

    cross <- crossing(
        measure = c("correct", "incorrect", "unassigned"),
        sample = s,
        bin_truth = bin_rows
    ) %>%
        arrange(sample, bin_truth)

    confusion_df <- d_unique %>%
        filter(sample == s) %>%
        group_by(measure, sample, bin_truth) %>%
        summarise(
            count = n()
        ) %>%
        right_join(cross) %>%
        arrange(sample, bin_truth, measure) %>%
        mutate(
            count = ifelse(is.na(count), 0, count)
        ) %>%
        rename(bin = bin_truth)

    # any na in confusion_df

    MIN <- 0.00
    MAX <- 1.00
    metrics <- confusion_df %>%
        ungroup() %>%
        pivot_wider(names_from = measure, values_from = count) %>%
        # fill NA with 0
        mutate_all(~ replace(., is.na(.), 0)) %>%
        mutate(
            misclassification_rate = incorrect / (correct + incorrect),
            precision = correct / (correct + incorrect),
            recall = correct / (correct + unassigned),
        ) %>%
        pivot_longer(
            cols = correct:recall,
            names_to = "metric",
            values_to = "value"
        ) %>%
        group_by(bin, metric) %>%
        summarise(
            value = mean(value, na.rm = TRUE)
        ) %>%
        mutate(
            group = "metric",
            label = paste0(round(value * 100, 0)),
            value = (value - MIN) / (MAX - MIN)
        ) %>%
        rename(x_axis = metric) %>%
        filter(x_axis %in% METRICS) %>%
        filter(!is.na(value))


    motifs_scored_f %>%
        rename(x_axis = motif_mod, value = median) %>%
        bind_rows(metrics) %>%
        mutate(
            bin = factor(bin, levels = bin_rows),
            x_axis = factor(x_axis, levels = c(mod_cols, METRICS))
        ) %>%
        ggplot(aes(x = x_axis, y = fct_rev(bin), fill = value)) +
        geom_tile(color = "gray30") +
        geom_text(data = metrics, aes(label = label)) +
        scale_fill_gradient2(low = "white", high = "blue") +
        theme_minimal() +
        facet_grid(~group, scales = "free", space = "free") +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1)
        )

    ggsave(file.path(args$output, paste0(s, "_motif_score.png")), width = 35, height = 12)
}


# false_positive_contigs <- d_combined %>%
#     filter(sample == s, str_detect(contig, "DSMZ11109")) %>%
#     filter(actual_contaminant == "contamination") %>%
#     pull(contig) %>%
#     unique()

# motifs_scored %>%
#     left_join(bin_truth %>% filter(sample == s)) %>%
#     rename(bin = bin_truth) %>%
#     filter(bin == "DSMZ11109-Desulfobacca_acetoxidans") %>%
#     filter(contig %in% false_positive_contigs) %>%
#     mutate(
#         motif_mod = paste0(motif, "_", mod_type, "_", mod_position),
#         contig = case_when(
#             contig %in% false_positive_contigs ~ paste0("<span style='color:red'>", contig, "</span>"),
#             TRUE ~ contig
#         )
#     ) %>%
#     filter(mean_read_cov > 2) %>%
#     filter((N_motif_obs * mean_read_cov) >= 40) %>%
#     ggplot(aes(x = motif_mod, y = contig, fill = median)) +
#     geom_tile(color = "gray30") +
#     scale_fill_gradient2(low = "white", high = "blue") +
#     theme_minimal() +
#     theme(
#         axis.text.x = element_text(angle = 90, hjust = 1),
#         axis.text.y = element_markdown()
#     )

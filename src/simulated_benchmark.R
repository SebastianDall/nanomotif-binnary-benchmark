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
    benchmarkdir = "output/benchmarks/detect_contamination/2024-11-07_read-level_minMeth0.5_n_motif_obs=8",
    benchmarks = ".",
    samples = ".",
    output = "analysis/2024-11-07_read-level_simulated_4_lognormal_20-100x_minMeth0.5_n_motif_obs=8"
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
        mode = "developement benchmark",
        contig_bin = map2(sample, benchmark, ~ load_file(file.path(args$benchmarkdir, .x, .y, "contig_bin.tsv"), has_header = FALSE)),
        bin_contamination = map2(sample, benchmark, ~ load_file(file.path(args$benchmarkdir, .x, .y, "bin_contamination.tsv")) %>% mutate(bin = as.character(bin))),
    )

df_baseline <- crossing(
    sample = samples,
    benchmark = benchmarks
) %>%
    mutate(
        mode = "nanomotif",
        contig_bin = map2(sample, benchmark, ~ load_file(file.path(args$baselinedir, .x, .y, "contig_bin.tsv"), has_header = FALSE)),
        bin_contamination = map2(sample, benchmark, ~ load_file(file.path(args$baselinedir, .x, .y, "bin_contamination.tsv")) %>% mutate(bin = as.character(bin))),
    )

df <- bind_rows(df_benchmark, df_baseline)


# Create bin truth
bin_truth <- df %>%
    filter(benchmark == "original_contig_bin") %>%
    select(sample, contig_bin) %>%
    unnest(contig_bin) %>%
    rename(
        contig = V1,
        bin_truth = V2
    ) %>%
    distinct()

# Ecoli was not assembled well and ended up in two bins which in this case would be considered contamination if placed in each others bin.
contig_bin <- df %>%
    select(sample, benchmark, mode, contig_bin) %>%
    unnest(contig_bin) %>%
    rename(
        contig = V1,
        bin = V2
    ) %>%
    filter(!str_detect(contig, "coli"))

contamination <- df %>%
    select(sample, benchmark, mode, bin_contamination) %>%
    mutate(
        tib_len = map_int(bin_contamination, ~ nrow(.x))
    ) %>%
    filter(tib_len > 0) %>%
    select(-tib_len) %>%
    unnest(bin_contamination) %>%
    mutate(
        binary_methylation_mismatch_score = case_when(
            is.na(binary_methylation_mismatch_score) ~ binary_methylation_missmatch_score,
            TRUE ~ binary_methylation_mismatch_score
        )
    ) %>%
    select(sample, benchmark, mode, contig, binary_methylation_mismatch_score, non_na_comparisons) %>%
    mutate(
        prediction = "contamination"
    )


# Combine data
d_combined <- contig_bin %>%
    left_join(bin_truth) %>%
    left_join(contamination) %>%
    mutate(
        prediction = ifelse(is.na(prediction), "not contamination", prediction),
        measure = case_when(
            prediction == "not contamination" & bin == bin_truth ~ "true negative",
            prediction == "not contamination" & bin != bin_truth ~ "false negative",
            prediction == "contamination" & bin == bin_truth ~ "false positive",
            prediction == "contamination" & bin != bin_truth ~ "true positive"
        )
    )



# Plot distribution of mismatch score of contaminants
d_combined %>%
    filter(prediction == "contamination") %>%
    group_by(sample, benchmark, mode, binary_methylation_mismatch_score, measure) %>%
    summarise(
        count = n()
    ) %>%
    ggplot(aes(x = as.factor(binary_methylation_mismatch_score), y = count, color = measure)) +
    geom_boxplot(outlier.alpha = 0) +
    geom_point(position = position_jitterdodge(), alpha = 0.5) +
    facet_grid(. ~ mode) +
    custom_theme +
    labs(
        x = "Binary methylation mismatch score",
        y = "Count"
    )




# Save plot
ggsave(file.path(args$output, "contamination_distribution.png"), width = 10, height = 5)


# Plot confusion matrix
confusion_df <- d_combined %>%
    group_by(measure, sample, benchmark, mode) %>%
    summarise(
        count = n()
    ) %>%
    mutate(
        count = ifelse(is.na(count), 0, count)
    ) %>%
    right_join(
        tibble(
            actual = c("contamination", "contamination", "not contamination", "not contamination"),
            predicted = c("contamination", "not contamination", "contamination", "not contamination"),
            measure = c("true positive", "false negative", "false positive", "true negative")
        )
    ) %>%
    mutate(
        count_label = paste0(measure, " (", count, ")")
    )

metrics <- confusion_df %>%
    select(sample, benchmark, mode, measure, count) %>%
    pivot_wider(names_from = measure, values_from = count) %>%
    clean_names() %>%
    mutate(
        sensitivity = true_positive / (true_positive + false_negative),
        specificity = true_negative / (true_negative + false_positive),
        accuracy = (true_positive + true_negative) / (true_positive + true_negative + false_positive + false_negative)
    ) %>%
    janitor::clean_names() %>%
    mutate(
        accuracy = (true_positive + true_negative) / (true_positive + true_negative + false_positive + false_negative),
        sensitivity = true_positive / (true_positive + false_negative),
        specificity = true_negative / (true_negative + false_positive),
        precision = true_positive / (true_positive + false_positive),
        F1 = 2 * (precision * sensitivity) / (precision + sensitivity),
        F0_5 = (1 + 0.5^2) * (precision * sensitivity) / (0.5^2 * precision + sensitivity)
    ) %>%
    pivot_longer(
        cols = false_negative:F0_5,
        names_to = "metric",
        values_to = "value"
    )


metrics %>%
    ggplot(aes(x = mode, y = value, color = mode)) +
    geom_boxplot(outlier.alpha = 0) +
    geom_point(position = position_jitterdodge()) +
    facet_wrap(~metric, scales = "free") +
    custom_theme +
    labs(
        x = "Sample",
        y = "Value",
        fill = "Metric"
    )
METRICS <- c("sensitivity", "specificity", "precision", "accuracy", "F1", "F0_5")

m1 <- metrics %>%
    filter(metric %in% METRICS) %>%
    mutate(
        metric = factor(metric, levels = METRICS)
    ) %>%
    ggplot(aes(x = metric, y = value, color = mode)) +
    geom_boxplot(outlier.alpha = 0) +
    geom_point(position = position_jitterdodge()) +
    facet_wrap(~sample, scales = "free") +
    # limit y axis
    scale_y_continuous(limits = c(0, 1)) +
    scale_color_manual(values = c("nanomotif" = "#3ab7ff", "developement benchmark" = "#196900")) +
    custom_theme +
    labs(
        y = "Value"
    )

m2 <- metrics %>%
    filter(!metric %in% METRICS) %>%
    ggplot(aes(x = metric, y = value, color = mode)) +
    geom_boxplot(outlier.alpha = 0) +
    geom_point(position = position_jitterdodge()) +
    facet_wrap(~sample, scales = "free") +
    scale_color_manual(values = c("nanomotif" = "#3ab7ff", "developement benchmark" = "#196900")) +
    # limit y axis
    custom_theme +
    labs(
        y = "Count"
    )

ggarrange(m1, m2, ncol = 2, nrow = 1, common.legend = TRUE, legend = "bottom")

ggsave(file.path(args$output, "metrics.png"), width = 10, height = 5)

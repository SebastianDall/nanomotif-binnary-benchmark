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
    benchmarkdir = "output/benchmarks/detect_contamination/2024-11-25_clustering_contamination_init",
    benchmarks = ".",
    samples = ".",
    output = "analysis/2024-11-25_clustering_contamination_init"
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

benchmarks <- benchmarks[benchmarks != "benchmarks.done"]
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
        bin_contamination = map2(sample, benchmark, ~ load_file(file.path(args$benchmarkdir, .x, .y, "bin_contamination.tsv")) %>% mutate(bin = as.character(bin))),
    )


# Create bin truth
bin_truth <- df_benchmark %>%
    group_by(sample) %>%
    filter(row_number() == 1) %>%
    select(sample, contig_bin) %>%
    unnest(contig_bin) %>%
    rename(
        contig = V1,
        bin_truth = V2
    ) %>%
    mutate(
        bin_truth = case_when(
            str_detect(contig, "contamination") ~ str_remove(contig, "contamination_"),
            TRUE ~ bin_truth
        ),
        bin_truth = str_remove(bin_truth, "_contig.*"),
        contig = str_remove(contig, "contamination_")
    ) %>%
    arrange(sample, contig)

contig_bin <- df_benchmark %>%
    select(sample, benchmark, contig_bin) %>%
    unnest(contig_bin) %>%
    rename(
        contig = V1,
        bin = V2
    ) %>%
    mutate(
        actual_contaminant = case_when(
            str_detect(contig, "contamination_") ~ "contamination",
            TRUE ~ "not contamination"
        ),
        contig = str_remove(contig, "contamination_")
    )

contamination <- df_benchmark %>%
    select(sample, benchmark, bin_contamination) %>%
    mutate(
        tib_len = map_int(bin_contamination, ~ nrow(.x))
    ) %>%
    filter(tib_len > 0) %>%
    select(-tib_len) %>%
    unnest(bin_contamination) %>%
    select(sample, benchmark, bin, contig, binary_methylation_mismatch_score, non_na_comparisons) %>%
    mutate(
        prediction = "contamination",
        contig = str_remove(contig, "contamination_")
    )

# Combine data
d_combined <- contig_bin %>%
    left_join(bin_truth, by = c("sample", "contig")) %>%
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
    group_by(sample, benchmark, binary_methylation_mismatch_score, measure) %>%
    summarise(
        count = n()
    ) %>%
    ggplot(aes(x = as.factor(binary_methylation_mismatch_score), y = count, color = measure)) +
    geom_boxplot(outlier.alpha = 0) +
    geom_point(position = position_jitterdodge(), alpha = 0.5) +
    facet_grid(. ~ sample) +
    custom_theme +
    labs(
        x = "Binary methylation mismatch score",
        y = "Count"
    )




# Save plot
ggsave(file.path(args$output, "contamination_distribution_allowed_mismatch1.png"), width = 10, height = 5)


# Plot confusion matrix
confusion_df <- d_combined %>%
    group_by(measure, sample, benchmark) %>%
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
    select(sample, benchmark, measure, count) %>%
    pivot_wider(names_from = measure, values_from = count) %>%
    # fill NA with 0
    mutate_all(~ replace(., is.na(.), 0)) %>%
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


METRICS <- c("sensitivity", "specificity", "precision", "accuracy", "F1", "F0_5")

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
    geom_point(position = position_jitterdodge()) +
    geom_text(data = metric_labels, aes(label = label, y = 1.05), position = position_dodge(width = 1)) +
    facet_wrap(~sample, scales = "free") +
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
    geom_point(position = position_jitterdodge()) +
    facet_wrap(~sample, scales = "free") +
    # scale_color_manual(values = c("nanomotif" = "#3ab7ff", "developement benchmark" = "#196900")) +
    # limit y axis
    custom_theme +
    labs(
        y = "Count"
    )

ggarrange(m1, m2, ncol = 2, nrow = 1, common.legend = TRUE, legend = "bottom")

ggsave(file.path(args$output, "metrics_dev_allowed_mismatch1.png"), width = 10, height = 5)


for (s in samples) {
    motifs_scored <- read_delim(file.path("output/baseline/contamination_files", s, "motifs-scored-read-methylation-3.tsv"), "\t")


    motifs_scored_f <- motifs_scored %>%
        filter(N_motif_obs >= 8) %>%
        left_join(bin_truth) %>%
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

    confusion_df <- d_combined %>%
        group_by(measure, sample, bin, benchmark) %>%
        summarise(
            count = n()
        ) %>%
        mutate(
            count = ifelse(is.na(count), 0, count)
        )

    # any na in confusion_df

    MIN <- 0.00
    MAX <- 1.00
    METRICS <- c("sensitivity", "specificity", "precision", "accuracy", "F1", "F0_5")
    metrics <- confusion_df %>%
        ungroup() %>%
        pivot_wider(names_from = measure, values_from = count) %>%
        # fill NA with 0
        mutate_all(~ replace(., is.na(.), 0)) %>%
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

    ggsave(file.path(args$output, "motif_score.png"), width = 20, height = 18)
}

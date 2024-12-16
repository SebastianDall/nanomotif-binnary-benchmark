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
    benchmarkdir = "output/benchmarks/detect_contamination/2024-12-12_clustering_contamination_ensemble_spectral_gmm_hdbscan_agg_frac0.85_filt24_pca90",
    benchmarks = ".",
    samples = ".",
    output = "analysis/2024-12-12_clustering_contamination_ensemble_spectral_gmm_hdbscan_agg_frac0.85_filt24_pca90"
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
source("src/colors.R")

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
    ) %>%
    mutate(
        sample = str_remove(sample, "fragmentation_")
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
    select(sample, benchmark, contig, bin) %>%
    distinct() %>%
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
            prediction == "not contamination" & actual_contaminant == "not contamination" ~ "true negative",
            prediction == "not contamination" & actual_contaminant == "contamination" ~ "false negative",
            prediction == "contamination" & actual_contaminant == "not contamination" ~ "false positive",
            prediction == "contamination" & actual_contaminant == "contamination" ~ "true positive"
        )
    )


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


METRICS <- c("sensitivity", "specificity", "precision", "accuracy", "F1")

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
    # geom_text(data = metric_labels, aes(label = label, y = 1.05), position = position_dodge(width = 1)) +
    # scale_y_continuous(limits = c(0, 1.05)) +
    labs(
        y = "",
        x = "",
        title = "Contamination detection",
        subtitle = "Benchmark metrics"
    ) +
    scale_y_continuous(limits = c(0, 1)) +
    # facet_wrap(~sample, scales = "free") +
    # limit y axis
    # scale_color_manual(values = c("nanomotif" = "#3ab7ff", "developement benchmark" = "#196900")) +
    custom_theme +
    guides(
        color = guide_legend(title = "", position = "inside")
    ) +
    theme(
        legend.position.inside = c(0.4, 0.25),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "transparent"),
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

ggsave(file.path(args$output, "metrics_boxplot_and_counts.png"), width = 15, height = 5)



metric_labels <- metrics %>%
    filter(metric %in% METRICS) %>%
    filter(sample == "benchmark_20kb_1600kb_w_motif_discovery") %>%
    group_by(sample, metric) %>%
    summarise(
        value = mean(value, na.rm = TRUE)
    ) %>%
    mutate(
        label = paste0(round(value, 2))
    )

metrics %>%
    filter(metric %in% METRICS) %>%
    filter(sample == "benchmark_20kb_1600kb_w_motif_discovery") %>%
    mutate(
        metric = factor(metric, levels = METRICS)
    ) %>%
    ggplot(aes(x = metric, y = value)) +
    geom_boxplot(outlier.alpha = 0) +
    geom_point(position = position_jitter(width = 0.15), alpha = 0.6) +
    geom_text(data = metric_labels, aes(label = label, y = 1.05), position = position_dodge(width = 1)) +
    scale_y_continuous(limits = c(0, 1.05)) +
    labs(
        y = "",
        x = "",
        title = "Contamination detection benchmark"
    ) +
    custom_theme +
    guides(
        color = guide_legend(title = "", position = "inside")
    ) +
    theme(
        legend.position.inside = c(0.4, 0.25),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "transparent"),
    )
ggsave(file.path(args$output, "figure1_contamination_benchmark.png"), width = 4, height = 4, bg = "transparent")

samples <- str_remove(samples, "fragmentation_")
for (s in samples) {
    motifs_scored <- read_delim(file.path("output/baseline/contamination_files", paste0("fragmentation_", s), "motifs-scored-read-methylation-3.tsv"), "\t")


    motifs_scored_f <- motifs_scored %>%
        filter((N_motif_obs * mean_read_cov) >= 24) %>%
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
        ) %>%
        separate(motif_mod, into = c("motif", "mod_type", "mod_position"), convert = T) %>%
        mutate(
            ns = str_extract(motif, "NNN*"),
            ns = str_length(ns),
            motif_axis = paste0(
                str_sub(motif, 1, mod_position),
                "<strong>",
                str_sub(motif, mod_position + 1, mod_position + 1),
                "<sub>",
                map(mod_type, function(x) MOD_TYPE_PRETTY[[x]]) %>% unlist(),
                "</sup>",
                "</strong>",
                str_sub(motif, mod_position + 2, -1)
            ),
            motif_mod = str_replace(motif_axis, "NNN*", paste0("(N)<sub>", ns, "</sub>"))
        )

    motifs <- motifs_scored_f %>%
        filter(!is.na(median)) %>%
        filter(median > 0.5) %>%
        pull(motif_mod) %>%
        unique()

    df_wide <- motifs_scored_f %>%
        filter(motif_mod %in% motifs) %>%
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

    measures <- crossing(
        measure = c("true positive", "false negative", "false positive", "true negative"),
        sample = s,
        bin = bin_rows,
        benchmark = benchmarks
    )


    sample_data <- d_combined %>%
        filter(sample == s) %>%
        group_by(measure, sample, bin_truth, benchmark) %>%
        summarise(
            count = n()
        ) %>%
        arrange(bin_truth, benchmark) %>%
        rename(bin = bin_truth)

    confusion_df <- measures %>%
        left_join(sample_data) %>%
        arrange(benchmark, bin) %>%
        mutate(
            count = ifelse(is.na(count), 0, count)
        )

    # any na in confusion_df

    MIN <- 0.00
    MAX <- 1.00
    METRICS <- c("sensitivity", "specificity", "precision", "accuracy", "F1")
    metrics <- confusion_df %>%
        ungroup() %>%
        pivot_wider(names_from = measure, values_from = count) %>%
        # fill NA with 0
        mutate_all(~ replace(., is.na(.), 0)) %>%
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
            # label = paste0(round(value * 100, 0)),
            value = (value - MIN) / (MAX - MIN)
        ) %>%
        rename(x_axis = metric) %>%
        filter(x_axis %in% METRICS) %>%
        filter(!is.na(value))


    df_cross <- crossing(motif_mod = motifs_scored_f$motif_mod, bin = motifs_scored_f$bin) %>%
        filter(motif_mod %in% motifs) %>%
        left_join(motifs_scored_f, by = c("bin", "motif_mod")) %>%
        mutate(
            label = case_when(
                is.na(median) ~ "\u00B7",
                TRUE ~ NA
            ),
            median = ifelse(is.na(median), 0, median),
            group = "methylation"
        ) %>%
        rename(x_axis = motif_mod, value = median)


    df_cross %>%
        filter(x_axis %in% motifs) %>%
        bind_rows(metrics) %>%
        mutate(
            bin = factor(bin, levels = bin_rows),
            x_axis = factor(x_axis, levels = c(mod_cols, METRICS))
        ) %>%
        ggplot(aes(x = x_axis, y = fct_rev(bin), fill = value)) +
        geom_tile(color = "gray30") +
        # geom_text(data = metrics, aes(label = label)) +
        geom_text(data = df_cross, aes(label = label), size = 10) +
        scale_fill_gradientn(
            labels = scales::percent,
            limits = c(0, 1),
            colours = c("white", PLOT_COLORS[[4]], PLOT_COLORS[[5]]), na.value = "white",
            values = c(0, 0.5, 1)
        ) +
        theme_minimal() +
        labs(
            x = "",
            y = "",
            title = "Methylation Pattern and Benchmark Metrics",
        ) +
        facet_grid(~group, scales = "free", space = "free") +
        theme(
            axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5),
            # add bg to strio
            strip.background = element_rect(fill = "gray70"),
        )

    ggsave(file.path(args$output, paste0(s, "_motif_score.png")), width = 30, height = 12)
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

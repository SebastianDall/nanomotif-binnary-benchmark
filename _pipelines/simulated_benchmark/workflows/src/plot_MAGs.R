if (!require("argparse")) install.packages("argparse")
library(argparse)

# Create argument parser
parser <- ArgumentParser(description = "Plot MAGs script")

# Define arguments
parser$add_argument("--motifs_scored", required = TRUE, help = "Path to motifs-scored.tsv file")
parser$add_argument("--bin_motifs", required = TRUE, help = "Path to bin-motifs.tsv file")
parser$add_argument("--contig_bins", required = TRUE, help = "Path to contig_bins.tsv file")
parser$add_argument("--contig_bins_truth", required = TRUE, help = "Path to contig_bins_truth.tsv file")
parser$add_argument("--bin_contamination", required = TRUE, help = "Path to bin_contamination.tsv file")
parser$add_argument("--mean_methylation_cutoff", required = TRUE, help = "mean_methylation_cutoff parameter")
parser$add_argument("--n_motif_contig_cutoff", required = TRUE, help = "n_motif_contig_cutoff parameter")
parser$add_argument("--n_motif_bin_cutoff", required = TRUE, help = "n_motif_bin_cutoff parameter")
parser$add_argument("--output", required = TRUE, help = "Path to the output file")


# Parse arguments
args <- parser$parse_args()
args <- list(
    motifs_scored = "output/baseline/motif_discovery/original_no_shuffling/motifs-scored.tsv",
    bin_motifs = "output/baseline/motif_discovery/original_no_shuffling/bin-motifs.tsv",
    contig_bins = "files/benchmarks/benchmark_0_shuffle_1_contig.tsv", # "data/datasets/simulated_3_lognormal/contig_mapping/mapped_contig_bin.tsv",
    contig_bins_truth = "data/datasets/simulated_3_lognormal/contig_mapping/mapped_contig_bin.tsv",
    bin_contamination = "output/benchmarks/test/benchmark_0_shuffle_1_contig/bin_contamination.tsv",
    mean_methylation_cutoff = 0.25,
    n_motif_contig_cutoff = 10,
    n_motif_bin_cutoff = 500,
    output = "MAGs_warp"
)

# Load libraries
if (!require("data.table")) install.packages("data.table")
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("ggtext")) install.packages("ggtext")
if (!require("janitor")) install.packages("janitor")
library(data.table)
library(tidyverse)
library(ggtext)
library(janitor)


# Load data
motifs_scored <- fread(args$motifs_scored)

bin_motifs <- fread(args$bin_motifs) %>%
    mutate(
        motif_mod = paste0(motif, "_", mod_type, "-", mod_position),
        n_motifs = n_mod_bin + n_nomod_bin,
        mean = n_mod_bin / n_motifs
    ) %>%
    filter(
        n_motifs >= args$n_motif_bin_cutoff,
        mean >= args$mean_methylation_cutoff
    ) %>%
    pull(motif_mod) %>%
    unique()

contig_bins <- fread(args$contig_bins, header = FALSE) %>%
    rename(
        contig = 1,
        bin = 2
    )

contig_bins_truth <- fread(args$contig_bins_truth) %>% 
    rename(
        bin_truth = bin
    )

bin_contamination <- fread(args$bin_contamination) %>%
    select(bin, contig) %>%
    mutate(
        prediction = "contamination"
    )

## Data processing
motifs_scored <- motifs_scored %>%
    mutate(
        motif_mod = paste0(motif, "_", mod_type, "-", mod_position),
        n_motifs = n_mod + n_nomod,
        mean = n_mod / n_motifs
    ) %>%
    filter(
        motif_mod %in% bin_motifs,
        n_motifs >= args$n_motif_contig_cutoff
    ) %>%
    left_join(contig_bins, by = "contig") %>%
    left_join(contig_bins_truth, by = "contig") %>%
    filter(!is.na(bin)) %>%
    left_join(bin_contamination, by = c("bin", "contig")) %>%
    mutate(
        prediction = ifelse(is.na(prediction), "not contamination", prediction),
    ) %>%
    mutate(
        measure = case_when(
            prediction == "not contamination" & bin == bin_truth ~ "true negative",
            prediction == "not contamination" & bin != bin_truth ~ "false negative",
            prediction == "contamination" & bin == bin_truth ~ "false positive",
            prediction == "contamination" & bin != bin_truth ~ "true positive"
        )
    )



motifs_scored_for_plotting <- motifs_scored %>%
    mutate(
        contig_highlighting = case_when(
            measure == "true positive" ~ paste0("<span style='color:red'>", contig, " - ", measure, "</span>"),
            measure == "true negative" ~ paste0("<span style='color:green'>", contig, " - ", measure, "</span>"),
            measure == "false positive" ~ paste0("<span style='color:blue'>", contig, " - ", measure, "</span>"),
            measure == "false negative" ~ paste0("<span style='color:orange'>", contig, " - ", measure, "</span>")
        )
    )
# Plot each bin
for (BIN in motifs_scored$bin %>% unique()) {
    motifs <- motifs_scored_for_plotting %>%
        filter(bin == BIN) %>%
        group_by(motif_mod) %>%
        summarise(
            mean = mean(mean)
        ) %>%
        arrange(desc(mean))

    motifs_high <- motifs %>%
        filter(mean > 0.15) %>%
        pull(motif_mod)

    motifs_low <- motifs %>%
        filter(mean < 0.15) %>%
        pull(motif_mod)

    motifs_scored_for_plotting %>%
        filter(motif_mod %in% c(motifs_high, motifs_low[1:10])) %>%
        filter(bin == BIN) %>%
        ggplot(aes(x = motif_mod, y = contig_highlighting, fill = mean)) +
        geom_tile() +
        scale_fill_gradient2(limits = c(0, 1), low = "white", high = "darkblue", na.value = "white") +
        coord_equal() +
        theme(
            axis.text.y = element_markdown(vjust = 0.5, hjust = 0, size = 5), # Rotate x-axis labels to 90 degrees
            axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5, size = 5),
        ) +
        labs(
            title = paste0("Methylation pattern of ", BIN),
            x = "Motif",
            y = ""
        ) +
        theme(
            legend.position = "bottom"
        )

    ggsave(paste(args$output, "/", BIN, "_motifs.png", sep = ""), width = 10, height = 12)
}


# Plot confusion matrix
confusion_matrix <- motifs_scored %>%
    group_by(measure) %>%
    summarise(
        count = n()
    ) %>%
    right_join(
        tibble(measure = c("true positive", "true negative", "false positive", "false negative")),
        by = "measure"
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

metrics <- confusion_matrix %>%
    select(measure, count) %>%
    pivot_wider(names_from = measure, values_from = count) %>%
    clean_names() %>%
    mutate(
        sensitivity = true_positive / (true_positive + false_negative),
        specificity = true_negative / (true_negative + false_positive),
        accuracy = (true_positive + true_negative) / (true_positive + true_negative + false_positive + false_negative)
    )

confusion_matrix %>%
    ggplot(aes(x = actual, y = predicted, label = count_label)) +
    geom_tile(fill = "white", color = "black") +
    geom_text() +
    theme(
        legend.position = "bottom",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    ) +
    labs(
        title = "Confusion matrix",
        subtitle = paste0("Sensitivity: ", round(metrics$sensitivity, 2), " | Specificity: ", round(metrics$specificity, 2), " | Accuracy: ", round(metrics$accuracy, 2)),
        x = "Actual",
        y = "Predicted"
    )

ggsave(paste(args$output, "/confusion_matrix.png", sep = ""), width = 10, height = 10)

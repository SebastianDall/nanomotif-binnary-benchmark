library(tidyverse)


# Load data
contig_bin <- read_delim("data/datasets/simulated_4_lognormal_20-100x/contig_mapping/mapped_contig_bin.tsv", "\t", col_names = FALSE) %>%
    rename(
        contig = X1,
        bin = X2
    )
contig_mapping <- read_delim("data/datasets/simulated_4_lognormal_20-100x/contig_mapping/coverage.csv", ",")

motifs_scored <- read_delim("output/benchmarks/detect_contamination/2024-11-07_read-level/simulated_4_lognormal_20-100x/original_contig_bin/motifs-scored-read-methylation.tsv", "\t")

# load assembly
library(seqinr)
assembly <- read.fasta(file = "data/datasets/simulated_4_lognormal_20-100x/eukfilt_assembly.fasta", seqtype = "DNA")

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
cov <- read_delim("data/datasets/simulated_4_lognormal_20-100x/1_cov.bam_0_data_cov.txt", "\t") %>%
    rename(
        contig = 1,
        start = 2,
        length = 3,
        read_coverage_mean = 4
    )

motifs_scored_f <- motifs_scored %>%
    filter(N_motif_obs >= 10) %>%
    left_join(contig_bin) %>%
    left_join(assembly_info %>% select(contig, length)) %>%
    mutate(
        motif_mod = paste0(motif, "_", mod_type, "_", mod_position),
        contig_l = paste0(contig, " - ", length)
    ) %>%
    filter(!is.na(bin))


# join data


contig_bin_mapping <- contig_bin %>%
    left_join(contig_mapping) %>%
    left_join(assembly_info) %>%
    left_join(cov)

# assert all contig_mappings are unique
contig_mapping %>%
    group_by(contig) %>%
    filter(n() > 1) %>%
    arrange(contig)

contig_mapping %>% filter(coverage < 0.8)


contig_bin_mapping %>%
    group_by(contig) %>%
    filter(coverage == max(coverage)) %>%
    arrange(coverage)

contig_bin_mapping %>% filter(bin == -1)

# Plot contig mapping coverage
contig_bin_mapping %>%
    group_by(contig) %>%
    filter(coverage == max(coverage)) %>%
    ggplot(aes(x = coverage, y = bin, color = as.factor(bin))) +
    # dot plot
    geom_point(position = position_dodge(), size = 2)

# number of bins spanning multiple samples
x <- contig_bin_mapping %>%
    filter(bin != -1) %>%
    group_by(sample) %>%
    filter(n_distinct(bin) > 1) %>%
    arrange(sample, bin)



# GC / coverage plot
contig_bin_mapping %>%
    distinct(contig, .keep_all = TRUE) %>%
    # filter(bin != -1) %>%
    ggplot(aes(x = gc, y = read_coverage_mean, color = as.factor(sample), size = length)) +
    geom_point() +
    scale_y_log10() +
    scale_color_discrete() +
    theme_bw() +
    theme(
        legend.position = "none"
    )

ggsave("analysis/figures/simulated_4_lognormal_20-100x/gc_cov.png", width = 7, height = 7, create.dir = TRUE)



# sample / coverage plot
sample_mean_sorted <- contig_bin_mapping %>%
    group_by(sample) %>%
    summarise(
        mean_coverage = mean(read_coverage_mean),
    ) %>%
    arrange(mean_coverage)


contig_bin_mapping %>%
    distinct(contig, .keep_all = TRUE) %>%
    filter(!is.na(sample)) %>%
    mutate(
        sample = factor(sample, levels = sample_mean_sorted$sample)
    ) %>%
    ggplot(aes(x = read_coverage_mean, y = sample, color = as.factor(sample), size = length)) +
    scale_x_log10() +
    geom_point() +
    theme_bw() +
    theme(
        legend.position = "none"
    )

ggsave("analysis/figures/simulated_4_lognormal_20-100x/sample_coverage.png", width = 7, height = 10)


motifs_scored_f %>%
    distinct(contig, .keep_all = TRUE) %>%
    group_by(bin) %>%
    summarise(
        bin_l = sum(length)
    ) %>%
    arrange(bin_l)

motifs_scored_f %>%
    group_by(bin, motif_mod) %>%
    summarise(
        median = mean(median)
    ) %>%
    ggplot(aes(x = motif_mod, y = bin, fill = median)) +
    geom_tile(color = "gray30") +
    scale_fill_gradient2(low = "white", high = "blue") +
    theme_minimal() +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1)
    )
ggsave("analysis/figures/simulated_4_lognormal_20-100x/methylation_pattern.png", width = 10, height = 10)

library(tidyverse)


# Load data
contig_bin <- read_delim("data/datasets/simulated_3_lognormal/contig_bins.tsv", "\t")
contig_mapping <- read_delim("data/datasets/simulated_3_lognormal/contig_mapping/coverage.csv", ",")

# load assembly
library(seqinr)
assembly <- read.fasta(file = "data/datasets/simulated_3_lognormal/eukfilt_assembly.fasta", seqtype = "DNA")

assembly_info <- data.frame(
    contig = names(assembly),
    length = sapply(assembly, function(seq) {
        length(seq)
    }),
    gc = sapply(assembly, function(seq) {
        seqinr::GC(seq)
    })
)

cov <- read_delim("data/datasets/simulated_3_lognormal/1_cov.bam_0_data_cov.csv", ",") %>%
    rename(
        contig = 1,
        read_coverage_mean = 2,
        read_coverage_var = 3
    )

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

contig_bin_mapping %>%
    filter(bin == 6)



# GC / coverage plot
contig_bin_mapping %>%
    distinct(contig, .keep_all = TRUE) %>%
    # filter(bin != -1) %>%
    ggplot(aes(x = gc, y = read_coverage_mean, color = as.factor(sample), size = length)) +
    geom_point() +
    scale_y_log10() +
    scale_color_discrete()

# sample / coverage plot
contig_bin_mapping %>%
    distinct(contig, .keep_all = TRUE) %>%
    filter(!is.na(sample)) %>%
    ggplot(aes(x = read_coverage_mean, y = sample, color = as.factor(sample), size = length)) +
    scale_x_log10() +
    geom_point()

contig_mapping <- read_delim("data/datasets/simulated_3_lognormal/contig_mapping/coverage.csv", ",")

mapped_contig_bin <- contig_mapping %>%
    filter(coverage > 0.9) %>%
    group_by(contig) %>%
    filter(n() == 1) %>%
    rename(
        bin = sample
    ) %>%
    select(-coverage) %>%
    arrange(sample)

write_delim(mapped_contig_bin, "data/datasets/simulated_3_lognormal/mapped_contig_bin.tsv", delim = "\t", col_names = FALSE)

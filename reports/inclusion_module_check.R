library(tidyverse)
library(data.table)
library(ggtext)
source("src/colors.R")


samples <- c("fecal_simple_recovered", "ZymoFecal", "fecal_inhouse", "anaerobic_digester", "mfd02199_backup")
labels <- c("Fecal Simple", "ZymoFecal", "Fecal Complex", "Anaerobic Digester", "Soil")

df <- tibble(
    sample = samples,
    label = labels
) %>%
    mutate(
        include = map(sample, ~ fread(paste0("data/real_communities/", .x, "/nanomotif_binnary/include_contigs.tsv"))),
        contig_bin = map(sample, ~ fread(paste0("data/real_communities/", .x, "/nanomotif_binnary/new_contig_bin.tsv"))),
    )


df_include <- df %>%
    select(sample, label, include) %>%
    unnest(include) %>%
    filter(confidence == "high_confidence") %>%
    mutate(
        included = TRUE
    ) %>%
    select(sample, label, contig, included, assigned_bin) %>%
    distinct()

df_contig_bin <- df %>%
    select(sample, label, contig_bin) %>%
    unnest(contig_bin) %>%
    left_join(df_include %>% select(-assigned_bin)) %>%
    mutate(
        included = ifelse(is.na(included), FALSE, included)
    )

mags_of_interest <- df_include %>%
    mutate(sb = paste0(sample, "_", assigned_bin)) %>%
    pull(sb) %>%
    unique()

df_contig_bin_f <- df_contig_bin %>%
    mutate(
        sb = paste0(sample, "_", bin)
    ) %>%
    filter(sb %in% mags_of_interest) %>%
    nest(.by = c(sample, label))

for (s in df_contig_bin_f$sample) {
    motifs_scored <- read_delim(file.path("data/real_communities/", s, "nanomotif_binnary/motifs-scored-read-methylation.tsv"), "\t")
    contig_bin <- df_contig_bin_f %>%
        filter(sample == s) %>%
        unnest(data)

    bins <- contig_bin %>%
        pull(bin) %>%
        unique()

    for (b in bins) {
        contigs <- contig_bin %>%
            filter(bin == b) %>%
            pull(contig)
        included <- contig_bin %>%
            filter(bin == b, included) %>%
            pull(contig)

        motifs_scored_f <- motifs_scored %>%
            filter(contig %in% contigs) %>%
            filter((N_motif_obs * mean_read_cov) >= 24) %>%
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
                motif_mod = str_replace(motif_axis, "NNN*", paste0("(N)<sub>", ns, "</sub>")),
                contig_l = case_when(
                    contig %in% included ~ paste0("<span style='color:#404080'>", contig, "</span>"),
                    TRUE ~ contig
                )
            )

        motifs <- motifs_scored_f %>%
            filter(!is.na(median)) %>%
            # filter(median >= 0.1) %>%
            pull(motif_mod) %>%
            unique()

        df_wide <- motifs_scored_f %>%
            filter(motif_mod %in% motifs) %>%
            select(contig_l, motif_mod, median) %>%
            pivot_wider(names_from = motif_mod, values_from = median, values_fill = list(median = 0))

        # Convert to matrix
        data_matrix <- as.matrix(df_wide[, -1])
        rownames(data_matrix) <- df_wide$contig_l

        # Perform clustering
        dist_rows <- dist(data_matrix, method = "euclidean")
        dist_cols <- dist(t(data_matrix), method = "euclidean")
        hc_rows <- hclust(dist_rows, method = "average")
        hc_cols <- hclust(dist_cols, method = "average")

        mod_cols <- colnames(data_matrix)[hc_cols$order]
        bin_rows <- rownames(data_matrix)[hc_rows$order]

        df_cross <- crossing(
            motif_mod = motifs,
            contig_l = bin_rows
        ) %>%
            left_join(motifs_scored_f) %>%
            mutate(
                contigl_l = factor(contig_l, levels = bin_rows),
                motif_mod = factor(motif_mod, levels = c(mod_cols)),
                txt = case_when(
                    is.na(median) ~ "\u00B7",
                    TRUE ~ NA
                ),
                median = ifelse(is.na(median), 0, median),
            )

        df_cross %>%
            ggplot(aes(x = motif_mod, y = contig_l, fill = median)) +
            geom_tile(color = "gray30") +
            geom_text(aes(label = txt), size = 10) +
            scale_fill_gradientn(
                labels = scales::percent,
                limits = c(0, 1),
                colours = c("white", PLOT_COLORS[[4]], PLOT_COLORS[[5]]), na.value = "white",
                values = c(0, 0.5, 1)
            ) +
            theme_minimal() +
            labs(
                x = "",
                y = ""
            ) +
            theme(
                axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5),
                axis.text.y = element_markdown(),
                # add bg to strio
                strip.background = element_rect(fill = "gray70"),
            )

        ggsave(file.path("analysis/include_story/included_contigs_all/", paste0(s, "_", b, "_motif_score.png")), width = 30, height = 12)
    }
}



include <- fread("data/real_communities/mfd02199_backup/nanomotif_binnary/include_contigs.tsv")


include %>%
    filter(confidence == "high_confidence") %>%
    distinct(contig, assigned_bin) %>%
    group_by(assigned_bin) %>%
    summarise(
        n = n()
    ) %>%
    arrange(desc(n)) %>%
    summarise(
        median = median(n)
    )

include %>%
    filter(confidence == "high_confidence") %>%
    distinct(contig, assigned_bin) %>%
    group_by(assigned_bin) %>%
    summarise(
        n = n()
    ) %>%
    ggplot(aes(x = n)) +
    geom_histogram(binwidth = 10) +
    theme_minimal() +
    geom_vline(xintercept = 100, color = "red") +
    annotate("text", x = 125, y = 10, label = "100", color = "red")

ggsave("analysis/include_check/include_contigs_hist.png", width = 6, height = 5, dpi = 300)

motifs_scored <- fread("data/real_communities/mfd02199_backup/nanomotif_binnary/motifs-scored-read-methylation.tsv")

include <- fread("data/real_communities/mfd02199_backup/nanomotif_binnary/include_contigs.tsv") %>%
    mutate(
        included = TRUE
    )

# assembly <- read.fasta("data/real_communities/anaerobic_digester/mmlong2_lite/results/mmlong2_assembly.fasta", seqtype = "DNA")

# assembly_info <- data.frame(
#     contig = names(assembly),
#     length = sapply(assembly, function(seq) {
#         length(seq)
#     }),
#     gc = sapply(assembly, function(seq) {
#         seqinr::GC(seq)
#     })
# )
# rm(assembly)
# assembly_info <- fread("reports/anaerobic_digester_assembly_info.txt")
assembly_info <- fread("data/real_communities/mfd02199_backup/mmlong2_lite/tmp/flye/assembly_info.txt") %>%
    rename(
        contig = 1,
        length = 2,
        cov = 3
    )
cov <- fread("data/real_communities/mfd02199_backup/mmlong2_lite/tmp/binning/cov_tmp/1_cov.tsv") %>% # "data/real_communities/mfd02199_backup/mmlong2_lite/tmp/binning/mapping/cov_all.tsv"
    rename(
        contig = 1,
        length = 2,
        read_coverage_mean = 3
    ) %>%
    select(contig, length, read_coverage_mean)

checkm <- fread("data/real_communities/mfd02199_backup/mmlong2_lite/results/mmlong2_lite_bins.tsv") %>%
    rename(
        bin = 1,
        completeness = 2,
        contamination = 3
    ) %>%
    select(bin, completeness, contamination) %>%
    mutate(
        mag_qual = case_when(
            completeness >= 90 & contamination < 5 ~ "HQ",
            completeness >= 50 & contamination < 10 ~ "MQ",
            TRUE ~ "low_quality"
        )
    )

contig_bin <- fread("data/real_communities/mfd02199_backup/nanomotif_binnary/new_contig_bin.tsv") %>%
    mutate(
        included = case_when(
            contig %in% filter(include, confidence == "high_confidence")$contig ~ TRUE,
            TRUE ~ FALSE
        )
    ) %>%
    distinct(contig, bin, included) # %>%
# filter(bin %in% filter(checkm, mag_qual == "HQ")$bin)

# all contig in contig_bin are unique()
length(contig_bin$contig) == length(unique(contig_bin$contig))

genomad <- fread("data/real_communities/mfd02199_backup/genomad/mmlong2_lite_assembly_aggregated_classification/mmlong2_lite_assembly_aggregated_classification.tsv") %>%
    rename(
        contig = 1
    ) %>%
    pivot_longer(cols = -contig, names_to = "mge", values_to = "count") %>%
    group_by(contig) %>%
    filter(count == max(count)) %>%
    select(-count) %>%
    mutate(
        mge = str_remove(mge, "_score")
    )

BIN <- "mmlong2_lite.bin.3.1026"
contigs_in_bin <- contig_bin %>%
    filter(bin == BIN) %>%
    pull(contig)

binstats <- checkm %>%
    filter(bin == BIN) %>%
    mutate(
        bin = str_remove(bin, "mmlong2_lite.")
    )

motifs_scored_f <- motifs_scored %>%
    filter(contig %in% contigs_in_bin) %>%
    filter(N_motif_obs * mean_read_cov >= 24) %>%
    left_join(contig_bin) %>%
    left_join(assembly_info %>% select(contig, length)) %>%
    left_join(genomad) %>%
    left_join(checkm) %>%
    mutate(
        motif_mod = paste0(motif, "_", mod_type, "_", mod_position),
        contig_l = paste0(contig, " - ", length, " - ", mge),
        bin = str_remove(bin, "mmlong2_lite."),
        bin_label = paste0(bin, " - ", mag_qual, " - ", completeness, " - ", contamination)
    ) %>%
    mutate(bin = ifelse(is.na(bin), "unbinned", bin)) %>%
    mutate(
        included = ifelse(is.na(included), FALSE, included)
    )


gc_cov <- assembly_info %>%
    left_join(cov) %>%
    left_join(contig_bin) %>%
    left_join(genomad) %>%
    mutate(
        bin = ifelse(is.na(bin), "unbinned", bin)
    ) %>%
    mutate(
        included = ifelse(is.na(included), FALSE, included)
    ) %>%
    mutate(
        bin = str_remove(bin, "mmlong2_lite.")
    )

BIN <- str_remove(BIN, "mmlong2_lite.")

bin_gc_cov <- gc_cov %>%
    filter(bin == BIN) %>%
    arrange(desc(read_coverage_mean))

contigs_of_interest <- bin_gc_cov %>%
    filter(included) %>%
    pull(contig)

contig_gc_cov <- gc_cov %>%
    filter(contig %in% contigs_of_interest)

gc_cov %>%
    filter(bin != BIN) %>%
    filter(read_coverage_mean > 5) %>%
    filter(length > 10000) %>%
    ggplot(aes(x = gc, y = read_coverage_mean, shape = mge, size = length)) +
    geom_point(color = "gray50", alpha = 0.6) +
    geom_point(data = bin_gc_cov, alpha = 0.8, color = PLOT_COLORS[[2]]) +
    geom_point(data = contig_gc_cov, alpha = 1, aes(color = mge)) +
    scale_color_manual(values = c(PLOT_COLORS[[4]], "#e9c600", PLOT_COLORS[[2]]), guide = "none") +
    guides(
        shape = guide_legend(title = "MGE"),
        size = guide_legend(title = "Contig length")
    ) +
    theme_minimal() +
    labs(
        x = "GC content",
        y = "Coverage",
        color = "Included",
        title = binstats$bin[1],
        subtitle = paste0(binstats$mag_qual[1], " - Comp. ", binstats$completeness[1], " - Cont. ", binstats$contamination[1])
    ) +
    scale_y_log10() +
    theme(
        # square plot
        aspect.ratio = 1,
    )

ggsave(paste0("analysis/include_check/", binstats$mag_qual[1], "_", binstats$bin[1], "_gc_cov.png"), width = 6, height = 5, dpi = 300)


motifs <- motifs_scored_f %>%
    filter(bin == BIN) %>%
    filter(median > 0.3) %>%
    pull(motif_mod) %>%
    unique()

p <- motifs_scored_f %>%
    filter(bin == BIN) %>%
    filter(motif_mod %in% motifs) %>%
    mutate(
        contig_l = case_when(
            contig %in% contigs_of_interest & mge == "plasmid" ~ paste0("<span style='color:", PLOT_COLORS[[4]], "'>", contig_l, "</span>"),
            contig %in% contigs_of_interest & mge == "virus" ~ paste0("<span style='color:", "#e9c600", "'>", contig_l, "</span>"),
            TRUE ~ contig_l
        )
    ) %>%
    select(median, contig_l, motif, mod_type, mod_position, included) %>%
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
        motif_axis = str_replace(motif_axis, "NNN*", paste0("(N)<sub>", ns, "</sub>")),
    ) %>%
    select(contig_l, motif_axis, median, included) %>%
    pivot_wider(names_from = motif_axis, values_from = median) %>%
    pivot_longer(cols = -c(contig_l, included), names_to = "motif_axis", values_to = "median") %>%
    mutate(
        txt = case_when(
            is.na(median) ~ "\u00B7",
            TRUE ~ NA
        ),
        median = ifelse(is.na(median), 0, median)
    ) %>%
    ggplot(aes(x = motif_axis, y = contig_l, fill = median, label = txt)) +
    geom_tile(color = "gray50") +
    geom_text(aes(label = txt), size = 15) +
    scale_fill_gradientn(
        labels = scales::percent,
        limits = c(0, 1),
        colours = c("white", PLOT_COLORS[[4]], PLOT_COLORS[[5]]), na.value = "white",
        values = c(0, 0.5, 1)
    ) +
    facet_grid(included ~ ., scales = "free", space = "free") +
    theme_minimal() +
    labs(
        x = "",
        y = "",
        fill = "Methylation degree",
        title = binstats$bin[1],
        subtitle = paste0(binstats$mag_qual[1], " - Comp. ", binstats$completeness[1], " - Cont. ", binstats$contamination[1])
    ) +
    theme(
        axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_markdown(),
        legend.position = "none"
    )
p
ggsave(paste0("analysis/include_check/soil_", binstats$mag_qual[1], "_", binstats$bin[1], "_motifs.png"), width = 30, height = 30, dpi = 300, create.dir = TRUE, plot = p)

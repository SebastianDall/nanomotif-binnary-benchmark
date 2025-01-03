library(tidyverse)
library(data.table)
library(seqinr)
library(ggtext)
source("src/colors.R")

motifs_scored <- fread("data/real_communities/fecal_simple_recovered/nanomotif_binnary/motifs-scored-read-methylation.tsv")

include <- fread("data/real_communities/fecal_simple_recovered/nanomotif_binnary/include_contigs.tsv") %>%
    mutate(
        included = TRUE
    )

assembly <- read.fasta("data/real_communities/fecal_simple_recovered/mmlong2_lite/results/mmlong2_assembly.fasta", seqtype = "DNA")

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

cov <- fread("data/real_communities/fecal_simple_recovered/mmlong2_lite/tmp/binning/mapping/cov_all.tsv") %>%
    rename(
        contig = 1,
        length = 2,
        read_coverage_mean = 3
    ) %>%
    select(contig, length, read_coverage_mean)

checkm <- fread("data/real_communities/fecal_simple_recovered/mmlong2_lite/results/mmlong2_bins.tsv") %>%
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

contig_bin <- fread("data/real_communities/fecal_simple_recovered/nanomotif_binnary/new_contig_bin.tsv") %>%
    mutate(
        included = case_when(
            contig %in% filter(include, confidence == "high_confidence")$contig ~ TRUE,
            TRUE ~ FALSE
        )
    ) %>%
    distinct(contig, bin, included) %>%
    filter(bin %in% filter(checkm, mag_qual == "HQ")$bin)

# all contig in contig_bin are unique()
length(contig_bin$contig) == length(unique(contig_bin$contig))

genomad <- fread("data/real_communities/fecal_simple_recovered/genomad/mmlong2_assembly_aggregated_classification/mmlong2_assembly_aggregated_classification.tsv") %>%
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

motifs_scored_f <- motifs_scored %>%
    filter(N_motif_obs * mean_read_cov >= 24) %>%
    left_join(contig_bin) %>%
    left_join(assembly_info %>% select(contig, length)) %>%
    left_join(genomad) %>%
    left_join(checkm) %>%
    mutate(
        motif_mod = paste0(motif, "_", mod_type, "_", mod_position),
        contig_l = paste0(contig, " - ", round(length / 1000, 0), " kbp - ", mge),
        bin = str_remove(bin, "mmlong2."),
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
        bin = str_remove(bin, "mmlong2.")
    )

gc_cov$bin %>% unique()

BIN <- "bin.1.7"

bin_gc_cov <- gc_cov %>%
    filter(bin == BIN) %>%
    arrange(desc(read_coverage_mean))

# contigs_of_interest <- bin_gc_cov %>%
#     filter(included) %>%
#     pull(contig)

contigs_of_interest <- c("contig_11", "contig_318", "contig_1858")

contig_gc_cov <- gc_cov %>%
    filter(contig %in% contigs_of_interest)

gc_cov %>%
    filter(bin != BIN) %>%
    ggplot(aes(x = gc, y = read_coverage_mean, shape = mge, size = length)) +
    geom_point(color = "gray50", alpha = 0.6) +
    geom_point(data = bin_gc_cov, alpha = 0.8, color = PLOT_COLORS[[2]]) +
    geom_point(data = contig_gc_cov, alpha = 1, aes(color = mge)) +
    scale_color_manual(values = c(PLOT_COLORS[[4]], "#e9c600"), guide = "none") +
    guides(
        shape = guide_legend(title = "MGE"),
        size = guide_legend(title = "Contig length")
    ) +
    theme_minimal() +
    labs(
        x = "GC content",
        y = "Coverage",
        color = "Included",
        title = "Simple Fecal",
        subtitle = "Contigs added to bin.1.7"
    ) +
    scale_y_log10() +
    theme(
        # square plot
        aspect.ratio = 1,
    )

ggsave("analysis/include_story/GC_cov.png", width = 6, height = 5, dpi = 300)


motifs <- c(
    "AAGRAG_a_4",
    "CAACAT_a_4",
    "CATANNNNNNTTC_a_3",
    "CATATG_a_3",
    "CGGGAG_a_4",
    "CGWAAT_a_4",
    "GACGTC_a_1",
    "GATC_a_1",
    "GCATC_a_2",
    "GGANNNNNNTATC_a_2",
    "RGATCY_a_2",
    "TTCGAA_a_5",
    "CCATGG_m_0",
    "CCWGG_m_1",
    "GCNGC_m_1",
    "GATC_m_3"
)

bin_included <- motifs_scored_f %>%
    filter(included, bin == BIN) %>%
    filter(contig %in% contigs_of_interest) %>%
    filter(motif_mod %in% motifs) %>%
    mutate(
        contig_l = case_when(
            contig %in% contigs_of_interest & mge == "plasmid" ~ paste0("<span style='color:", PLOT_COLORS[[4]], "'>", contig_l, "</span>"),
            contig %in% contigs_of_interest & mge == "virus" ~ paste0("<span style='color:", "#e9c600", "'>", contig_l, "</span>"),
            TRUE ~ contig_l
        )
    ) %>%
    select(-bin_label, median, contig_l, motif_mod, included) %>%
    rename(
        bin_label = contig_l,
        mean = median
    ) %>%
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
    select(bin_label, motif_axis, mean, included) %>%
    pivot_wider(names_from = motif_axis, values_from = mean) %>%
    pivot_longer(cols = -c(bin_label, included), names_to = "motif_axis", values_to = "mean") %>%
    mutate(
        txt = case_when(
            is.na(mean) ~ "\u00B7",
            TRUE ~ NA
        ),
        mean = ifelse(is.na(mean), 0, mean)
    )

bin_scores <- motifs_scored_f %>%
    filter(!included, bin != "unbinned")

bins_fct <- bin_scores %>%
    pull(bin_label) %>%
    unique()
# make bin 1.7 the first
bins_fct <- c(bins_fct[bins_fct != "bin.1.7 - HQ - 99.86 - 0.2"], "bin.1.7 - HQ - 99.86 - 0.2", bin_included$bin_label %>% unique())

bin_scores %>%
    filter(motif_mod %in% motifs) %>%
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
        motif_axis = str_replace(motif_axis, "NNN*", paste0("(N)<sub>", ns, "</sub>"))
    ) %>%
    group_by(bin_label, motif_axis, included) %>%
    summarise(
        mean = mean(median),
    ) %>%
    mutate(
        txt = case_when(
            is.na(mean) ~ "\u00B7",
            TRUE ~ NA
        ),
        mean = ifelse(is.na(mean), 0, mean),
        bin_label = case_when(
            bin_label == "bin.1.7 - HQ - 99.86 - 0.2" ~ paste0("<strong><span style='color:", PLOT_COLORS[[2]], "'>", bin_label, "</span></strong>"),
            TRUE ~ bin_label
        )
    ) %>%
    bind_rows(bin_included) %>%
    mutate(
        included = ifelse(included, paste0("Added to ", BIN), "Bin") # ,
        # bin_label = factor(bin_label, levels = bins_fct, labels = bins_fct)
    ) %>%
    ggplot(aes(x = motif_axis, y = fct_rev(bin_label), fill = mean, label = txt)) +
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
        title = "Simple Fecal",
        subtitle = "Contigs added to bin.1.7"
    ) +
    theme(
        axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_markdown(),
        # remove strips
        strip.text = element_blank(),
        strip.background = element_blank(),
        legend.position = "none"
    )

ggsave("analysis/include_story/bin_scores.png", width = 6, height = 5, dpi = 300)
